#' Import counts from either HDF5 or MTX files.
#'
#' @note Updated 2022-06-07.
#' @noRd
.importCounts <-
    function(matrixFiles) {
        assert(allAreFiles(matrixFiles))
        if (all(grepl("\\.h5", matrixFiles, ignore.case = TRUE))) {
            fun <- .importCountsFromHdf5
        } else if (all(grepl("\\.mtx", matrixFiles, ignore.case = TRUE))) {
            fun <- .importCountsFromMtx
        } else {
            abort("Unexpected import failure.") # nocov
        }
        alert(sprintf(
            fmt = "Importing counts from {.file %s} %s.",
            basename(matrixFiles[[1L]]),
            ngettext(n = length(matrixFiles), msg1 = "file", msg2 = "files")
        ))
        list <- Map(
            sampleId = names(matrixFiles),
            file = matrixFiles,
            f = function(sampleId, file) {
                counts <- fun(file)
                ## Strip index when all barcodes end with "-1".
                if (
                    allAreMatchingRegex(x = colnames(counts), pattern = "-1$")
                ) {
                    colnames(counts) <- sub("-1", "", colnames(counts))
                }
                # Now move the multiplexed index name/number to the beginning,
                # for more logical sorting and consistency with bcbioSingleCell.
                colnames(counts) <- sub(
                    pattern = "^([ACGT]+)-(.+)$",
                    replacement = "\\2-\\1",
                    x = colnames(counts)
                )
                # Prefix cell barcodes with sample identifier when we're loading
                # counts from multiple samples.
                if (
                    length(matrixFiles) > 1L ||
                        allAreMatchingRegex(
                            x = colnames(counts),
                            pattern = "^([[:digit:]]+)-([ACGT]+)$"
                        )
                ) {
                    colnames(counts) <-
                        paste(sampleId, colnames(counts), sep = "-")
                }
                ## Ensure dimnames are valid. Note that this may sanitize gene
                ## names for some model systems, and or transgenes.
                counts <- makeDimnames(counts)
                counts
            }
        )
        # Bind the matrices.
        do.call(what = cbind, args = list)
    }



#' Import Cell Ranger count matrix from HDF5 file
#'
#' @note Updated 2019-08-21.
#' @noRd
#'
#' @seealso `cellrangerRkit::get_matrix_from_h5()`
#'
#' @return `sparseMatrix`.
#' Cell barcodes in the columns, features (i.e. genes) in the rows.
#'
#' @examples
#' ## > x <- .importCountsFromHdf5(file = "filtered_feature_bc_matrix.h5")
#' ## > dim(x)
.importCountsFromHdf5 <- # nolint
    function(file) {
        assert(
            isAFile(file),
            grepl(pattern = "\\.h5", x = file, ignore.case = TRUE)
        )
        names <- names(h5dump(file, load = FALSE))
        assert(isString(names))
        ## Import HDF5 data.
        h5 <- h5read(file = file, name = names[[1L]])
        ## v3 names:
        ## - "barcodes"
        ## - "data"
        ## - "features"
        ## - "indices"
        ## - "indptr"
        counts <- sparseMatrix(
            i = h5[["indices"]] + 1L,
            p = h5[["indptr"]],
            x = as.numeric(h5[["data"]]),
            dims = h5[["shape"]],
            giveCsparse = TRUE
        )
        ## Row names.
        if ("features" %in% names(h5)) {
            ## > names(h5[["features"]])
            ## [1] "_all_tag_keys" "feature_type"
            ## [3] "genome"        "id"
            ## [5] "name"          "pattern"
            ## [7] "read"          "sequence"
            ## Stable gene identifiers are stored in "id".
            ## Gene symbols are stored in "name".
            rownames <- h5[["features"]][["id"]]
        } else if ("genes" %in% names(h5)) {
            ## Older H5 objects (v2) use "genes" instead of "features".
            rownames <- h5[["genes"]]
        }
        ## Column names.
        colnames <- h5[["barcodes"]]
        assert(
            identical(length(rownames), nrow(counts)),
            identical(length(colnames), ncol(counts))
        )
        ## Return.
        rownames(counts) <- as.character(rownames)
        colnames(counts) <- as.character(colnames)
        counts
    }



#' Import Cell Ranger count matrix from MTX file
#'
#' @note Data import using HDF5 file is now recommended over this approach.
#' @note Updated 2022-06-07.
#' @noRd
#'
#' @section Matrix Market Exchange (MEX/MTX) format:
#'
#' Loading from this matrix requires sidecar files containing cell barcodes and
#' feature (i.e. gene) identifiers.
#'
#' Cell Ranger v3:
#'
#' - `barcodes.tsv.gz`: Cell barcodes.
#' - `features.tsv.gz`: Feature identifiers.
#'
#' Cell Ranger v2:
#'
#' - `barcodes.tsv`: Cell barcodes.
#' - `genes.tsv`: Gene identifiers.
#'
#' @examples
#' ## > x <- .importCountsFromMtx(file = "matrix.mtx.gz")
#' ## > dim(x)
.importCountsFromMtx <-
    function(file) {
        assert(
            isAFile(file),
            grepl(pattern = "\\.mtx", x = file, ignore.case = TRUE)
        )
        path <- dirname(realpath(file))
        files <- sort(list.files(
            path = path,
            pattern = "*.tsv*",
            full.names = FALSE,
            recursive = FALSE,
            ignore.case = TRUE
        ))
        ## Count matrix.
        counts <- import(file)
        assert(is(counts, "sparseMatrix"))
        ## Row names.
        ## Legacy v2 output uses "genes" instead of "features".
        ## Also, older Cell Ranger output doesn't compress these files.
        rownamesFile <- grep(
            pattern = "^(features|genes)\\.tsv(\\.gz)?$",
            x = files,
            value = TRUE,
            ignore.case = TRUE
        )
        rownamesFile <- file.path(path, rownamesFile)
        assert(isAFile(rownamesFile))
        ## Note that `features.tsv` is tab delimited.
        ## v2: id, name
        ## v3: id, name, expression type
        rownames <- import(
            con = rownamesFile,
            format = "tsv",
            colnames = FALSE
        )
        rownames <- rownames[[1L]]
        ## Column names.
        colnamesFile <- grep(
            pattern = "^barcodes\\.tsv(\\.gz)?$",
            x = files,
            value = TRUE,
            ignore.case = TRUE
        )
        colnamesFile <- file.path(path, colnamesFile)
        assert(isAFile(colnamesFile))
        ## Note that `barcodes.tsv` is source code lines, not tab delimited.
        colnames <- import(
            con = colnamesFile,
            format = "tsv",
            colnames = FALSE
        )
        colnames <- colnames[[1L]]
        ## Return.
        assert(
            identical(length(rownames), nrow(counts)),
            identical(length(colnames), ncol(counts))
        )
        rownames(counts) <- as.character(rownames)
        colnames(counts) <- as.character(colnames)
        counts
    }



#' Import sample-level metrics
#'
#' @note Updated 2023-09-28.
#' @noRd
.importSampleMetrics <-
    function(sampleDirs) {
        files <- file.path(sampleDirs, "outs", "metrics_summary.csv")
        assert(allAreFiles(files))
        list <- lapply(X = files, FUN = import)
        out <- do.call(what = rbind, args = list)
        ## Ensure we sanitize marked UTF-8 strings.
        out <- lapply(
            X = out,
            FUN = gsub,
            pattern = "[^.[:alnum:]]",
            replacement = ""
        )
        out <- do.call(what = DataFrame, args = out)
        out <- camelCase(out, strict = TRUE)
        rownames(out) <- names(sampleDirs)
        out
    }
