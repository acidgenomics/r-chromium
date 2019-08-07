#' Import counts from either HDF5 or MTX files.
#' @note Updated 2019-08-01.
#' @noRd
.importCounts <- function(
    sampleDirs,
    filtered,
    BPPARAM = BiocParallel::SerialParam()
) {
    assert(
        allAreDirectories(sampleDirs),
        hasValidNames(sampleDirs),
        isFlag(filtered)
    )
    matrixFiles <- unlist(bplapply(
        X = sampleDirs,
        FUN = .findCountMatrix,
        filtered = filtered,
        BPPARAM = BPPARAM
    ))
    assert(
        hasLength(unique(basename(matrixFiles)), n = 1L),
        hasValidNames(matrixFiles)
    )
    if (all(grepl("\\.H5", matrixFiles, ignore.case = TRUE))) {
        fun <- .importCountsFromHDF5
    } else if (all(grepl("\\.MTX", matrixFiles, ignore.case = TRUE))) {
        fun <- .importCountsFromMTX
    } else {
        stop(sprintf("Unexpected import failure.", file))
    }
    message(sprintf(
        "Importing counts from '%s' file.",
        basename(matrixFiles[[1L]])
    ))
    ## This step seems to have issues when parsing HDF5 files in parallel
    ## on an Azure Files mount over CIFS.
    list <- bpmapply(
        sampleID = names(matrixFiles),
        file = matrixFile,
        FUN = function(sampleID, file) {
            counts <- fun(file)
            ## Strip index when all barcodes end with "-1".
            if (all(grepl("-1$", colnames(counts)))) {
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
                grepl("^([[:digit:]]+)-([ACGT]+)$", colnames(counts))
            ) {
                colnames(counts) <- paste(sampleID, colnames(counts), sep = "-")
            }
            ## Ensure names are valid. Don't sanitize the gene identifier rows,
            ## which can contain some invalid names due to gene symbols. This is
            ## an edge case that happens with non-standard genomes and/or
            ## spike-in names.
            colnames(counts) <- makeNames(colnames(counts))
            counts
        },
        SIMPLIFY = FALSE,
        USE.NAMES = TRUE,
        BPPARAM = BPPARAM
    )
    # Bind the matrices.
    do.call(what = cbind, args = list)
}





#' Import Cell Ranger count matrix from HDF5 file
#' 
#' @note Updated 2019-07-31.
#' @noRd
#' 
#' @seealso `cellrangerRkit::get_matrix_from_h5()`
#' 
#' @return `sparseMatrix`.
#'   Cell barcodes in the columns, features (i.e. genes) in the rows.
#'   
#' @examples
#' ## > x <- importCountsHDF5(file = "filtered_feature_bc_matrix.h5")
#' ## > dim(x)
.importCountsFromHDF5 <-  # nolint
    function(file) {
        assert(
            isAFile(file),
            grepl(pattern = "\\.H5", x = file, ignore.case = TRUE)
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
        
        ## Want `Csparse` not `Tsparse` matrix.
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
        
        rownames(counts) <- rownames
        colnames(counts) <- colnames
        
        counts
    }



#' Import Cell Ranger count matrix from MTX file
#' 
#' @note Data import using HDF5 file is now recommended over this approach.
#' @note Updated 2019-08-01.
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
#' ## > x <- importCountsMTX(file = "matrix.mtx.gz")
#' ## > dim(x)
.importCountsFromMTX <-  # nolint
    function(file) {
        assert(
            isAFile(file),
            grepl(pattern = "\\.MTX", x = file, ignore.case = TRUE)
        )
        path <- dirname(realpath(file))
        files <- list.files(
            path = path,
            pattern = "*.tsv*",
            full.names = FALSE,
            recursive = FALSE,
            ignore.case = TRUE
        )
        
        counts <- readMM(file)
        assert(is(counts, "sparseMatrix"))
        
        ## Row names.
        ## Legacy v2 output uses "genes" instead of "features".
        ## Also, older Cell Ranger output doesn't compress these files.
        rownamesFile <- grep(
            pattern = "^(features|genes)\\.tsv(\\.gz)?$",
            x = files,
            value = TRUE
        )
        assert(isString(rownamesFile))
        ## Note that `features.tsv` is tab delimited.
        ## v2: id, name
        ## v3: id, name, expression type
        rownames <- read_tsv(
            file = rownamesFile,
            col_names = FALSE
        )
        rownames <- rownames[[1L]]
        
        ## Column names.
        colnamesFile <- grep(
            pattern = "^barcodes\\.tsv(\\.gz)?$",
            x = files,
            value = TRUE
        )
        assert(isString(colnamesFiles))
        colnamesFile <- file.path(path, colnamesFile)
        ## Note that `barcodes.tsv` is NOT tab delimited.
        colnames <- read_lines(barcodesFile)
        
        assert(
            identical(length(rownames), nrow(counts)),
            identical(length(colnames), ncol(counts))
        )
        
        rownames(counts) <- rownames
        colnames(counts) <- colnames
        
        counts
    }
