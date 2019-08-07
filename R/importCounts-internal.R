#' Find the count matrix inside a Cell Ranger sample directory.
#' 
#' This function will automatically move inside an `outs/` directory inside
#' the path if it is detected.
#' 
#' Currently preferring HDF5 over MTX.
#' 
#' @section File name history:
#' 
#' Cell Ranger v3
#' 
#' - H5: `filtered_feature_bc_matrix.h5`.
#' - MTX: `filtered_feature_bc_matrix/matrix.mtx.gz`.
#' 
#' Cell Ranger v2
#' 
#' - H5: `filtered_gene_bc_matrices_h5.h5`
#' - MTX: `filtered_gene_bc_matrices/<genomeBuild>/matrix.mtx`
#' 
#' @param dir Sample directory.
#' @param filtered `logical(1)`.
#'   - `TRUE`: Look for `filtered_*` matrix.
#'   - `FALSE`: Look for `raw_*` matrix.
#'   
#'   Doesn't apply if there's only a single matrix file in the directory.
#' 
#' @note Updated 2019-08-07.
#' @noRd
.findCountMatrix <- function(dir, filtered = TRUE) {
    assert(
        isADirectory(dir),
        isFlag(filtered)
    )
    dir <- realpath(dir)
    
    ## Simple mode -------------------------------------------------------------
    ## For minimal examples and data downloaded from 10X website.
    if (!dir.exists(file.path(dir, "outs"))) {
        file <- list.files(
            path = dir,
            pattern = "matrix\\.(h5|mtx)(\\.gz)?",
            full.names = TRUE
        )
        if (isAFile(file)) return(file)
    }
    
    ## Standard Cell Ranger output ---------------------------------------------
    ## Recurse into `outs/` directory by default.
    dir <- file.path(dir, "outs")
    assert(isADirectory(dir))
    
    if (isTRUE(filtered)) {
        prefix <- "filtered"
    } else {
        prefix <- "raw"
    }
    
    files <- list.files(
        path = dir,
        pattern = paste0("^", prefix, "_"),
        recursive = FALSE,
        full.names = FALSE
    )
    assert(hasLength(files))
    
    ## Get the Cell Ranger version, based on the file names.
    if (isTRUE(any(grepl(
        pattern = paste0("^", prefix, "_feature_bc_matrix$"),
        x = files
    )))) {
        version <- "3"
        filestem <- paste0(prefix, "_feature_bc_matrix")
    } else if (isTRUE(any(grepl(
        pattern = paste0("^", prefix, "_gene_bc_matrices$"),
        x = files
    )))) {
        version <- "2"
        filestem <- paste0(prefix, "_gene_bc_matrices")
    } else {
        stop("Failed to detect Cell Ranger version based on file names.")
    }
    version <- numeric_version(version)
    
    ## Currently preferring HDF5 over MTX.
    if (isTRUE(
        file.exists(file.path(dir, paste0(filestem, ".h5")))
    )) {
        ## v3 HDF5
        file <- file.path(dir, paste0(filestem, ".h5"))
    } else if (isTRUE(
        file.exists(file.path(dir, paste0(filestem, "_h5.h5")))
    )) {
        ## v2 HDF5
        file <- file.path(dir, paste0(filestem, "_h5.h5"))
    } else if (isTRUE(
        file.exists(file.path(dir, filestem, "matrix.mtx.gz"))
    )) {
        ## v3 MTX
        file <- file.path(dir, filestem, "matrix.mtx.gz")
    } else if (isTRUE(
        dir.exists(file.path(dir, filestem))
    )) {
        ## v2 MTX
        ## Get the genome build from the first sample directory.
        genomeBuild <- list.dirs(
            path = dir[[1L]],
            full.names = FALSE,
            recursive = FALSE
        )
        assert(isString(genomeBuild))
        file <- file.path(dir, genomeBuild, "matrix.mtx")
    }
    
    assert(isAFile(file))
    file
}



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
        fmt = "Importing counts from '%s' %s.",
        basename(matrixFiles[[1L]]),
        ngettext(n = length(matrixFiles), msg1 = "file", msg2 = "files")
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
