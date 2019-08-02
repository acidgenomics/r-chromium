#' Import Cell Ranger count matrix from MTX file
#' 
#' @note Data import using HDF5 file is now recommended over this approach.
#' @note Updated 2019-08-01.
#' @export
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
#' @inheritParams basejump::params
#' 
#' @examples
#' ## > x <- importCountsMTX(file = "matrix.mtx.gz")
#' ## > dim(x)
importCountsMTX <-  # nolint
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
