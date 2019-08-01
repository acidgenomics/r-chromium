#' Import Cell Ranger count matrix from HDF5 file
#' 
#' @note Updated 2019-07-31.
#' @export
#' 
#' @inheritParams basejump::params
#' 
#' @seealso `cellrangerRkit::get_matrix_from_h5()`
#' 
#' @return `sparseMatrix`.
#'   Cell barcodes in the columns, features (i.e. genes) in the rows.
#'   
#' @examples
#' ## > x <- importCountsHDF5(file = "filtered_feature_bc_matrix.h5")
#' ## > dim(x)
importCountsHDF5 <-  # nolint
    function(file) {
        assert(
            isAFile(file),
            grepl(pattern = "\\.H5", x = file, ignore.case = TRUE)
        )
        
        name <- names(h5dump(file, load = FALSE))
        assert(isString(name))
        
        ## Import HDF5 data.
        h5 <- h5read(file = file, name = name)
        
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
