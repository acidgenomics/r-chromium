## FIXME Rethink this approach.



# Error in .validate_names(rownames, ans_rownames, "assay rownames()",
# "rowData rownames() / rowRanges names()") : assay rownames() must be NULL or
# a character vector with no attributes



#' Import Cell Ranger simple mode
#' 
#' @param matrixFile `character(1)`.
#'   File path to:
#'   - HDF5 matrix file (e.g. `filtered_feature_bc_matrix.h5`)
#'   - MTX sparse matrix file (e.g. `matrix.mtx.gz`).
#'   
#'   Filtered (**recommended**) or unfiltered matrix is supported.
#'   
#' @return `SingleCellExperiment`.
CellRangerRNASeqSimple <- function(matrixFile) {
    counts <- .importCounts(matrixFile)
    counts <- makeDimnames(counts)
    assert(hasValidDimnames(counts))
    assays <- list(counts = counts)
    `new,CellRangerRNASeq`(assays = assays)
}
