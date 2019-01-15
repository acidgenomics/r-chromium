#' Chromium
#' 
#' Toolkit for 10X Genomics single-cell RNA-seq data.
#' 
#' @aliases NULL
#' @keywords internal
"_PACKAGE"

globalVariables(".")

# FIXME Need to fix these functions shared with bcbioSingleCell:
# - .updateMetadata
# - calculateMetrics
# - readSampleData
# Look for `requireNamespace` calls and improve...
