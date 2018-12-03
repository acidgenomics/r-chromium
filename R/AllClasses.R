#' 10X Genomics Chromium Single-Cell RNA-Seq Data Set
#' 
#' Contains single-cell RNA-seq counts processed by the Cell Ranger pipeline.
#'
#' Extends `SingleCellExperiment`, with additional validity checks on the
#' `metadata()` slot.
#'
#' @family S4 classes
#' @author Michael Steinbaugh
#' @export
#'
#' @seealso `Chromium()`.
setClass(
    Class = "Chromium",
    contains = "SingleCellExperiment"
)
# FIXME Work on setting actual validity checks here.
setValidity(
    Class = "Chromium",
    method = function(object) {
        TRUE
    }
)
