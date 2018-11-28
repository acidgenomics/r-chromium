#' 10X Genomics Cell Ranger Data Set
#'
#' Extends `SingleCellExperiment`, with additional validity checks on the
#' `metadata()` slot.
#'
#' @family S4 classes
#' @author Michael Steinbaugh
#' @export
#'
#' @seealso `CellRanger()`.
setClass(
    Class = "CellRanger",
    contains = "SingleCellExperiment"
)
setValidity(
    Class = "CellRanger",
    method = function(object) {
        TRUE
    }
)
