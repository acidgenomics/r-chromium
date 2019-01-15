#' 10X Genomics Chromium data set
#'
#' Extends `SingleCellExperiment`, with additional validity checks on the
#' [`metadata()`][S4Vectors::metadata] slot.
#'
#' @author Michael Steinbaugh
#' @export
#'
#' @seealso `Chromium`.
setClass(
    Class = "Chromium",
    contains = "SingleCellExperiment"
)
setValidity(
    Class = "Chromium",
    method = function(object) {
        TRUE
    }
)
