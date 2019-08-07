#' Cell Ranger RNA-seq data set
#' 
#' Contains UMI droplet-based single-cell RNA-seq data.
#'
#' @author Michael Steinbaugh
#' @export
#' 
#' @examples
#' showClass("CellRanger")
setClass(
    Class = "CellRanger",
    contains = "SingleCellExperiment"
)
setValidity(
    Class = "CellRanger",
    method = function(object) {
        ## FIXME
        TRUE
    }
)



#' Cell Ranger ATAC-seq data set
#' 
#' Contains single-cell ATAC-seq data.
#'
#' @author Michael Steinbaugh
#' @export
#' 
#' @examples
#' showClass("CellRangerATAC")
setClass(
    Class = "CellRangerATAC",
    contains = "SingleCellExperiment",
    validity = function(object) {
        ## FIXME
        TRUE
    }
)
