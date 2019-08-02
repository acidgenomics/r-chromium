#' 10X Genomics Chromium
#' 
#' Extends `SingleCellExperiment`.
#' 
#' @author Michael Steinbaugh
#' @export
#' 
#' @examples
#' showClass("Chromium")
setClass(
    Class = "Chromium",
    contains = "SingleCellExperiment"
)
setValidity(
    Class = "Chromium",
    method = function(object) {
        ## FIXME
        TRUE
    }
)



#' Cell Ranger RNA-seq data set
#' 
#' Contains UMI droplet-based single-cell RNA-seq data.
#'
#' @author Michael Steinbaugh
#' @export
#' 
#' @examples
#' showClass("CellRangerRNASeq")
setClass(
    Class = "CellRangerRNASeq",
    contains = "Chromium"
)
setValidity(
    Class = "CellRangerRNASeq",
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
#' showClass("CellRangerATACSeq")
setClass(
    Class = "CellRangerATACSeq",
    contains = "Chromium",
    validity = function(object) {
        ## FIXME
        TRUE
    }
)
