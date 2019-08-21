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
        colData <- colData(object)
        metadata <- metadata(object)
        sampleData <- sampleData(object)

        ## Assays --------------------------------------------------------------
        ok <- validate(isSubset("counts", names(assays(object))))
        if (!isTRUE(ok)) return(ok)

        ## Row data ------------------------------------------------------------
        ok <- validate(
            is(rowRanges(object), "GRanges"),
            is(rowData(object), "DataFrame")
        )
        if (!isTRUE(ok)) return(ok)

        ## Column data ---------------------------------------------------------
        ok <- validate(
            ## Require that metrics columns are defined.
            isSubset(metricsCols, colnames(colData)),
            ## Ensure that `interestingGroups` isn't slotted in colData.
            areDisjointSets("interestingGroups", colnames(colData))
        )
        if (!isTRUE(ok)) return(ok)

        ## Metadata ------------------------------------------------------------
        ok <- validateClasses(
            object = metadata,
            expected = list(
                allSamples = "logical",
                call = "call",
                dir = "character",
                date = "Date",
                ensemblRelease = "integer",
                genomeBuild = "character",
                gffFile = "character",
                interestingGroups = "character",
                lanes = "integer",
                level = "character",
                matrixFiles = "character",
                organism = "character",
                pipeline = "character",
                refJSON = "list",
                refdataDir = "character",
                sampleDirs = "character",
                sampleMetadataFile = "character",
                sessionInfo = "session_info",
                umiType = "character",
                version = "package_version",
                wd = "character"
            ),
            subset = TRUE
        )
        if (!isTRUE(ok)) return(ok)

        ok <- validate(
            !isSubset("sampleName", names(metadata)),
            isSubset(metadata[["level"]], c("genes", "transcripts"))
        )
        if (!isTRUE(ok)) return(ok)

        TRUE
    }
)



## nolint start
##
## #' Cell Ranger ATAC-seq data set
## #'
## #' Contains single-cell ATAC-seq data.
## #'
## #' @author Michael Steinbaugh
## #' @export
## #'
## #' @examples
## #' showClass("CellRangerATAC")
## setClass(
##     Class = "CellRangerATAC",
##     contains = "SingleCellExperiment"
## )
## setValidity(
##     Class = "CellRangerATAC",
##     method = function(object) {
##         TRUE
##     }
## )
##
## nolint end
