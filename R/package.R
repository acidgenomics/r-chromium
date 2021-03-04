#' Chromium
#'
#' Toolkit for 10X Genomics Chromium single cell data.
#'
#' @aliases NULL
#' @keywords internal
#' 
#' @importFrom BiocParallel SerialParam bplapply bpmapply bpparam
#' @importFrom basejump DataFrame SimpleList SingleCellExperiment assay assays
#'   calculateMetrics camelCase detectLanes droplevels emptyRanges import
#'   importSampleData leftJoin makeDimnames makeGRangesFromEnsembl
#'   makeGRangesFromGFF makeDimnames makeNames makeSingleCellExperiment
#'   mapCellsToSamples mcols metadata metadata<- metricsCols minimalSampleData
#'   packageName packageVersion printString realpath sparseMatrix
#'   standardizeCall
#' @importFrom cli cat_line cli_alert cli_alert_info cli_alert_success
#'   cli_alert_warning cli_h1 cli_h2 cli_text cli_ul
#' @importFrom goalie allAreDirectories allAreFiles allAreMatchingRegex
#'   areDisjointSets areSetEqual assert hasLength hasNames hasRownames
#'   hasValidDimnames hasValidNames isADirectory isAFile isAny isCharacter
#'   isFlag isInt isScalar isString isSubset validNames validate validateClasses
#' @importFrom methods as is new setClass setMethod setValidity validObject
#' @importFrom rhdf5 h5dump h5read
#' @importFrom stringr str_match str_split str_trunc
"_PACKAGE"
