#' Chromium
#'
#' Toolkit for 10X Genomics Chromium single cell data.
#'
#' @aliases NULL
#' @keywords internal
#'
#' @importFrom BiocParallel SerialParam bplapply bpmapply bpparam
#' @importFrom basejump DataFrame SimpleList SingleCellExperiment alert
#'   alertInfo alertSuccess alertWarning assay assays calculateMetrics camelCase
#'   detectLanes droplevels emptyRanges h1 h2 import importSampleData leftJoin
#'   makeDimnames makeGRangesFromEnsembl makeGRangesFromGFF makeDimnames
#'   makeNames makeSingleCellExperiment mapCellsToSamples mcols metadata
#'   metadata<- metricsCols minimalSampleData packageName packageVersion
#'   printString realpath sparseMatrix standardizeCall txt ul
#' @importFrom goalie allAreDirectories allAreFiles allAreMatchingRegex
#'   areDisjointSets areSetEqual assert hasLength hasNames hasRownames
#'   hasValidDimnames hasValidNames isADirectory isAFile isAny isCharacter
#'   isFlag isInt isScalar isString isSubset validNames validate validateClasses
#' @importFrom methods as is new setClass setMethod setValidity validObject
#' @importFrom rhdf5 h5dump h5read
#' @importFrom stringr str_match str_split str_trunc
"_PACKAGE"
