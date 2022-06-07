## S4 generics and methods =====================================================

#' @importFrom AcidGenerics calculateMetrics camelCase droplevels2 leftJoin
#' makeDimnames makeNames sampleData
#' @importFrom S4Vectors mcols metadata metadata<-
#' @importFrom SummarizedExperiment assay assays colData rowData rowRanges
#' @importFrom pipette import
#'
#' @importMethodsFrom AcidExperiment calculateMetrics droplevels2 sampleData
#' @importMethodsFrom AcidPlyr leftJoin
#' @importMethodsFrom pipette import
#' @importMethodsFrom syntactic camelCase makeDimnames makeNames
NULL



## Standard functions ==========================================================

#' @importFrom AcidBase metricsCols printString realpath standardizeCall
#' @importFrom AcidCLI abort alert alertInfo alertSuccess alertWarning txt ul
#' @importFrom AcidExperiment detectLanes importSampleData minimalSampleData
#' @importFrom AcidGenomes emptyRanges makeGRangesFromEnsembl makeGRangesFromGFF
#' @importFrom AcidSingleCell makeSingleCellExperiment mapCellsToSamples
#' @importFrom Matrix sparseMatrix
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom goalie allAreDirectories allAreFiles allAreMatchingRegex
#' areDisjointSets areSetEqual assert hasLength hasNames hasRownames
#' hasValidDimnames hasValidNames isADirectory isAFile isAny isCharacter
#' isFlag isInt isScalar isString isSubset validNames validate validateClasses
#' @importFrom methods as is new setClass setMethod setValidity validObject
#' @importFrom rhdf5 h5dump h5read
#' @importFrom utils packageName packageVersion
NULL



## FIXME Switch this to stringi package.

#' @importFrom stringr str_match str_split
NULL
