#' Chromium
#'
#' Toolkit for 10X Genomics Chromium single cell data.
#'
#' @aliases NULL
#' @keywords internal
"_PACKAGE"



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
#' strMatch
#' @importFrom AcidCLI abort alert alertInfo alertSuccess alertWarning txt ul
#' @importFrom AcidExperiment detectLanes importSampleData minimalSampleData
#' @importFrom AcidGenomes emptyRanges makeGRangesFromEnsembl makeGRangesFromGff
#' @importFrom AcidSingleCell makeSingleCellExperiment mapCellsToSamples
#' @importFrom Matrix sparseMatrix
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom goalie allAreDirs allAreFiles allAreMatchingRegex
#' areDisjointSets areSetEqual assert hasLength hasNames hasRownames
#' hasValidDimnames hasValidNames isADir isAFile isAny isCharacter isFlag isInt
#' isScalar isString isSubset validNames validate validateClasses
#' @importFrom methods as is new setClass setMethod setValidity validObject
#' @importFrom rhdf5 h5dump h5read
#' @importFrom utils packageName packageVersion
NULL
