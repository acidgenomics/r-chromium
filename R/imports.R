#' @importFrom Matrix readMM sparseMatrix
#' @importFrom SingleCellExperiment SingleCellExperiment isSpike
#' @importFrom SummarizedExperiment assays
#' @importFrom S4Vectors SimpleList mcols
#' @importFrom basejump camel detectLanes emptyRanges import makeDimnames
#'   makeGRangesFromEnsembl makeGRangesFromGFF makeNames
#'   makeSingleCellExperiment mapCellsToSamples minimalSampleData realpath
#'   relevelColData relevelRowRanges
#' @importFrom bcbioBase readSampleData
#' @importFrom bcbioSingleCell calculateMetrics
#' @importFrom dplyr mutate_if pull
#' @importFrom goalie allAreDirectories allAreFiles assert hasLength hasNames
#'   hasValidDimnames isADirectory isAFile isCharacter isFlag isInt isNonEmpty
#'   isScalar isString isSubset
#' @importFrom magrittr %>%
#' @importFrom methods as is new setClass setValidity validObject
#' @importFrom readr read_lines read_tsv
#' @importFrom rhdf5 h5dump h5read
#' @importFrom stringr str_match str_split
#' @importFrom utils globalVariables packageVersion
NULL

globalVariables(".")
