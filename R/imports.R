#' @importFrom Matrix readMM sparseMatrix
#' @importFrom SingleCellExperiment isSpike
#' @importFrom SummarizedExperiment assays
#' @importFrom basejump camel detectLanes emptyRanges import makeDimnames
#'   makeGRangesFromEnsembl makeGRangesFromGFF makeNames
#'   makeSingleCellExperiment mapCellsToSamples minimalSampleData realpath
#' @importFrom bcbioBase readSampleData
#' @importFrom bcbioSingleCell calculateMetrics
#' @importFrom dplyr mutate_if pull
#' @importFrom goalie assert isADirectory
#' @importFrom magrittr %>%
#' @importFrom methods as new validObject
#' @importFrom readr read_lines read_tsv
#' @importFrom rhdf5 h5dump h5read
#' @importFrom stringr str_match str_split
#' @importFrom utils globalVariables packageVersion
NULL
