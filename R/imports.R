#' @importFrom BiocParallel SerialParam bplapply bpmapply bpparam
#' @importFrom Matrix sparseMatrix
#' @importFrom SingleCellExperiment SingleCellExperiment isSpike
#' @importFrom SummarizedExperiment assays
#' @importFrom S4Vectors SimpleList mcols
#' @importFrom basejump calculateMetrics camel detectLanes droplevels
#'   emptyRanges import makeDimnames makeGRangesFromEnsembl makeGRangesFromGFF
#'   makeDimnames makeNames makeSingleCellExperiment mapCellsToSamples
#'   minimalSampleData readSampleData realpath
#' @importFrom goalie allAreDirectories allAreFiles assert hasLength hasNames
#'   hasValidDimnames hasValidNames isADirectory isAFile isAny isCharacter
#'   isFlag isInt isNonEmpty isScalar isString isSubset validNames
#' @importFrom methods as is new setClass setMethod setValidity validObject
#' @importFrom rhdf5 h5dump h5read
#' @importFrom stringr str_match str_split str_trunc
#' @importFrom utils globalVariables packageVersion
NULL
