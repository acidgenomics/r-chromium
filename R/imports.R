#' @importFrom BiocParallel SerialParam bplapply bpmapply bpparam
#' @importFrom Matrix sparseMatrix
#' @importFrom SingleCellExperiment SingleCellExperiment isSpike
#' @importFrom SummarizedExperiment assay assays
#' @importFrom S4Vectors DataFrame SimpleList mcols metadata metadata<-
#' @importFrom basejump calculateMetrics camelCase detectLanes droplevels
#'   emptyRanges import importSampleData leftJoin makeDimnames
#'   makeGRangesFromEnsembl makeGRangesFromGFF makeDimnames makeNames
#'   makeSingleCellExperiment mapCellsToSamples metricsCols minimalSampleData
#'   printString realpath standardizeCall
#' @importFrom goalie allAreDirectories allAreFiles allAreMatchingRegex
#'   areDisjointSets areSetEqual assert hasLength hasNames hasRownames
#'   hasValidDimnames hasValidNames isADirectory isAFile isAny isCharacter
#'   isFlag isInt isNonEmpty isScalar isString isSubset validNames validate
#'   validateClasses
#' @importFrom methods as is new setClass setMethod setValidity validObject
#' @importFrom rhdf5 h5dump h5read
#' @importFrom stringr str_match str_split str_trunc
#' @importFrom utils globalVariables packageVersion
NULL
