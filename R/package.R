#' Chromium
#' 
#' Toolkit for 10X Genomics Chromium single-cell RNA-seq data.
#' 
#' @aliases NULL
#' @keywords internal
#' 
#' @importFrom Matrix readMM sparseMatrix
#' @importFrom SingleCellExperiment isSpike
#' @importFrom SummarizedExperiment assays
#' @importFrom assertive.base assert_are_identical
#' @importFrom assertive.files assert_all_are_dirs assert_all_are_existing_files
#' @importFrom assertive.properties assert_has_names assert_is_non_empty
#' @importFrom assertive.sets assert_is_subset
#' @importFrom assertive.types assert_is_a_bool assert_is_a_string
#'   assert_is_all_of assert_is_any_of assert_is_character is_a_string
#' @importFrom basejump camel detectLanes emptyRanges import makeDimnames
#'   makeGRangesFromEnsembl makeGRangesFromGFF makeNames
#'   makeSingleCellExperiment mapCellsToSamples minimalSampleData realpath
#' @importFrom dplyr mutate_if pull
#' @importFrom goalie assertIsAnImplicitIntegerOrNULL assertIsStringOrNULL
#' @importFrom magrittr %>%
#' @importFrom methods as new validObject
#' @importFrom readr read_lines read_tsv
#' @importFrom rhdf5 h5dump h5read
#' @importFrom rlang has_length
#' @importFrom stringr str_match str_split
#' @importFrom utils globalVariables packageVersion
"_PACKAGE"

globalVariables(".")

# FIXME Need to fix these functions shared with bcbioSingleCell:
# - .updateMetadata
# - calculateMetrics
# - readSampleData
# Look for `requireNamespace()` calls and improve...
