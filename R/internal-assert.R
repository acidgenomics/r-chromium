#' Does the directory contain Cell Ranger aggr output?
#'
#' @note Updated 2023-09-28.
#' @noRd
.isAggregate <- function(dir) {
    isADir(file.path(dir, "SC_RNA_AGGREGATOR_CS"))
}



#' Does a given directory path contain a subdirectory name?
#'
#' @note Updated 2023-09-28.
#'
#' @param paths `character`.
#' Directory paths. Parameterized.
#'
#' @param name `character(1)`.
#' Subdirectory name.
#'
#' @examples
#' .hasSubdir(path = "sample_0001", name = "SC_RNA_COUNTER_CS")
#'
#' @noRd
.hasSubdir <- function(paths, name) {
    assert(
        allAreDirs(paths),
        isString(name)
    )
    out <- vapply(
        X = paths,
        FUN = function(path, name) {
            isADir(file.path(path, name))
        },
        FUN.VALUE = logical(1L),
        USE.NAMES = FALSE,
        name = name
    )
    names(out) <- basename(paths)
    out
}



#' Does the dataset contain a minimal, single directory structure?
#'
#' @note Updated 2023-09-28.
#' @noRd
.isMinimalSample <- function(dir) {
    ok <- isAFile(file.path(dir, "outs", "metrics_summary.csv"))
    unname(!ok)
}
