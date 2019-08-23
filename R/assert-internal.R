#' Does the directory contain Cell Ranger aggr output?
#' @note Updated 2019-08-22.
#' @noRd
.isAggregate <- function(dir) {
    isADirectory(file.path(dir, "SC_RNA_AGGREGATOR_CS"))
}



#' Does a given directory path contain a subdirectory name?
#'
#' @note Updated 2019-08-22.
#'
#' @param paths `character`.
#'   Directory paths. Parameterized.
#' @param name `character(1)`.
#'   Subdirectory name.
#'
#' @examples
#' .hasSubdir(path = "sample_0001", name = "SC_RNA_COUNTER_CS")
#'
#' @noRd
.hasSubdir <- function(paths, name) {
    assert(
        allAreDirectories(paths),
        isString(name)
    )
    out <- vapply(
        X = paths,
        FUN = function(path, name) {
            dir.exists(file.path(path, name))
        },
        FUN.VALUE = logical(1L),
        USE.NAMES = FALSE,
        name = name
    )
    names(out) <- basename(paths)
    out
}



#' Does the dataset contain a minimal, single directory structure?
#' @note Updated 2019-08-22.
#' @noRd
.isMinimalSample <- function(dir) {
    files <- list.files(
        path = dir,
        pattern = "\\.(h5|mtx)(\\.gz)?",
        recursive = FALSE
    )
    isString(files)
}



#' Does the dataset contain a single sample?
#'
#' This is a utility function intended to make loading of example datasets
#' from the 10X Genomics website easier.
#'
#' @note Updated 2019-08-22.
#'
#' @param dir Cell Ranger output directory.
#'
#' @note Updated 2019-08-01.
#' @noRd
.isSingleSample <- function(dir) {
    assert(isADirectory(dir))
    ## Check for minimal, single sample directory structure.
    ok <- .isMinimalSample(dir)
    if (isTRUE(ok)) return(ok)
    ## Check matrix files in `outs/` subdirectory.
    outsDir <- file.path(dir, "outs")
    if (dir.exists(outsDir)) {
        files <- list.files(path = outsDir, recursive = FALSE)
        ## Note that Cell Ranger v2 outputs "matrices" instead of "matrix".
        ok <- any(grepl("^(filtered|raw)_", files))
        if (isTRUE(ok)) return(ok)
    }
    FALSE
}
