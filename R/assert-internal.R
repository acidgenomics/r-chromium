#' Does a given directory path contain a subdirectory name?
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



#' Does the dataset contain a single sample?
#'
#' This is a utility function intended to make loading of example datasets
#' from the 10X Genomics website easier.
#'
#' @param dir Cell Ranger output directory.
#'
#' @note Updated 2019-08-01.
#' @noRd
.isSingleSample <- function(dir) {
    assert(isADirectory(dir))

    ## Check for matrix in the top level of the directory.
    files <- list.files(path = dir, recursive = FALSE)
    ok <- any(grepl(pattern = "matrix\\.(h5|mtx)(\\.gz)?$", x = files))
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
