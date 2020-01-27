#' Determine which subdirectories contain sample files.
#'
#' Checks for the presence of nested `SC_RNA_COUNTER_CS` directories.
#'
#' @note `aggr` returns `SC_RNA_AGGREGATOR_CS/` directory.
#' @note Updated 2020-01-26.
#'
#' @param dir Cell Ranger output directory.
#'
#' @return `character`.
#'   Directory paths that contain individual scRNA-seq sample files.
#'   Note that aggregate directories from `aggr` are excluded here.
#'
#' @noRd
.sampleDirs <- function(dir) {
    ## Check for single sample mode, used for 10X example datasets.
    if (isAFile(tryCatch(
        expr = .findMatrixFile(dir),
        error = function(e) NULL
    ))) {
        dir <- realpath(dir)
        names(dir) <- makeNames(basename(dir))
        return(dir)
    }
    dirs <- sort(list.dirs(path = dir, full.names = TRUE, recursive = FALSE))
    assert(hasLength(dirs))
    ## Must contain `SC_RNA_COUNTER_CS` subdirectory.
    subdir <- "SC_RNA_COUNTER_CS"
    keep <- .hasSubdir(paths = dirs, name = subdir)
    if (!any(keep)) {
        stop(sprintf(
            fmt = "No sample subdirectories containing '%s'.",
            subdir
        ))
    }
    dirs <- dirs[keep]
    ## Must contain `outs` subdirectory.
    subdir <- "outs"
    keep <- .hasSubdir(paths = dirs, name = subdir)
    if (!any(keep)) {
        stop(sprintf(
            fmt = "No sample subdirectories containing '%s'.",
            subdir
        ))
    }
    dirs <- dirs[keep]
    assert(allAreDirectories(dirs))
    names(dirs) <- makeNames(basename(dirs))
    cli_text(sprintf(
        fmt = "%d %s detected:",
        length(dirs),
        ngettext(n = length(dirs), msg1 = "sample", msg2 = "samples")
    ))
    cli_ul(names(dirs))
    dirs
}
