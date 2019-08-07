#' Determine which subdirectories contain sample files.
#' 
#' Checks for the presence of nested `SC_RNA_COUNTER_CS` directories.
#' 
#' @note `aggr` returns `SC_RNA_AGGREGATOR_CS/` directory.
#' @note Updated 2019-08-07.
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
        expr = .findCountMatrix(dir),
        error = function(e) NULL
    ))) {
        dir <- realpath(dir)
        names(dir) <- makeNames(basename(dir))
        return(dir)
    }
    
    dirs <- sort(list.dirs(path = dir, full.names = TRUE, recursive = FALSE))
    assert(hasLength(dirs))
    
    ## Must contain `SC_RNA_COUNTER_CS` subdirectory.
    keep <- .hasSubdir(paths = dirs, name = "SC_RNA_COUNTER_CS")
    assert(any(keep))
    dirs <- dirs[keep]
    
    ## Must contain `outs` subdirectory.
    keep <- .hasSubdir(paths = dirs, name = "outs")
    assert(any(keep))
    dirs <- dirs[keep]
    names(dirs) <- makeNames(basename(dirs))
    
    assert(allAreDirectories(dirs))
    message(sprintf(
        fmt = "%d %s detected:\n%s",
        length(dirs),
        ngettext(n = length(dirs), msg1 = "sample", msg2 = "samples"),
        printString(names(dirs))
    ))
    dirs
}
