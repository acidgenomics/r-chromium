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
    
    assert(allAreDirectories(dirs))
    dirs
}
