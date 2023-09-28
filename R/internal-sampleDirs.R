#' Determine which subdirectories contain sample files.
#'
#' Checks for the presence of nested `SC_RNA_COUNTER_CS` directories.
#'
#' @note `aggr` returns `SC_RNA_AGGREGATOR_CS/` directory.
#' @note Updated 2023-09-28.
#'
#' @param dir Cell Ranger output directory.
#'
#' @return `character`.
#' Directory paths that contain individual scRNA-seq sample files.
#' Note that aggregate directories from `aggr` are excluded here.
#'
#' @noRd
.sampleDirs <-
    function(dir, filtered) {
        assert(
            isADir(dir),
            isFlag(filtered)
        )
        ## Check for single sample mode, used for 10X example datasets.
        if (isAFile(tryCatch(
            expr = .findMatrixFile(dir = dir, filtered = filtered),
            error = function(e) {
                NULL
            }
        ))) {
            dir <- realpath(dir)
            names(dir) <- makeNames(basename(dir))
            return(dir)
        }
        dirs <- sort(list.dirs(
            path = dir,
            full.names = TRUE,
            recursive = FALSE
        ))
        assert(
            hasLength(dirs),
            msg = sprintf("Failed to detect samples in {.dir %s}.", dir)
        )
        ## Must contain `SC_RNA_COUNTER_CS` subdirectory.
        subdir <- "SC_RNA_COUNTER_CS"
        keep <- .hasSubdir(paths = dirs, name = subdir)
        assert(
            any(keep),
            msg = sprintf(
                "No sample subdirectories containing {.var %s}.",
                subdir
            )
        )
        dirs <- dirs[keep]
        ## Must contain `outs` subdirectory.
        subdir <- "outs"
        keep <- .hasSubdir(paths = dirs, name = subdir)
        assert(
            any(keep),
            msg = sprintf(
                "No sample subdirectories containing {.var %s}.",
                subdir
            )
        )
        dirs <- dirs[keep]
        assert(allAreDirs(dirs))
        names(dirs) <- makeNames(basename(dirs))
        alertInfo(sprintf(
            fmt = "%d %s detected:",
            length(dirs),
            ngettext(n = length(dirs), msg1 = "sample", msg2 = "samples")
        ))
        ul(names(dirs))
        dirs
    }
