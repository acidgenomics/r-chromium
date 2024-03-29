#' Find the count matrix inside a Cell Ranger sample directory.
#'
#' This function will automatically move inside an `outs/` directory inside
#' the path if it is detected.
#'
#' Currently preferring HDF5 over MTX.
#'
#' @section File name history:
#'
#' Cell Ranger v3
#'
#' - H5: `filtered_feature_bc_matrix.h5`.
#' - MTX: `filtered_feature_bc_matrix/matrix.mtx.gz`.
#'
#' Cell Ranger v2
#'
#' - H5: `filtered_gene_bc_matrices_h5.h5`
#' - MTX: `filtered_gene_bc_matrices/<genomeBuild>/matrix.mtx`
#'
#' @param dir Sample directory.
#'
#' @param filtered `logical(1)`.
#' - `TRUE`: Look for `filtered_*` matrix.
#' - `FALSE`: Look for `raw_*` matrix.
#'
#' Doesn't apply if there's only a single matrix file in the directory.
#'
#' @note Updated 2023-09-28.
#' @noRd
.findMatrixFile <-
    function(dir, filtered) {
        assert(
            isADir(dir),
            isFlag(filtered)
        )
        dir <- realpath(dir)
        outsDir <- file.path(dir, "outs")
        ## Simple mode ---------------------------------------------------------
        ## For minimal examples and data downloaded from 10X website.
        if (!isADir(outsDir)) {
            ## First, look to see if the input contains a supported HDF5 or MTX
            ## file directly.
            matFilePattern <- paste0(
                "^.+_",
                ifelse(
                    test = filtered,
                    yes = "filtered",
                    no = "raw"
                ),
                "_.+\\.(h5|mtx)(\\.gz)?$"
            )
            matFile <- list.files(
                path = dir,
                pattern = matFilePattern,
                full.names = TRUE,
                recursive = FALSE,
                ignore.case = TRUE
            )
            if (isAFile(matFile)) {
                return(matFile)
            }
            matDirPattern <- paste0(
                "^",
                ifelse(
                    test = filtered,
                    yes = "filtered",
                    no = "raw"
                ),
                "_(",
                paste(
                    "feature_bc_matrix",
                    "gene_bc_matrices",
                    sep = "|"
                ),
                ")$"
            )
            matDir <- list.files(
                path = dir,
                pattern = matDirPattern,
                full.names = TRUE,
                recursive = FALSE,
                ignore.case = FALSE
            )
            if (isADir(matDir)) {
                matFile <- list.files(
                    path = matDir,
                    pattern = "^matrix\\.mtx(\\.gz)?$",
                    full.names = TRUE,
                    recursive = TRUE,
                    ignore.case = TRUE
                )
                if (isAFile(matFile)) {
                    return(matFile)
                }
            }
            abort(sprintf("Failed to detect matrix in {.dir %s}.", dir))
        }
        ## Standard Cell Ranger output -----------------------------------------
        ## Recurse into `outs/` directory by default.
        dir <- outsDir
        if (isTRUE(filtered)) {
            prefix <- "filtered"
        } else {
            prefix <- "raw"
        }
        files <- list.files(
            path = dir,
            pattern = paste0("^", prefix, "_"),
            recursive = FALSE,
            full.names = FALSE,
            ignore.case = FALSE
        )
        assert(hasLength(files))
        ## Get the Cell Ranger version, based on the file names.
        if (isTRUE(any(grepl(
            pattern = paste0("^", prefix, "_feature_bc_matrix$"),
            x = files
        )))) {
            version <- "3"
            filestem <- paste0(prefix, "_feature_bc_matrix")
        } else if (isTRUE(any(grepl(
            pattern = paste0("^", prefix, "_gene_bc_matrices$"),
            x = files
        )))) {
            version <- "2"
            filestem <- paste0(prefix, "_gene_bc_matrices")
        } else {
            abort("Failed to detect Cell Ranger version based on file names.")
        }
        version <- numeric_version(version)
        ## Currently preferring HDF5 over MTX.
        if (isTRUE(
            file.exists(file.path(dir, paste0(filestem, ".h5")))
        )) {
            file <- file.path(dir, paste0(filestem, ".h5"))
            attr(file, "pipeline") <- "Cell Ranger v3 HDF5"
        } else if (isTRUE(
            file.exists(file.path(dir, paste0(filestem, "_h5.h5")))
        )) {
            file <- file.path(dir, paste0(filestem, "_h5.h5"))
            attr(file, "pipeline") <- "Cell Ranger v2 HDF5"
        } else if (isTRUE(
            file.exists(file.path(dir, filestem, "matrix.mtx.gz"))
        )) {
            file <- file.path(dir, filestem, "matrix.mtx.gz")
            attr(file, "pipeline") <- "Cell Ranger v3 MTX"
        } else if (isADir(file.path(dir, filestem))) {
            ## Get the genome build from the first sample directory.
            genomeBuild <- list.dirs(
                path = file.path(dir, filestem),
                full.names = FALSE,
                recursive = FALSE
            )
            assert(isString(genomeBuild))
            file <- file.path(dir, filestem, genomeBuild, "matrix.mtx")
            attr(file, "pipeline") <- "Cell Ranger v2 MTX"
        } else {
            abort("Failed to locate count matrix.")
        }
        assert(
            isAFile(file),
            isString(attr(file, "pipeline"))
        )
        file
    }



#' Find all matrix files for a data set
#'
#' @note Updated 2022-06-07.
#' @noRd
.matrixFiles <-
    function(sampleDirs,
             filtered) {
        list <- lapply(
            X = sampleDirs,
            FUN = .findMatrixFile,
            filtered = filtered
        )
        ## Note that `unlist()` call drops attributes on the file paths.
        pipeline <- vapply(
            X = list,
            FUN = function(x) {
                attr <- attr(x, "pipeline")
                ## Handle samples loaded from simple mode.
                if (!hasLength(attr)) {
                    attr <- NA_character_
                }
                attr
            },
            FUN.VALUE = character(1L),
            USE.NAMES = TRUE
        )
        files <- unlist(list, use.names = TRUE)
        assert(
            allAreFiles(files),
            hasLength(unique(basename(files)), n = 1L),
            hasLength(unique(pipeline), n = 1L),
            hasValidNames(files)
        )
        attr(files, "pipeline") <- unique(pipeline)
        files
    }
