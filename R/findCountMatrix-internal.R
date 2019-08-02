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
#' @param filtered `logical(1)`.
#'   - `TRUE`: Look for `filtered_*` matrix.
#'   - `FALSE`: Look for `raw_*` matrix.
#'   
#'   Doesn't apply if there's only a single matrix file in the directory.
#' 
#' @note Updated 2019-08-01.
#' @noRd
.findCountMatrix <- function(dir, filtered = TRUE) {
    assert(
        isADirectory(dir),
        isFlag(filtered)
    )
    path <- realpath(dir)
    
    ## Simple mode -------------------------------------------------------------
    ## For minimal examples and data downloaded from 10X website.
    if (!dir.exists(file.path(path, "outs"))) {
        file <- list.files(
            path = path,
            pattern = "matrix\\.(h5|mtx)(\\.gz)?",
            full.names = TRUE
        )
        assert(isAFile(file))
        return(file)
    }
    
    ## Standard Cell Ranger output ---------------------------------------------
    ## Recurse into `outs/` directory by default.
    path <- file.path(path, "outs")
    assert(isADirectory(path))
    
    if (isTRUE(filtered)) {
        prefix <- "filtered"
    } else {
        prefix <- "raw"
    }
    
    files <- list.files(
        path = path,
        pattern = paste0("^", prefix, "_"),
        recursive = FALSE,
        full.names = FALSE
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
        stop("Failed to detect Cell Ranger version based on file names.")
    }
    version <- numeric_version(version)
    
    ## Currently preferring HDF5 over MTX.
    if (isTRUE(
        file.exists(file.path(path, paste0(filestem, ".h5")))
    )) {
        ## v3 HDF5
        file <- file.path(path, paste0(filestem, ".h5"))
    } else if (isTRUE(
        file.exists(file.path(path, paste0(filestem, "_h5.h5")))
    )) {
        ## v2 HDF5
        file <- file.path(path, paste0(filestem, "_h5.h5"))
    } else if (isTRUE(
        file.exists(file.path(path, filestem, "matrix.mtx.gz"))
    )) {
        ## v3 MTX
        file <- file.path(path, filestem, "matrix.mtx.gz")
    } else if (isTRUE(
        dir.exists(file.path(path, filestem))
    )) {
        ## v2 MTX
        ## Get the genome build from the first sample directory.
        genomeBuild <- list.dirs(
            path = path[[1L]],
            full.names = FALSE,
            recursive = FALSE
        )
        assert(isString(genomeBuild))
        file <- file.path(path, genomeBuild, "matrix.mtx")
    }
    
    assert(isAFile(file))
    file
}
