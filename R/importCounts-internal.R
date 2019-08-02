#' Import counts from either HDF5 or MTX file.
#' @note Updated 2019-08-01.
#' @noRd
.importCounts <- function(file, ...) {
    assert(isAFile(file))
    if (grepl("\\.H5", file, ignore.case = TRUE)) {
        fun <- importCountsHDF5
    } else if (grepl("\\.MTX", file, ignore.case = TRUE)) {
        fun <- importCountsMTX
    } else {
        stop(sprintf("Failed to import '%s'.", file))
    }
    fun(file, ...)
}
