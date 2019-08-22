#' Import Cell Ranger sample files.
#' @note Updated 2019-08-22.
#' @param sampleFiles Count matrix files.
#' @noRd
.importSamples <-  # nolint
    function(sampleFiles) {
        assert(
            allAreFiles(sampleFiles),
            hasNames(sampleFiles)
        )
        message("Importing counts.")
        if (all(grepl("\\.h5$", sampleFiles))) {
            fun <- .importCountsFromHDF5
        } else if (all(grepl("\\.mtx$", sampleFiles))) {
            fun <- .importCountsFromMTX
        } else {
            stop(
                "Failed to determine which file extension (e.g. H5, MTX) ",
                "to use for count matrix import."
            )
        }
        list <- mapply(
            sampleID = names(sampleFiles),
            file = sampleFiles,
            FUN = function(sampleID, file) {
                counts <- fun(file)
                ## Strip index when all barcodes end with "-1".
                if (all(grepl("-1$", colnames(counts)))) {
                    colnames(counts) <- sub("-1", "", colnames(counts))
                }
                ## Now move the multiplexed index name/number to the beginning,
                ## for more logical sorting and consistency with bcbio approach.
                colnames(counts) <- sub(
                    pattern = "^([ACGT]+)-(.+)$",
                    replacement = "\\2-\\1",
                    x = colnames(counts)
                )
                ## Prefix cell barcodes with sample identifier when we're
                ## loading counts from multiple samples.
                if (
                    length(sampleFiles) > 1L ||
                    grepl("^([[:digit:]]+)-([ACGT]+)$", colnames(counts))
                ) {
                    colnames(counts) <-
                        paste(sampleID, colnames(counts), sep = "_")
                }
                ## Ensure names are valid.
                counts <- makeDimnames(counts)
                counts
            },
            SIMPLIFY = FALSE,
            USE.NAMES = TRUE
        )
        ## Bind the matrices.
        do.call(cbind, list)
    }
