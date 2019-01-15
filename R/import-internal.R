.import <-  # nolint
    function(sampleFiles) {
        assert(
            allAreFiles(sampleFiles),
            hasNames(sampleFiles)
        )

        message("Importing counts.")

        if (all(grepl("\\.mtx$", sampleFiles))) {
            fun <- .import.mtx
        } else if (all(grepl("\\.h5$", sampleFiles))) {
            fun <- .import.h5
        }

        list <- mapply(
            sampleID = names(sampleFiles),
            file = sampleFiles,
            FUN = function(sampleID, file) {
                counts <- fun(file)

                # Strip index when all barcodes end with "-1".
                if (all(grepl("-1$", colnames(counts)))) {
                    colnames(counts) <- sub("-1", "", colnames(counts))
                }

                # Now move the multiplexed index name/number to the beginning,
                # for more logical sorting and consistency with bcbio approach.
                colnames(counts) <- sub(
                    pattern = "^([ACGT]+)-(.+)$",
                    replacement = "\\2-\\1",
                    x = colnames(counts)
                )

                # Prefix cell barcodes with sample identifier when we're loading
                # counts from multiple samples.
                if (
                    length(sampleFiles) > 1L ||
                    grepl("^([[:digit:]]+)-([ACGT]+)$", colnames(counts))
                ) {
                    colnames(counts) <-
                        paste(sampleID, colnames(counts), sep = "_")
                }
                # Ensure names are valid.
                counts <- makeDimnames(counts)
                counts
            },
            SIMPLIFY = FALSE,
            USE.NAMES = TRUE
        )

        # Bind the matrices.
        do.call(cbind, list)
    }



# Import Cell Ranger HDF5 Counts
# @seealso `cellrangerRkit::get_matrix_from_h5`
.import.h5 <-  # nolint
    function(file) {
        assert(isAFile(file))

        # Get the genome build, which we need to pass as "name".
        genomeBuild <- names(h5dump(file, load = FALSE))
        assert(isString(genomeBuild))

        # Use genome build name (e.g. "/GRCh38/data").
        h5 <- h5read(file = file, name = genomeBuild)

        # Want `Csparse` not `Tsparse` matrix.
        counts <- sparseMatrix(
            i = h5[["indices"]] + 1L,
            p = h5[["indptr"]],
            x = as.numeric(h5[["data"]]),
            dims = h5[["shape"]],
            giveCsparse = TRUE
        )

        rownames <- h5[["genes"]]
        colnames <- h5[["barcodes"]]

        assert(
            identical(length(rownames), nrow(counts)),
            identical(length(colnames), ncol(counts))
        )

        rownames(counts) <- rownames
        colnames(counts) <- colnames

        counts
    }



# Import Cell Ranger Sparse Counts
# Matrix Market Exchange (MEX/MTX) format.
.import.mtx <-  # nolint
    function(file) {
        assert(isAFile(file))

        # Locate required sidecar files.
        barcodesFile <- file.path(dirname(file), "barcodes.tsv")
        genesFile <- file.path(dirname(file), "genes.tsv")
        assert(allAreFiles(c(barcodesFile, genesFile)))

        # `genes.tsv` is tab delimited.
        rownames <- read_tsv(
            file = genesFile,
            col_names = c("geneID", "geneName"),
            col_types = "cc"
        ) %>%
            pull("geneID")

        # `barcodes.tsv` is not tab delimited.
        colnames <- read_lines(barcodesFile)

        counts <- readMM(file)

        assert(
            identical(length(rownames), nrow(counts)),
            identical(length(colnames), ncol(counts))
        )

        rownames(counts) <- rownames
        colnames(counts) <- colnames

        counts
    }
