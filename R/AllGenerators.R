#' @export
#' @inherit CellRanger-class title description
#' @note Updated 2023-09-28.
#'
#' @details
#' Read [10x Genomics Cell Ranger](https://www.10xgenomics.com/software/) output
#' for a Chromium data set into a `SingleCellExperiment` object.
#'
#' Currently supports loading of a single genome.
#'
#' @section Directory structure for multiple samples:
#'
#' Cell Ranger can vary in its output directory structure, but we're requiring a
#' single, consistent directory structure for datasets containing multiple
#' samples that have not been aggregated into a single matrix with `aggr`.
#'
#' Cell Ranger v3 output:
#'
#' ```
#' | <dir>/
#' |-- <sampleName>/
#' |---- SC_RNA_COUNTER_CS/
#' |---- outs/
#' |------ filtered_feature_bc_matrix/
#' |-------- barcodes.tsv.gz
#' |-------- features.tsv.gz
#' |-------- matrix.mtx.gz
#' |------ filtered_feature_bc_matrix.h5
#' |------ metrics_summary.csv
#' |------ molecule_info.h5
#' |------ possorted_genome_bam.bam
#' |------ possorted_genome_bam.bam.bai
#' |------ raw_feature_bc_matrix/
#' |-------- barcodes.tsv.gz
#' |-------- features.tsv.gz
#' |-------- matrix.mtx.gz
#' |------ raw_feature_bc_matrix.h5
#' |------ web_summary.html
#' ```
#'
#' Cell Ranger v2 output:
#'
#' ```
#' | <dir>/
#' |-- <sampleName>/
#' |---- SC_RNA_COUNTER_CS/
#' |---- outs/
#' |------ filtered_gene_bc_matrices/
#' |-------- <genomeBuild>/
#' |---------- barcodes.tsv
#' |---------- genes.tsv
#' |---------- matrix.mtx
#' |------ filtered_gene_bc_matrices_h5.h5
#' |------ metrics_summary.csv
#' |------ molecule_info.h5
#' |------ possorted_genome_bam.bam
#' |------ possorted_genome_bam.bam.bai
#' |------ raw_gene_bc_matrices/
#' |-------- <genomeBuild>/
#' |---------- barcodes.tsv
#' |---------- genes.tsv
#' |---------- matrix.mtx
#' |------ raw_gene_bc_matrices_h5.h5
#' ```
#'
#' @section Sample metadata:
#'
#' A user-supplied sample metadata file defined by `sampleMetadataFile` is
#' required for multiplexed datasets. Otherwise this can be left `NULL`, and
#' minimal sample data will be used, based on the directory names.
#'
#' @section Reference data:
#'
#' We strongly recommend supplying the corresponding reference data required for
#' Cell Ranger with the `refdataDir` argument. It will convert the gene
#' annotations defined in the GTF file into a `GRanges` object, which get
#' slotted in [`rowRanges()`][SummarizedExperiment::rowRanges]. Otherwise, the
#' function will attempt to use the most current annotations available from
#' Ensembl, and some gene IDs may not match, due to deprecation in the current
#' Ensembl release.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @param dir `character(1)`.
#' Directory path to Cell Ranger output.
#'
#' @param filtered `logical(1)`.
#' Use filtered (recommended) or raw counts. Note that raw counts still
#' contain only whitelisted cellular barcodes.
#'
#' @param refdataDir `character(1)` or `NULL`.
#' Directory path to Cell Ranger reference annotation data.
#'
#' @return `CellRanger`.
#'
#' @seealso
#' - https://support.10xgenomics.com/single-cell-gene-expression/
#'
#' @examples
#' dir <- system.file("extdata", "cellranger_v3", package = "Chromium")
#' x <- CellRanger(dir)
#' print(x)
CellRanger <- # nolint
    function(dir,
             filtered = TRUE,
             organism = NULL,
             ensemblRelease = NULL,
             genomeBuild = NULL,
             gffFile = NULL,
             refdataDir = NULL,
             samples = NULL,
             censorSamples = NULL,
             sampleMetadataFile = NULL,
             transgeneNames = NULL,
             interestingGroups = "sampleName") {
        assert(
            isADir(dir),
            isFlag(filtered),
            isString(organism, nullOk = TRUE),
            isInt(ensemblRelease, nullOk = TRUE),
            isString(genomeBuild, nullOk = TRUE),
            isString(gffFile, nullOk = TRUE),
            isADir(refdataDir, nullOk = TRUE),
            isAny(samples, classes = c("character", "NULL")),
            isAny(censorSamples, classes = c("character", "NULL")),
            isAFile(sampleMetadataFile, nullOk = TRUE),
            isCharacter(transgeneNames, nullOk = TRUE),
            isCharacter(interestingGroups)
        )
        alert("Importing Chromium single-cell RNA-seq run.")
        ## Run info ------------------------------------------------------------
        level <- "genes"
        dir <- realpath(dir)
        if (isADir(refdataDir)) {
            refdataDir <- realpath(refdataDir) ## nocov
        }
        sampleDirs <- .sampleDirs(dir = dir, filtered = filtered)
        lanes <- detectLanes(sampleDirs)
        assert(isInt(lanes) || identical(lanes, integer()))
        ## Sample metadata -----------------------------------------------------
        allSamples <- TRUE
        sampleData <- NULL
        ## Get the sample data.
        if (isString(sampleMetadataFile)) {
            ## Normalize path of local file.
            if (file.exists(sampleMetadataFile)) {
                sampleMetadataFile <- realpath(sampleMetadataFile)
            }
            ## Note that URL input is also supported here.
            sampleData <- importSampleData(
                file = sampleMetadataFile,
                lanes = lanes,
                pipeline = "cellranger"
            )
            assert(isSubset(rownames(sampleData), names(sampleDirs)))
            sampleIds <- rownames(sampleData)
        } else {
            sampleIds <- names(sampleDirs)
        }
        ## Subset the sample directories, if necessary.
        if (is.character(samples) || is.character(censorSamples)) {
            if (is.character(samples)) {
                samples <- makeNames(samples)
                assert(isSubset(samples, sampleIds))
                sampleIds <- samples
            }
            if (is.character(censorSamples)) {
                censorSamples <- makeNames(censorSamples)
                assert(isSubset(censorSamples, sampleIds))
                sampleIds <- setdiff(sampleIds, censorSamples)
            }
            assert(
                isCharacter(sampleIds),
                isSubset(sampleIds, names(sampleDirs))
            )
        }
        assert(
            hasLength(sampleIds),
            isSubset(sampleIds, names(sampleDirs)),
            validNames(sampleIds)
        )
        if (length(sampleIds) < length(sampleDirs)) {
            sampleDirs <- sampleDirs[sampleIds]
            txt("Loading a subset of samples:")
            ul(basename(sampleDirs))
            ## Subset the user-defined sample metadata to match, if necessary.
            if (!is.null(sampleData)) {
                keep <- rownames(sampleData) %in% sampleIds
                sampleData <- sampleData[keep, , drop = FALSE]
            }
            allSamples <- FALSE
        }
        ## Assays (counts) -----------------------------------------------------
        matrixFiles <- .matrixFiles(
            sampleDirs = sampleDirs,
            filtered = filtered
        )
        ## Get the pipeline from the matrix file attributes.
        pipeline <- attr(matrixFiles, "pipeline")
        assert(isString(pipeline) || identical(pipeline, NA_character_))
        attr(matrixFiles, "pipeline") <- NULL
        counts <- .importCounts(matrixFiles)
        assert(hasValidDimnames(counts))
        ## Row data (genes/transcripts) ----------------------------------------
        refJson <- NULL
        ## Prepare gene annotations as GRanges.
        if (isADir(refdataDir)) {
            ## nocov start
            alertInfo(sprintf(
                fmt = paste0(
                    "Using 10X Genomics reference data ",
                    "for feature annotations: %s"
                ),
                basename(refdataDir)
            ))
            ## JSON data.
            refJsonFile <- file.path(refdataDir, "reference.json")
            assert(isAFile(refJsonFile))
            refJson <- import(refJsonFile)
            ## Get the genome build from JSON metadata.
            genomeBuild <- unlist(refJson[["genomes"]])
            assert(isString(genomeBuild))
            ## Get the Ensembl release version from JSON metadata.
            ## e.g. "Homo_sapiens.GRCh38.93.filtered.gtf"
            ensemblRelease <-
                as.integer(strsplit(
                    x = refJson[["input_gtf_files"]][[1L]],
                    split = ".",
                    fixed = TRUE
                )[[1L]][[3L]])
            assert(isInt(ensemblRelease))
            ## Convert the GTF file to GRanges.
            gffFile <- file.path(refdataDir, "genes", "genes.gtf")
            assert(isString(gffFile))
            rowRanges <- makeGRangesFromGff(gffFile)
            ## nocov end
        } else if (isString(gffFile)) {
            ## This step is necessary for generating v2 working example. Note
            ## that this works with a remote URL.
            rowRanges <- makeGRangesFromGff(
                file = gffFile,
                level = "genes",
                ignoreVersion = TRUE
            )
        } else if (isString(organism)) {
            ## Cell Ranger uses Ensembl refdata internally. Here we're fetching
            ## the annotations with AnnotationHub rather than pulling from the
            ## GTF file in the refdata directory. It will also drop genes that
            ## are now dead in the current Ensembl release. Don't warn about old
            ## Ensembl release version.
            rowRanges <- makeGRangesFromEnsembl(
                organism = organism,
                level = level,
                genomeBuild = genomeBuild,
                release = ensemblRelease,
                ignoreVersion = TRUE
            )
            if (is.null(genomeBuild)) {
                genomeBuild <- metadata(rowRanges)[["genomeBuild"]]
            }
            if (is.null(ensemblRelease)) {
                ensemblRelease <- metadata(rowRanges)[["ensemblRelease"]]
            }
        } else {
            alertWarning(sprintf(
                "Slotting empty ranges into {.fun %s}.",
                "rowRanges"
            ))
            rowRanges <- emptyRanges(rownames(counts))
        }
        assert(is(rowRanges, "GenomicRanges"))
        ## Metrics -------------------------------------------------------------
        ## Note that "molecule_info.h5" file contains additional information
        ## that may be useful for quality control metric calculations.
        aggregation <- NULL
        sampleMetrics <- NULL
        summary <- NULL
        if (.isAggregate(dir)) {
            aggregation <- import(file.path(dir, "outs", "aggregation.csv"))
            aggregation <- as(aggregation, "DataFrame")
            summary <- import(file.path(dir, "outs", "summary.json"))
            summary <- as(summary, "SimpleList")
        } else if (!.isMinimalSample(dir)) {
            sampleMetrics <- .importSampleMetrics(sampleDirs)
        }
        ## Column data ---------------------------------------------------------
        colData <- DataFrame(row.names = colnames(counts))
        ## Generate automatic sample metadata, if necessary.
        if (is.null(sampleData)) {
            ## Define the grep pattern to use for sample ID extraction.
            pattern <- "^(.+)_[ACGT]+$"
            if (all(grepl(pattern, colnames(counts)))) {
                match <- strMatch(x = colnames(counts), pattern = pattern)
                samples <- unique(match[, 2L, drop = TRUE])
            } else if (hasLength(sampleDirs, n = 1L)) {
                samples <- names(sampleDirs)
            }
            sampleData <- minimalSampleData(samples)
        }
        ## Join `sampleData` into cell-level `colData`.
        if (identical(nrow(sampleData), 1L)) {
            colData[["sampleId"]] <- as.factor(rownames(sampleData))
        } else {
            colData[["sampleId"]] <- mapCellsToSamples(
                cells = rownames(colData),
                samples = rownames(sampleData)
            )
        }
        sampleData[["sampleId"]] <- as.factor(rownames(sampleData))
        ## Need to ensure the `sampleId` factor levels match up, otherwise we'll
        ## get a warning during the `leftJoin()` call below.
        assert(areSetEqual(
            x = levels(colData[["sampleId"]]),
            y = levels(sampleData[["sampleId"]])
        ))
        levels(sampleData[["sampleId"]]) <- levels(colData[["sampleId"]])
        colData <- leftJoin(colData, sampleData, by = "sampleId")
        assert(
            is(colData, "DataFrame"),
            hasRownames(colData)
        )
        ## Metadata ------------------------------------------------------------
        interestingGroups <- camelCase(interestingGroups, strict = TRUE)
        assert(isSubset(interestingGroups, colnames(sampleData)))
        metadata <- list(
            "aggregation" = aggregation,
            "allSamples" = allSamples,
            "call" = standardizeCall(),
            "dir" = dir,
            "ensemblRelease" = as.integer(ensemblRelease),
            "genomeBuild" = as.character(genomeBuild),
            "gffFile" = as.character(gffFile),
            "interestingGroups" = interestingGroups,
            "lanes" = lanes,
            "level" = level,
            "matrixFiles" = matrixFiles,
            "organism" = as.character(organism),
            "packageVersion" = .pkgVersion,
            "pipeline" = pipeline,
            "refJson" = as.list(refJson),
            "refdataDir" = as.character(refdataDir),
            "sampleDirs" = sampleDirs,
            "sampleMetadataFile" = as.character(sampleMetadataFile),
            "sampleMetrics" = sampleMetrics,
            "summary" = summary,
            "umiType" = "chromium"
        )
        ## SingleCellExperiment ------------------------------------------------
        object <- makeSingleCellExperiment(
            assays = SimpleList(counts = counts),
            rowRanges = rowRanges,
            colData = colData,
            metadata = metadata,
            transgeneNames = transgeneNames
        )
        ## Return --------------------------------------------------------------
        ## Always prefilter, removing very low quality cells and/or genes.
        object <- calculateMetrics(object = object, prefilter = TRUE)
        object <- new(Class = "CellRanger", object)
        alertSuccess("Chromium single-cell RNA-seq run imported successfully.")
        object
    }
