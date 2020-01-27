## CPI test data examples
##
## Single sample:
## > dir <- file.path(
## >     "",
## >     "mnt",
## >     "azbioinfoseq03",
## >     "scRNA-seq",
## >     "2019_06_CBPHAT_LNCaP_MCF7_CPI1612_scRNAseq",
## >     "cellranger",
## >     "MCF7_50nM_CPI1612"
## > x <- CellRanger(dir)
##
## Aggregation:
## > dir <- file.path(
## >     "",
## >     "mnt",
## >     "azbioinfoseq03",
## >     "scRNA-seq",
## >     "2019_06_CBPHAT_LNCaP_MCF7_CPI1612_scRNAseq",
## >     "cellranger",
## >     "MCF7"
## > x <- CellRanger(dir)



#' @inherit CellRanger-class title description
#' @note Currently supports loading of a single genome.
#' @note Updated 2020-01-26.
#' @export
#'
#' @details
#' Read [10x Genomics Cell Ranger](https://www.10xgenomics.com/software/) output
#' for a Chromium data set into a `SingleCellExperiment` object.
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
#' @inheritParams acidroxygen::params
#' @param dir `character(1)`.
#'   Directory path to Cell Ranger output.
#' @param filtered `logical(1)`.
#'   Use filtered (recommended) or raw counts. Note that raw counts still
#'   contain only whitelisted cellular barcodes.
#' @param refdataDir `character(1)` or `NULL`.
#'   Directory path to Cell Ranger reference annotation data.
#'
#' @return `CellRanger`.
#'
#' @examples
#' dir <- system.file("extdata/cellranger_v3", package = "Chromium")
#' x <- CellRanger(dir)
#' print(x)
CellRanger <- function(  # nolint
    dir,
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
    spikeNames = NULL,
    interestingGroups = "sampleName",
    BPPARAM = BiocParallel::bpparam()  # nolint
) {
    assert(
        isADirectory(dir),
        isFlag(filtered),
        isString(organism, nullOK = TRUE),
        isInt(ensemblRelease, nullOK = TRUE),
        isString(genomeBuild, nullOK = TRUE),
        isString(gffFile, nullOK = TRUE),
        isADirectory(refdataDir, nullOK = TRUE),
        isAny(samples, classes = c("character", "NULL")),
        isAny(censorSamples, classes = c("character", "NULL")),
        isAFile(sampleMetadataFile, nullOK = TRUE),
        isCharacter(transgeneNames, nullOK = TRUE),
        isCharacter(spikeNames, nullOK = TRUE),
        isCharacter(interestingGroups),
        identical(attr(class(BPPARAM), "package"), "BiocParallel")
    )
    
    cli_h1("Importing Chromium single-cell RNA-seq run")
    
    ## Run info ----------------------------------------------------------------
    level <- "genes"
    dir <- realpath(dir)
    if (isADirectory(refdataDir)) {
        refdataDir <- realpath(refdataDir)  ## nocov
    }
    sampleDirs <- .sampleDirs(dir)
    lanes <- detectLanes(sampleDirs)
    assert(isInt(lanes) || identical(lanes, integer()))

    ## Samples -----------------------------------------------------------------
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
        sampleIDs <- rownames(sampleData)
    } else {
        sampleIDs <- names(sampleDirs)
    }
    ## Subset the sample directories, if necessary.
    if (is.character(samples) || is.character(censorSamples)) {
        if (is.character(samples)) {
            samples <- makeNames(samples)
            assert(isSubset(samples, sampleIDs))
            sampleIDs <- samples
        }
        if (is.character(censorSamples)) {
            censorSamples <- makeNames(censorSamples)
            assert(isSubset(censorSamples, sampleIDs))
            sampleIDs <- setdiff(sampleIDs, censorSamples)
        }
        assert(
            isCharacter(sampleIDs),
            isSubset(sampleIDs, names(sampleDirs))
        )
    }
    assert(
        hasLength(sampleIDs),
        isSubset(sampleIDs, names(sampleDirs)),
        validNames(sampleIDs)
    )
    if (length(sampleIDs) < length(sampleDirs)) {
        sampleDirs <- sampleDirs[sampleIDs]
        message(sprintf(
            fmt = "Loading a subset of samples:\n%s",
            str_trunc(toString(basename(sampleDirs)), width = 80L)
        ))
        ## Subset the user-defined sample metadata to match, if necessary.
        if (!is.null(sampleData)) {
            keep <- rownames(sampleData) %in% sampleIDs
            sampleData <- sampleData[keep, , drop = FALSE]
        }
        allSamples <- FALSE
    }

    ## Assays ------------------------------------------------------------------
    matrixFiles <- .matrixFiles(
        sampleDirs = sampleDirs,
        filtered = filtered,
        BPPARAM = BPPARAM
    )
    ## Get the pipeline from the matrix file attributes.
    pipeline <- attr(matrixFiles, "pipeline")
    assert(isString(pipeline) || identical(pipeline, NA_character_))
    attr(matrixFiles, "pipeline") <- NULL
    counts <- .importCounts(
        matrixFiles = matrixFiles,
        BPPARAM = BPPARAM
    )
    assert(hasValidDimnames(counts))

    ## Row data ----------------------------------------------------------------
    refJSON <- NULL
    ## Prepare gene annotations as GRanges.
    if (isADirectory(refdataDir)) {
        ## nocov start
        message(sprintf(
            fmt = paste0(
                "Using 10X Genomics reference data ",
                "for feature annotations: %s"
            ),
            basename(refdataDir)
        ))
        ## JSON data.
        refJSONFile <- file.path(refdataDir, "reference.json")
        assert(isAFile(refJSONFile))
        refJSON <- import(refJSONFile)
        ## Get the genome build from JSON metadata.
        genomeBuild <- unlist(refJSON[["genomes"]])
        assert(isString(genomeBuild))
        ## Get the Ensembl release version from JSON metadata.
        ## e.g. "Homo_sapiens.GRCh38.93.filtered.gtf"
        ensemblRelease <- refJSON[["input_gtf_files"]][[1L]]
        ensemblRelease <- str_split(
            string = ensemblRelease,
            pattern = "\\.",
            simplify = TRUE
        )
        ensemblRelease <- as.integer(ensemblRelease[1L, 3L])
        assert(isInt(ensemblRelease))
        ## Convert the GTF file to GRanges.
        gffFile <- file.path(refdataDir, "genes", "genes.gtf")
        assert(isString(gffFile))
        rowRanges <- makeGRangesFromGFF(gffFile)
        ## nocov end
    } else if (isString(gffFile)) {
        ## This step is necessary for generating v2 working example.
        ## Note that this works with a remote URL.
        rowRanges <- makeGRangesFromGFF(gffFile, level = "genes")  # nocov
    } else if (isString(organism)) {
        ## Cell Ranger uses Ensembl refdata internally. Here we're fetching the
        ## annotations with AnnotationHub rather than pulling from the GTF file
        ## in the refdata directory. It will also drop genes that are now dead
        ## in the current Ensembl release. Don't warn about old Ensembl release
        ## version.
        message("Using makeGRangesFromEnsembl() for annotations.")
        rowRanges <- makeGRangesFromEnsembl(
            organism = organism,
            level = level,
            genomeBuild = genomeBuild,
            release = ensemblRelease
        )
        if (is.null(genomeBuild)) {
            genomeBuild <- metadata(rowRanges)[["genomeBuild"]]
        }
        if (is.null(ensemblRelease)) {
            ensemblRelease <- metadata(rowRanges)[["ensemblRelease"]]
        }
    } else {
        message("Unknown organism. Skipping annotations.")
        rowRanges <- emptyRanges(rownames(counts))
    }
    assert(is(rowRanges, "GRanges"))

    ## Metrics -----------------------------------------------------------------
    ## Note that "molecule_info.h5" file contains additional information that
    ## may be useful for quality control metric calculations.
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

    ## Column data -------------------------------------------------------------
    colData <- DataFrame(row.names = colnames(counts))
    ## Generate automatic sample metadata, if necessary.
    if (is.null(sampleData)) {
        ## Define the grep pattern to use for sample ID extraction.
        pattern <- "^(.+)_[ACGT]+$"
        if (all(grepl(pattern, colnames(counts)))) {
            match <- str_match(
                string = colnames(counts),
                pattern = pattern
            )
            samples <- unique(match[, 2L, drop = TRUE])
        } else if (hasLength(sampleDirs, n = 1L)) {
            samples <- names(sampleDirs)
        }
        sampleData <- minimalSampleData(samples)
    }
    ## Join `sampleData` into cell-level `colData`.
    if (identical(nrow(sampleData), 1L)) {
        colData[["sampleID"]] <- as.factor(rownames(sampleData))
    } else {
        colData[["sampleID"]] <- mapCellsToSamples(
            cells = rownames(colData),
            samples = rownames(sampleData)
        )
    }
    sampleData[["sampleID"]] <- as.factor(rownames(sampleData))
    ## Need to ensure the `sampleID` factor levels match up, otherwise we'll get
    ## a warning during the `leftJoin()` call below.
    assert(areSetEqual(
        x = levels(colData[["sampleID"]]),
        y = levels(sampleData[["sampleID"]])
    ))
    levels(sampleData[["sampleID"]]) <- levels(colData[["sampleID"]])
    colData <- leftJoin(colData, sampleData, by = "sampleID")
    assert(
        is(colData, "DataFrame"),
        hasRownames(colData)
    )

    ## Metadata ----------------------------------------------------------------
    interestingGroups <- camelCase(interestingGroups)
    assert(isSubset(interestingGroups, colnames(sampleData)))
    metadata <- list(
        aggregation = aggregation,
        allSamples = allSamples,
        call = standardizeCall(),
        dir = dir,
        ensemblRelease = as.integer(ensemblRelease),
        genomeBuild = as.character(genomeBuild),
        gffFile = as.character(gffFile),
        interestingGroups = interestingGroups,
        lanes = lanes,
        level = level,
        matrixFiles = matrixFiles,
        organism = as.character(organism),
        pipeline = pipeline,
        refJSON = as.list(refJSON),
        refdataDir = as.character(refdataDir),
        sampleDirs = sampleDirs,
        sampleMetadataFile = as.character(sampleMetadataFile),
        sampleMetrics = sampleMetrics,
        summary = summary,
        umiType = "chromium",
        version = .version
    )

    ## SingleCellExperiment ----------------------------------------------------
    object <- makeSingleCellExperiment(
        assays = SimpleList(counts = counts),
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        transgeneNames = transgeneNames,
        spikeNames = spikeNames
    )

    ## Return ------------------------------------------------------------------
    ## Always prefilter, removing very low quality cells and/or genes.
    object <- calculateMetrics(object = object, prefilter = TRUE)
    object <- new(Class = "CellRanger", object)
    cli_alert_success("Chromium single-cell RNA-seq run imported successfully.")
    object
}
