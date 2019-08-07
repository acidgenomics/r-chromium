## FIXME Add support for single sample mode.
## FIXME Consider moving calculate metrics to basejump?
## Make calculateMetrics an S4 generic...
## FIXME Consider using HDF5Array here instead
## FIXME Check for single genome.
## FIXME Need to rework approach to aggregate samples with "_1", "_2" suffix.



#' @inherit CellRanger-class title description
#' @note Updated 2019-08-07.
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
#' \preformatted{
#' | <dir>/
#' |-- <sampleName>/
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
#' )
#' }
#'
#' Cell Ranger v2 output:
#'
#' \preformatted{
#' | <dir>/
#' |-- <sampleName>/
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
#' )
#' }
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
#' @author Michael Steinbaugh
#' @note Updated 2019-08-07.
#' @export
#' 
#' @inheritParams acidroxygen::params
#' @param dir `character(1)`.
#'   Path to Cell Ranger output directory (final upload). This directory path
#'   must contain `filtered_gene_bc_matrices*` as a child directory.
#' @param filtered `logical(1)`.
#'   Use filtered (recommended) or raw counts. Note that raw counts still
#'   contain only whitelisted cellular barcodes.
#' @param refdataDir `character(1)` or `NULL`.
#'   Directory path to Cell Ranger reference annotation data.
#'
#' @return `CellRanger`.
#'
#' @examples
#' dir <- system.file("extdata/cellranger", package = "Chromium")
#' x <- CellRanger(dir)
#' print(x)
CellRanger <- function(
    dir,
    filtered = TRUE,
    organism = NULL,
    ensemblRelease = NULL,
    genomeBuild = NULL,
    refdataDir = NULL,
    samples = NULL,
    censorSamples = NULL,
    sampleMetadataFile = NULL,
    transgeneNames = NULL,
    spikeNames = NULL
) {
    assert(
        isADirectory(dir),
        isFlag(filtered),
        isString(organism, nullOK = TRUE),
        isInt(ensemblRelease, nullOK = TRUE),
        isString(genomeBuild, nullOK = TRUE),
        isADirectory(refdataDir, nullOK = TRUE),
        isAny(samples, classes = c("character", "NULL")),
        isAny(censorSamples, classes = c("character", "NULL")),
        isAFile(sampleMetadataFile, nullOK = TRUE),
        isCharacter(transgeneNames, nullOK = TRUE),
        isCharacter(spikeNames, nullOK = TRUE)
    )
    level <- "genes"
    
    ## Directory paths ---------------------------------------------------------
    dir <- realpath(dir)
    if (isADirectory(refdataDir)) {
        refdataDir <- realpath(refdataDir)
    }
    sampleDirs <- .sampleDirs(dir)
    
    ## Sequencing lanes --------------------------------------------------------
    lanes <- detectLanes(sampleDirs)
    assert(
        isInt(lanes) ||
            identical(lanes, integer())
    )
    
    ## Samples -----------------------------------------------------------------
    allSamples <- TRUE
    sampleData <- NULL
    
    ## Get the sample data.
    if (isString(sampleMetadataFile)) {
        ## Normalize path of local file.
        if (file.exists(sampleMetadataFile)) {
            sampleMetadataFile <- realpath(sampleMetadataFile)
        }
        ## Note that `readSampleData()` also supports URLs.
        sampleData <- readSampleData(file = sampleMetadataFile, lanes = lanes)
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
    counts <- .importCounts(sampleDirs = sampleDirs, filtered = filtered)
    assays <- SimpleList(counts = counts)

    ## Row data ----------------------------------------------------------------
    refJSON <- NULL
    
    ## Prepare gene annotations as GRanges.
    if (isADirectory(refdataDir)) {
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
        ensemblRelease <- refJSON %>%
            .[["input_gtf_files"]] %>%
            .[[1L]] %>%
            str_split("\\.", simplify = TRUE) %>%
            .[1L, 3L] %>%
            as.integer()
        assert(isInt(ensemblRelease))
        ## Convert the GTF file to GRanges.
        gffFile <- file.path(refdataDir, "genes", "genes.gtf")
        assert(isString(gffFile))
        rowRanges <- makeGRangesFromGFF(gffFile)
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
    
    ## Column data -------------------------------------------------------------
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
        } else if (length(sampleFiles) == 1L) {
            samples <- names(sampleFiles)
        }
        sampleData <- minimalSampleData(samples)
    }
    
    ## Always prefilter, removing very low quality cells with no UMIs or genes.
    ## FIXME Move this code to basejump to avoid bcbioSingleCell dependency.
    
    ## FIXME Seeing this warning pop up:
    ## Warning in NSBS(i, x, exact = exact, strict.upper.bound = !allow.append,
    ## : subscript is an array, passing it thru as.vector() first
    ## Calls: <Anonymous> ... extractROWS -> normalizeSingleBracketSubscript ->
    ## NSBS -> NSBS
    ## FIXME This error is only happening when we pass in rowRanges.
    
    colData <- bcbioSingleCell::calculateMetrics(
        counts = counts,
        rowRanges = rowRanges,
        prefilter = TRUE
    )
    
    ## Subset the counts to match the prefiltered metrics.
    assert(isSubset(rownames(colData), colnames(counts)))
    counts <- counts[, rownames(colData), drop = FALSE]
    
    ## Join `sampleData` into cell-level `colData`.
    if (length(nrow(sampleData)) == 1L) {
        colData[["sampleID"]] <- as.factor(rownames(sampleData))
    } else {
        colData[["sampleID"]] <- mapCellsToSamples(
            cells = rownames(colData),
            samples = rownames(sampleData)
        )
    }
    
    ## Metadata ----------------------------------------------------------------
    metadata <- list(
        version = .version,
        level = level,
        dir = dir,
        sampleMetadataFile = as.character(sampleMetadataFile),
        organism = organism,
        genomeBuild = as.character(genomeBuild),
        ensemblRelease = as.integer(ensemblRelease),
        allSamples = allSamples,
        lanes = lanes,
        pipeline = "cellranger",
        umiType = "chromium",
        ## cellranger-specific -------------------------------------------------
        refdataDir = refdataDir,
        refJSON = refJSON,
        call = match.call()
    )
    
    ## Return ------------------------------------------------------------------
    `new,CellRanger`(
        assays = assays = assays,
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        transgeneNames = transgeneNames,
        spikeNames = spikeNames
    )
}



## Updated 2019-08-01.
`new,CellRanger` <-  # nolint
    function(...) {
        new(Class = "CellRanger", makeSingleCellExperiment(...))
    }
