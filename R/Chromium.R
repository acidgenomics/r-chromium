# FIXME Break these back out into separate files.

# FIXME Either make `sampleName` required or strip it from minimal examples.

# FIXME Check to see if we can import tx2gene.csv

# FIXME Can we parse the Cell Ranger `runDate` from the refData YAML?

# FIXME Allow this function to work if the user points at dir containing matrix.

# FIXME Add documentation about simple mode.

# FIXME Consolidate with bcbioSingleCell
# Undocumented arguments in documentation object 'Chromium'
# ‘sampleMetadataFile’ ‘genomeBuild’ ‘gffFile’ ‘transgeneNames’
# ‘spikeNames’



#' @rdname Chromium-class
#' @export
#'
#' @details
#' Read [10x Genomics Cell Ranger](https://www.10xgenomics.com/software/) output
#' for a Chromium data set into a `SingleCellExperiment` object.
#'
#' @section Directory structure:
#' 
#' Cell Ranger can vary in its output directory structure, but we're requiring a
#' single, consistent directory structure for all datasets, even those that only
#' contain a single sample:
#'
#' \preformatted{
#' file.path(
#'     "<dir>",
#'     "<sampleName>",
#'     "outs",
#'     "filtered_gene_bc_matrices*",
#'     "<genomeBuild>",
#'     "matrix.mtx"
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
#' slotted in `rowRanges`. Otherwise, the function will attempt to use the
#' most current annotations available from Ensembl, and some gene IDs may not
#' match, due to deprecation in the current Ensembl release.
#'
#' @author Michael Steinbaugh
#' @export
#' @inheritParams basejump::params
#'
#' @param dir `character(1)`.
#'   Path to Cell Ranger output directory (final upload). This directory path
#'   must contain `filtered_gene_bc_matrices*` as a child directory.
#' @param filtered `logical(1)`.
#'   Use filtered (recommended) or raw counts. Note that raw counts still
#'   contain only whitelisted cellular barcodes.
#' @param format `character(1)`.
#'   Output format, either MatrixMarket ("`mtx`") or HDF5 ("`hdf5`").
#' @param refdataDir `character(1)` or `NULL`.
#'   Directory path to Cell Ranger reference annotation data.
#'
#' @return `Chromium`.
#'
#' @examples
#' dir <- system.file("extdata/cellranger", package = "bcbioSingleCell")
#' x <- Chromium(dir)
#' print(x)
Chromium <- function(  # nolint
    dir,
    format = c("mtx", "hdf5"),
    filtered = TRUE,
    sampleMetadataFile = NULL,
    organism = NULL,
    ensemblRelease = NULL,
    genomeBuild = NULL,
    refdataDir = NULL,
    gffFile = NULL,
    transgeneNames = NULL,
    spikeNames = NULL,
    interestingGroups = "sampleName"
) {
    assert(
        isADirectory(dir),
        isFlag(filtered),
        isAFile(sampleMetadataFile, nullOK = TRUE),
        isString(organism, nullOK = TRUE),
        isInt(ensemblRelease, nullOK = TRUE),
        isString(genomeBuild, nullOK = TRUE),
        isADirectory(refdataDir, nullOK = TRUE),
        isAFile(gffFile, nullOK = TRUE),
        isCharacter(transgeneNames, nullOK = TRUE),
        isCharacter(spikeNames, nullOK = TRUE),
        isCharacter(interestingGroups)
    )
    dir <- realpath(dir)
    format <- match.arg(format)
    if (isADirectory(refdataDir)) {
        refdataDir <- realpath(refdataDir)
    }
    
    pipeline <- "cellranger"
    level <- "genes"
    umiType <- "chromium"
    
    # Sample files -------------------------------------------------------------
    sampleFiles <- .sampleFiles(dir = dir, format = format, filtered = filtered)
    
    # Sequencing lanes ---------------------------------------------------------
    lanes <- detectLanes(sampleFiles)
    
    # Sample metadata ----------------------------------------------------------
    allSamples <- TRUE
    sampleData <- NULL
    
    if (isAFile(sampleMetadataFile)) {
        sampleData <- readSampleData(sampleMetadataFile)
        # Allow sample selection by with this file.
        if (nrow(sampleData) < length(sampleFiles)) {
            message("Loading a subset of samples, defined by the metadata.")
            allSamples <- FALSE
            sampleFiles <- sampleFiles[rownames(sampleData)]
            message(paste(length(sampleFiles), "samples matched by metadata."))
        }
    }
    
    # Counts -------------------------------------------------------------------
    counts <- .import(sampleFiles)
    
    # Row data -----------------------------------------------------------------
    refJSON <- NULL
    
    # Prepare gene annotations as GRanges.
    if (isADirectory(refdataDir)) {
        message("Using 10X Genomics reference data for annotations.")
        message(paste("refdataDir:", refdataDir))
        # JSON data.
        refJSONFile <- file.path(refdataDir, "reference.json")
        assert(allAreFiles(refJSONFile))
        refJSON <- import(refJSONFile)
        # Get the genome build from JSON metadata.
        genomeBuild <- unlist(refJSON[["genomes"]])
        assert(isString(genomeBuild))
        # Convert the GTF file to GRanges.
        gffFile <- file.path(refdataDir, "genes", "genes.gtf")
        assert(isString(gffFile))
        rowRanges <- makeGRangesFromGFF(gffFile)
        # Get the Ensembl version from the GTF file name.
        # Example: "Homo_sapiens.GRCh37.82.filtered.gtf"
        ensemblRelease <- gffFile %>%
            str_split("\\.", simplify = TRUE) %>%
            .[1L, 3L] %>%
            as.integer()
    } else if (isString(gffFile)) {
        # Note that this works with a remote URL.
        rowRanges <- makeGRangesFromGFF(gffFile, level = "genes")
    } else if (isString(organism)) {
        # Cell Ranger uses Ensembl refdata internally. Here we're fetching the
        # annotations with AnnotationHub rather than pulling from the GTF file
        # in the refdata directory. It will also drop genes that are now dead in
        # the current Ensembl release. Don't warn about old Ensembl release
        # version.
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
    
    # Column data --------------------------------------------------------------
    # Automatic sample metadata.
    if (is.null(sampleData)) {
        # Define the grep pattern to use for sample ID extraction.
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
    
    # Always prefilter, removing very low quality cells with no UMIs or genes.
    colData <- calculateMetrics(
        counts = counts,
        rowRanges = rowRanges,
        prefilter = TRUE
    )
    
    # Subset the counts to match the prefiltered metrics.
    assert(isSubset(rownames(colData), colnames(counts)))
    counts <- counts[, rownames(colData), drop = FALSE]
    
    # Join `sampleData` into cell-level `colData`.
    if (length(nrow(sampleData)) == 1L) {
        colData[["sampleID"]] <- as.factor(rownames(sampleData))
    } else {
        colData[["sampleID"]] <- mapCellsToSamples(
            cells = rownames(colData),
            samples = rownames(sampleData)
        )
    }
    
    # Metadata -----------------------------------------------------------------
    # Interesting groups.
    interestingGroups <- camel(interestingGroups)
    assert(isSubset(interestingGroups, colnames(sampleData)))
    
    metadata <- list(
        version = packageVersion,
        pipeline = pipeline,
        level = level,
        dir = dir,
        sampleMetadataFile = as.character(sampleMetadataFile),
        interestingGroups = interestingGroups,
        organism = organism,
        genomeBuild = as.character(genomeBuild),
        ensemblRelease = as.integer(ensemblRelease),
        umiType = umiType,
        allSamples = allSamples,
        lanes = lanes,
        # cellranger-specific --------------------------------------------------
        refdataDir = refdataDir,
        refJSON = refJSON,
        call = match.call()
    )
    
    # Return -------------------------------------------------------------------
    .new.Chromium(
        assays = list(counts = counts),
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        transgeneNames = transgeneNames,
        spikeNames = spikeNames
    )
}



.new.Chromium <-  # nolint
    function(...) {
        new(Class = "Chromium", makeSingleCellExperiment(...))
    }
