#' @name extract
#' @author Michael Steinbaugh
#' @inherit base::Extract title params references
#' @inheritParams basejump::params
#'
#' @description Extract genes by row and cells by column.
#'
#' @details
#' Refer to [`cell2sample()`][basejump::cell2sample] and
#' [`selectSamples()`][basejump::selectSamples] if sample-level extraction is
#' desired. Note that `sampleID` is slotted into
#' [`colData()`][SummarizedExperiment::colData] and defines the cell-to-sample
#' mappings.
#'
#' Unfiltered cellular barcode distributions for the entire dataset, including
#' cells not kept in the matrix will be dropped in favor of the `nCount` column
#' of `colData`.
#'
#' @return `Chromium`.
#'
#' @examples
#' data(pbmc)
#'
#' cells <- head(colnames(pbmc), 100L)
#' head(cells)
#' genes <- head(rownames(pbmc), 100L)
#' head(genes)
#'
#' ## Subset by cell identifiers.
#' pbmc[, cells]
#'
#' ## Subset by genes.
#' pbmc[genes, ]
#'
#' ## Subset by both genes and cells.
#' pbmc[genes, cells]
NULL



## Updated 2019-08-01.
`extract,Chromium` <-  # nolint
    function(x, i, j, ..., drop = FALSE) {
        validObject(x)
        
        ## Genes
        if (missing(i)) {
            i <- 1L:nrow(x)
        }
        ## Cells
        if (missing(j)) {
            j <- 1L:ncol(x)
        }
        
        if (identical(
            x = c(length(i), length(j)),
            y = dim(x)
        )) {
            return(x)
        }
        
        ## Subset using SCE method.
        sce <- as(x, "SingleCellExperiment")
        sce <- sce[i, j, drop = drop]
        
        genes <- rownames(sce)
        cells <- colnames(sce)
        
        ## Row data -------------------------------------------------------------
        ## Ensure factors get releveled, if necessary.
        rowRanges <- rowRanges(sce)
        if (
            ncol(mcols(rowRanges)) > 0L &&
            !identical(rownames(sce), rownames(x))
        ) {
            rowRanges <- relevelRowRanges(rowRanges)
        }
        
        ## Column data ----------------------------------------------------------
        ## Ensure factors get releveled, if necessary.
        colData <- colData(sce)
        if (
            ncol(colData) > 0L &&
            !identical(colnames(sce), colnames(x))
        ) {
            colData <- relevelColData(colData)
        }
        
        ## Metadata -------------------------------------------------------------
        metadata <- metadata(sce)
        metadata[["subset"]] <- TRUE
        
        ## Drop unfiltered cellular barcode list.
        metadata[["cellularBarcodes"]] <- NULL
        
        ## filterCells
        filterCells <- metadata[["filterCells"]]
        if (!is.null(filterCells)) {
            filterCells <- intersect(filterCells, cells)
            metadata[["filterCells"]] <- filterCells
        }
        
        ## filterGenes
        filterGenes <- metadata[["filterGenes"]]
        if (!is.null(filterGenes)) {
            filterGenes <- intersect(filterGenes, genes)
            metadata[["filterGenes"]] <- filterGenes
        }
        
        ## Return ---------------------------------------------------------------
        new(
            Class = class(x)[[1L]],
            makeSingleCellExperiment(
                assays = assays(sce),
                rowRanges <- rowRanges(sce),
                colData <- colData,
                metadata = metadata,
                spikeNames = rownames(sce)[isSpike(sce)]
            )
        )
    }



#' @rdname extract
#' @export
setMethod(
    f = "[",
    signature = signature(
        x = "Chromium",
        i = "ANY",
        j = "ANY",
        drop = "ANY"
    ),
    definition = `extract,Chromium`
)
