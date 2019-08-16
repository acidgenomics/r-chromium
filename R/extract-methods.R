#' @name extract
#' @author Michael Steinbaugh
#' @inherit base::Extract title params references
#' @inheritParams acidroxygen::params
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
#' @return `CellRanger`.
#'
#' @examples
#' data(pbmc4k_v2)
#' 
#' ## CellRanger ====
#' object <- pbmc4k_v2
#'
#' cells <- head(colnames(object), 100L)
#' head(cells)
#' genes <- head(rownames(object), 100L)
#' head(genes)
#'
#' ## Subset by cell identifiers.
#' object[, cells]
#'
#' ## Subset by genes.
#' object[genes, ]
#'
#' ## Subset by both genes and cells.
#' object[genes, cells]
NULL



## Updated 2019-08-07.
`extract,CellRanger` <-  # nolint
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
        ## FIXME Simplify
        rowRanges <- rowRanges(sce)
        
        ## Column data ----------------------------------------------------------
        ## FIXME Simplify
        colData <- colData(sce)
        
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
        ## FIXME Simplify this.
        ## FIXME Look at current bcbioSingleCell approach.
        ## FIXME droplevels call.
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
        x = "CellRanger",
        i = "ANY",
        j = "ANY",
        drop = "ANY"
    ),
    definition = `extract,CellRanger`
)
