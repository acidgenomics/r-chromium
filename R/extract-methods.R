#' @name extract
#' @author Michael Steinbaugh
#' @inherit base::Extract title params references
#' @note Updated 2022-06-07.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @description Extract genes by row and cells by column.
#'
#' @details
#' Refer to [`cell2sample()`][basejump::cell2sample] and
#' [`selectSamples()`][basejump::selectSamples] if sample-level extraction is
#' desired. Note that `sampleId` is slotted into
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
#' data(pbmc_v3)
#'
#' ## CellRanger ====
#' object <- pbmc_v3
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



## This approach is adapted from bcbioSingleCell method.
## Updated 2022-06-07.
`extract,CellRanger` <- # nolint
    function(x, i, j, ..., drop = FALSE) {
        validObject(x)
        ## Genes (rows).
        if (missing(i)) {
            i <- seq_len(nrow(x))
        }
        ## Cells (columns).
        if (missing(j)) {
            j <- seq_len(ncol(x))
        }
        ## Determine whether we should stash subset in metadata.
        if (identical(x = dim(x), y = c(length(i), length(j)))) {
            subset <- FALSE
        } else {
            subset <- TRUE
        }
        ## Subset using SCE method.
        sce <- as(x, "SingleCellExperiment")
        sce <- sce[i, j, drop = drop]
        ## Early return original object, if unmodified.
        if (identical(assay(sce), assay(x))) {
            return(x)
        }
        ## Metadata ------------------------------------------------------------
        metadata <- metadata(sce)
        if (isTRUE(subset)) {
            metadata[["filterGenes"]] <- NULL
            metadata[["subset"]] <- TRUE
        }
        metadata <- Filter(f = Negate(is.null), x = metadata)
        metadata(sce) <- metadata
        ## Return --------------------------------------------------------------
        sce <- droplevels2(sce)
        new(Class = "CellRanger", sce)
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
