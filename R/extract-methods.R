#' @rdname extract
#' @export
setMethod(
    "[",
    signature(
        x = "CellRanger",
        i = "ANY",
        j = "ANY",
        drop = "ANY"
    ),
    definition = extract.SingleCellExperiment
)
