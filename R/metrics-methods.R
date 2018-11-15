#' @rdname metrics
#' @export
setMethod(
    f = "metrics",
    signature = signature("CellRanger"),
    definition = metrics.SingleCellExperiment
)
