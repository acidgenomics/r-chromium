#' Import sample-level metrics
#' @note Updated 2019-08-22.
#' @noRd
.importSampleMetrics <- function(sampleDirs) {
    files <- file.path(sampleDirs, "outs", "metrics_summary.csv")
    list <- lapply(
        X = files,
        FUN = function(file) {
            data <- withCallingHandlers(
                expr = import(file),
                message = function(m) {
                    if (grepl("syntactic", m)) {
                        invokeRestart("muffleMessage")
                    }
                    m
                }
            )
        }
    )
    out <- DataFrame(do.call(what = rbind, args = list))
    out <- camelCase(out)
    rownames(out) <- names(sampleDirs)
    out
}
