globalVariables(".")



.version <- packageVersion("Chromium")

#' Chromium test data URL
#' @keywords internal
#' @export
#' @examples
#' ChromiumTestsURL
ChromiumTestsURL <-  # nolint
    paste0(
        "https://tests.acidgenomics.com/Chromium/",
        "v", .version$major, ".", .version$minor  # nolint
    )



## Trailing number is to match cellranger output.
barcodePattern <- ")_([ACGT_]{6,})(_[0-9]+)?$"
