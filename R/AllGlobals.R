.pkgName <- packageName()
.pkgVersion <- packageVersion(.pkgName)



#' Chromium test data URL
#' @keywords internal
#' @export
#' @examples
#' ChromiumTestsURL
ChromiumTestsURL <-  # nolint
    paste0(
        "https://r.acidgenomics.com/testdata/chromium/",
        "v", .pkgVersion$major, ".", .pkgVersion$minor  # nolint
    )



## Trailing number is to match cellranger output.
barcodePattern <- ")_([ACGT_]{6,})(_[0-9]+)?$"
