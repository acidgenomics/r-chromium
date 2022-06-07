.pkgName <- packageName()
.pkgVersion <- packageVersion(.pkgName)

#' Chromium test data URL
#' @keywords internal
#' @export
#' @examples
#' ChromiumTestsURL
ChromiumTestsURL <- # nolint
    paste0(
        "https://r.acidgenomics.com/testdata/chromium/",
        "v", .pkgVersion$major, ".", .pkgVersion$minor # nolint
    )
