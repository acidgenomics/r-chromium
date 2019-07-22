context("Chromium")

dir <- system.file("extdata/cellranger", package = "Chromium")

## Minimal mode, with no metadata or annotations.
## This is fast but doesn't slot a lot of useful info.
test_that("Minimal metadata", {
    x <- Chromium(dir = dir)
    expect_s4_class(x, "Chromium")
})

## Automatic organism annotations from AnnotationHub.
test_that("AnnotationHub metadata", {
    x <- Chromium(
        dir = dir,
        organism = "Homo sapiens",
        ensemblRelease = 87L
    )
    expect_s4_class(x, "Chromium")
})
