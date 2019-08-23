context("Cell Ranger v2")

dir <- system.file("extdata/cellranger_v2", package = "Chromium")

## Minimal mode, with no metadata or annotations.
## This is fast but doesn't slot a lot of useful info.
test_that("Minimal metadata", {
    x <- CellRanger(dir = dir)
    expect_s4_class(x, "CellRanger")
})

## Automatic organism annotations from AnnotationHub.
test_that("AnnotationHub metadata", {
    object <- CellRanger(
        dir = dir,
        organism = "Homo sapiens",
        ensemblRelease = 87L
    )
    expect_s4_class(object, "CellRanger")
})

test_that("HDF5", {
    dir <- file.path("cache", "pbmc_v2")
    object <- CellRanger(dir)
    expect_s4_class(object, "CellRanger")
})



context("Cell Ranger v3")

dir <- system.file("extdata/cellranger_v3", package = "Chromium")

## FIXME
test_that("User-defined sample metadata", {
})

test_that("HDF5", {
    dir <- file.path("cache", "pbmc_v3")
    object <- CellRanger(dir)
    expect_s4_class(object, "CellRanger")
})
