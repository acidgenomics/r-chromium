context("Cell Ranger v2")

dir <- system.file("extdata/cellranger_v2", package = "Chromium")

test_that("MTX : Fast mode", {
    x <- CellRanger(dir = dir)
    expect_s4_class(x, "CellRanger")
})

test_that("MTX : User-defined sample metadata", {
    sampleMetadataFile <-
        system.file("extdata/cellranger_v2.csv", package = "Chromium")
    object <- CellRanger(
        dir = dir,
        sampleMetadataFile = sampleMetadataFile
    )
    expect_s4_class(object, "CellRanger")
    expect_identical(
        levels(colData(object)[["sampleName"]]),
        "pbmc4k"
    )
})

test_that("MTX : AnnotationHub", {
    object <- CellRanger(
        dir = dir,
        organism = "Homo sapiens",
        ensemblRelease = 87L
    )
    expect_s4_class(object, "CellRanger")
})

test_that("HDF5", {
    skip_on_appveyor()
    dir <- file.path("cache", "pbmc_v2")
    object <- CellRanger(dir)
    expect_s4_class(object, "CellRanger")
})



context("Cell Ranger v3")

dir <- system.file("extdata/cellranger_v3", package = "Chromium")

test_that("MTX : Fast mode", {
    object <- CellRanger(dir)
    expect_s4_class(object, "CellRanger")
})

test_that("MTX : User-defined sample metadata", {
    sampleMetadataFile <-
        system.file("extdata/cellranger_v3.csv", package = "Chromium")
    object <- CellRanger(
        dir = dir,
        sampleMetadataFile = sampleMetadataFile
    )
    expect_s4_class(object, "CellRanger")
    expect_identical(
        levels(colData(object)[["sampleName"]]),
        "5k_pbmc_protein_v3"
    )
})

test_that("HDF5", {
    skip_on_appveyor()
    dir <- file.path("cache", "pbmc_v3")
    object <- CellRanger(dir)
    expect_s4_class(object, "CellRanger")
})
