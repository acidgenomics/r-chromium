## FIXME Rework this under long tests.
test_that("HDF5 v2", {
    dir <- file.path(cacheDir, "pbmc_v2")
    object <- CellRanger(dir)
    expect_s4_class(object, "CellRanger")
})

## FIXME Rework this under longtests.
test_that("v3 HDF5", {
    dir <- file.path(cacheDir, "pbmc_v3")
    object <- CellRanger(dir)
    expect_s4_class(object, "CellRanger")
})
