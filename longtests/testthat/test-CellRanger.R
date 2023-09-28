## FIXME Need to organize these better.
## Extract the tar.gz files.

test_that("v2 : HDF5", {
    dir <- pbmc2
    object <- CellRanger(dir)
    expect_s4_class(object, "CellRanger")
})

## FIXME Rework this under longtests.
test_that("v3 : HDF5", {
    dir <- pbmc3
    object <- CellRanger(dir)
    expect_s4_class(object, "CellRanger")
})
