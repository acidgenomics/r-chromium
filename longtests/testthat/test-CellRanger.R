test_that("v2 : HDF5", {
    dir <- dirs[["pbmc_v2_hdf5"]]
    object <- CellRanger(dir, filtered = FALSE)
    expect_s4_class(object, "CellRanger")
    expect_error(
        object = CellRanger(dir, filtered = TRUE),
        regexp = "Failed to detect samples"
    )
})

test_that("v2 : MTX", {
    dir <- dirs[["pbmc_v2_mtx"]]
    object <- CellRanger(dir, filtered = TRUE)
    expect_s4_class(object, "CellRanger")
    object <- CellRanger(dir, filtered = FALSE)
    expect_s4_class(object, "CellRanger")
})

test_that("v3 : HDF5", {
    dir <- dirs[["pbmc_v3_hdf5"]]
    object <- CellRanger(dir, filtered = TRUE)
    expect_s4_class(object, "CellRanger")
    object <- CellRanger(dir, filtered = FALSE)
    expect_s4_class(object, "CellRanger")
})

test_that("v3 : MTX", {
    dir <- dirs[["pbmc_v2_mtx"]]
    object <- CellRanger(dir, filtered = TRUE)
    expect_s4_class(object, "CellRanger")
    object <- CellRanger(dir, filtered = FALSE)
    expect_s4_class(object, "CellRanger")
})
