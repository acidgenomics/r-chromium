context("extract")

test_that("CellRanger", {
    object <- pbmc_v3
    object <- object[1L, 1L]
    expect_s4_class(object, "CellRanger")
})
