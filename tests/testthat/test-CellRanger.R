lst <- list(
    "v2" = list(
        "dir" = system.file(
            "extdata", "cellranger_v2",
            package = .pkgName,
            mustWork = TRUE
        ),
        "sampleMetadataFile" = system.file(
            "extdata", "cellranger_v2.csv",
            package = .pkgName,
            mustWork = TRUE
        )
    ),
    "v3" = list(
        "dir" = system.file(
            "extdata", "cellranger_v3",
            package = .pkgName,
            mustWork = TRUE
        ),
        "sampleMetadataFile" = system.file(
            "extdata", "cellranger_v3.csv",
            package = .pkgName,
            mustWork = TRUE
        )
    )
)

test_that("v2 : MTX : Fast mode", {
    x <- CellRanger(dir = lst[["v2"]][["dir"]])
    expect_s4_class(x, "CellRanger")
})

test_that("v2 : MTX : User-defined sample metadata", {
    object <- CellRanger(
        dir = lst[["v2"]][["dir"]],
        sampleMetadataFile = lst[["v2"]][["sampleMetadataFile"]]
    )
    expect_s4_class(object, "CellRanger")
    expect_identical(
        object = levels(colData(object)[["sampleName"]]),
        expected = "pbmc4k"
    )
})

test_that("v2 : MTX : AnnotationHub", {
    object <- CellRanger(
        dir = lst[["v2"]][["dir"]],
        organism = "Homo sapiens",
        ensemblRelease = 87L
    )
    expect_s4_class(object, "CellRanger")
})

test_that("v3 : MTX : Fast mode", {
    object <- CellRanger(lst[["v3"]][["dir"]])
    expect_s4_class(object, "CellRanger")
})

test_that("v3 MTX : User-defined sample metadata", {
    object <- CellRanger(
        dir = lst[["v3"]][["dir"]],
        sampleMetadataFile = lst[["v3"]][["sampleMetadataFile"]]
    )
    expect_s4_class(object, "CellRanger")
    expect_identical(
        object = levels(colData(object)[["sampleName"]]),
        expected = "5k_pbmc_protein_v3"
    )
})
