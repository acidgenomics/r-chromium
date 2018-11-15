# CellRanger ===================================================================
test_that("CellRanger", {
    uploadDir = system.file("extdata/cellranger", package = "bcbioSingleCell")

    # Minimal mode, with no metadata or annotations.
    # This is fast but doesn't slot a lot of useful info.
    x <- CellRanger(uploadDir = uploadDir)
    expect_s4_class(x, "CellRanger")

    # Automatic organism annotations from AnnotationHub.
    x <- CellRanger(
        uploadDir = uploadDir,
        organism = "Homo sapiens",
        ensemblRelease = 87L
    )
    expect_s4_class(x, "CellRanger")
})
