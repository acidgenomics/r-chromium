# Chromium =====================================================================
test_that("Chromium", {
    uploadDir = system.file("extdata/cellranger", package = "Chromium")

    # Minimal mode, with no metadata or annotations.
    # This is fast but doesn't slot a lot of useful info.
    x <- Chromium(uploadDir = uploadDir)
    expect_s4_class(x, "Chromium")

    # Automatic organism annotations from AnnotationHub.
    x <- Chromium(
        uploadDir = uploadDir,
        organism = "Homo sapiens",
        ensemblRelease = 87L
    )
    expect_s4_class(x, "Chromium")
})
