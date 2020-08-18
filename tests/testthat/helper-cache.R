if (!isTRUE(goalie::hasInternet())) {
    warning("No Internet connection detected.")
    return(invisible())
}

## pbmc_v2
dir.create(
    path = file.path("cache", "pbmc_v2"),
    showWarnings = FALSE,
    recursive = TRUE
)
file <- "pbmc4k_raw_gene_bc_matrices_h5.h5"
url <- paste(ChromiumTestsURL, file, sep = "/")
destfile <- file.path("cache", "pbmc_v2", file)
if (!file.exists(destfile)) {
    utils::download.file(url = url, destfile = destfile)
}

## pbmc_v3
dir.create(
    path = file.path("cache", "pbmc_v3"),
    showWarnings = FALSE,
    recursive = TRUE
)
file <- "5k_pbmc_protein_v3_filtered_feature_bc_matrix.h5"
url <- paste(ChromiumTestsURL, file, sep = "/")
destfile <- file.path("cache", "pbmc_v3", file)
if (!file.exists(destfile)) {
    utils::download.file(url = url, destfile = destfile)
}
