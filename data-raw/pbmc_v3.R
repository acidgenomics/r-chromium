## 10X Chromium Cell Ranger v3 example output.
## 5k Peripheral blood mononuclear cells (PBMCs) from a healthy donor with cell
## surface proteins (v3 chemistry).
## https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.1.0/5k_pbmc_protein_v3
## Updated 2019-08-22.

library(usethis)
library(pryr)
library(readr)
library(Matrix)
library(basejump)

dataset_name <- "pbmc_v3"
data_raw_dir <- "data-raw"

## Restrict to 2 MB.
## Use `pryr::object_size()` instead of `utils::object.size()`.
limit <- structure(2e6, class = "object_size")

## Complete dataset ============================================================
## Create the example dataset directory structure.
dir <- initDir(file.path(data_raw_dir, dataset_name))
unlink(dir, recursive = TRUE)
sample_dir <- initDir(file.path(dir, "pbmc"))
outs_dir <- initDir(file.path(sample_dir, "outs"))
## Touch an empty file in the counter directory.
counter_dir <- file.path(sample_dir, "SC_RNA_COUNTER_CS")
initDir(counter_dir)
file.create(file.path(counter_dir, "empty"))
## Download the example files.
prefix <- "5k_pbmc_protein_v3"
files <- paste0(
    prefix, "_",
    c(
        "filtered_feature_bc_matrix.h5",
        "filtered_feature_bc_matrix.tar.gz",
        "metrics_summary.csv",
        "molecule_info.h5",
        "raw_feature_bc_matrix.h5",
        "raw_feature_bc_matrix.tar.gz"
    )
)
url <- pasteURL(
    "cf.10xgenomics.com",
    "samples",
    "cell-exp",
    "3.1.0",
    prefix,
    protocol = "http"
)
invisible(lapply(
    X = files,
    FUN = function(file, dir, url) {
        destfile <- file.path(dir, file)
        if (!file.exists(destfile)) {
            download.file(
                url = file.path(url, file),
                destfile = destfile
            )
        }
    },
    dir = dir,
    url = url
))

## Extract the filtered MTX matrix files.
tarfile <- file.path(dir, paste0(prefix, "_filtered_feature_bc_matrix.tar.gz"))
untar(tarfile = tarfile, exdir = outs_dir)
stopifnot(identical(dir(outs_dir), "filtered_feature_bc_matrix"))

## Copy the example outs files.
files <- c(
    "filtered_feature_bc_matrix.h5",
    "metrics_summary.csv",
    "molecule_info.h5"
)
invisible(lapply(
    X = files,
    FUN = function(file) {
        file.copy(
            from = file.path(dir, paste0(prefix, "_", file)),
            to = file.path(outs_dir, file),
            overwrite = TRUE
        )

    }
))
