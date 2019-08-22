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

dataset_name <- "pbmc5k_v3"
data_raw_dir <- "data-raw"

## Restrict to 2 MB.
## Use `pryr::object_size()` instead of `utils::object.size()`.
limit <- structure(2e6, class = "object_size")





## Complete dataset ============================================================
dir <- initDir(file.path(data_raw_dir, dataset_name))
unlink(dir, recursive = TRUE)
sample_dir <- initDir(file.path(dir, "pbmc"))
outs_dir <- initDir(file.path(sample_dir, "outs"))
## Touch an empty file in the counter directory.
counter_dir <- file.path(sample_dir, "SC_RNA_COUNTER_CS")
initDir(counter_dir)
file.create(file.path(counter_dir, "empty"))









prefix <- pasteURL(
    "cf.10xgenomics.com",
    "samples",
    "cell-exp",
    "3.1.0",
    "5k_pbmc_protein_v3",
    protocol = "http"
)
files <- c(
    ## Per-molecule read information
    "molecule_info.h5",
    ## Feature / cell matrix HDF5 (filtered)
    "filtered_feature_bc_matrix.h5",
    ## Feature / cell matrix (filtered)
    "filtered_feature_bc_matrix.tar.gz"
)
invisible(lapply(
    X = files,
    FUN = function(file, dataset) {
        destfile <- file.path(dataset, file)
        if (!file.exists(destfile)) {
            download.file(
                url = file.path(prefix, paste0(dataset, "_", file)),
                destfile = destfile
            )
        }
    },
    dataset = dataset
))
