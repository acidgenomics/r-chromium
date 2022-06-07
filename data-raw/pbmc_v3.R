## 10X Chromium Cell Ranger v3 example output.
## 5k Peripheral blood mononuclear cells (PBMCs) from a healthy donor with cell
## surface proteins (v3 chemistry).
## https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.1.0/5k_pbmc_protein_v3
## Updated 2021-03-03.

## nolint start
suppressPackageStartupMessages({
    library(devtools)
    library(usethis)
    library(pryr)
    library(readr)
    library(Matrix)
    library(basejump)
})
## nolint end

load_all()
dataset_name <- "pbmc_v3"
data_raw_dir <- "data-raw"
## Restrict to 2 MB.
limit <- structure(2e6L, class = "object_size")

## Complete dataset ============================================================
## Create the example dataset directory structure.
dir <- initDir(file.path(data_raw_dir, dataset_name))
unlink(dir, recursive = TRUE)
sample_dir <- initDir(file.path(dir, "pbmc"))
outs_dir <- initDir(file.path(sample_dir, "outs"))
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
## Using Ensembl 93 GTF annotations.
## Alternatively, can use ensembldb here.
gffFile <- file.path(data_raw_dir, "Homo_sapiens.GRCh38.93.gtf.gz")
if (!file.exists(gffFile)) {
    download.file(
        url = paste(
            "ftp://ftp.ensembl.org",
            "pub",
            "release-93",
            "gtf",
            "homo_sapiens",
            basename(gffFile),
            sep = "/"
        ),
        destfile = gffFile
    )
}
object <- CellRanger(
    dir = dir,
    organism = "Homo sapiens",
    gffFile = gffFile
)
## We're using a subset of this object for our working example (see below).
assignAndSaveData(
    name = dataset_name,
    object = object,
    dir = data_raw_dir
)

## Example object ==============================================================
counts <- counts(object)
## Subset the matrix to include only the top genes and cells.
top_genes <-
    counts |>
    rowSums() |>
    sort(decreasing = TRUE) |>
    head(n = 500L)
genes <- sort(names(top_genes))
top_cells <-
    counts |>
    colSums() |>
    sort(decreasing = TRUE) |>
    head(n = 100L)
cells <- sort(names(top_cells))
## Subset the original object dataset to contain only top genes and cells.
object <- object[genes, cells]
## Report the size of each slot in bytes.
stopifnot(
    object.size(object) < limit,
    validObject(object)
)
pbmc_v3 <- object # nolint
usethis::use_data(pbmc_v3, compress = "xz", overwrite = TRUE)

## Example Cell Ranger v3 output ===============================================
input_dir <- dir
input_sample_dir <- sample_dir
input_outs_dir <- outs_dir
input_matrix_dir <- file.path(
    input_outs_dir,
    "filtered_feature_bc_matrix"
)
output_dir <- file.path(
    "inst",
    "extdata",
    "cellranger_v3"
)
unlink(output_dir, recursive = TRUE)
output_sample_dir <- initDir(file.path(output_dir, "pbmc"))
output_outs_dir <- initDir(file.path(output_sample_dir, "outs"))
output_counter_dir <- file.path(output_sample_dir, "SC_RNA_COUNTER_CS")
initDir(output_counter_dir)
file.create(file.path(output_counter_dir, "empty"))
output_matrix_dir <- file.path(
    output_outs_dir,
    "filtered_feature_bc_matrix"
)
initDir(output_matrix_dir)
## Copy metrics summary.
file.copy(
    from = file.path(input_outs_dir, "metrics_summary.csv"),
    to = file.path(output_outs_dir, "metrics_summary.csv"),
    overwrite = TRUE
)
## Prepare the sparse matrix.
counts <- Matrix::readMM(file = file.path(input_matrix_dir, "matrix.mtx.gz"))
counts <- counts[seq_len(100), seq_len(100)]
Matrix::writeMM(counts, file = file.path(output_matrix_dir, "matrix.mtx.gz"))
features <- read_tsv(
    file = file.path(input_matrix_dir, "features.tsv.gz"),
    col_names = FALSE,
    n_max = 100
)
write_tsv(
    x = features,
    file = file.path(output_matrix_dir, "features.tsv.gz"),
    col_names = FALSE
)
barcodes <- read_lines(
    file = file.path(input_matrix_dir, "barcodes.tsv.gz"),
    n_max = 100
)
write_lines(
    x = barcodes,
    file = file.path(output_matrix_dir, "barcodes.tsv.gz")
)
