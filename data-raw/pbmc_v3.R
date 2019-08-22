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

## Using Ensembl 93 GTF annotations.
gtf_file <- file.path(data_raw_dir, "Homo_sapiens.GRCh38.93.gtf.gz")
if (!file.exists(gtf_file)) {
    download.file(
        url = paste(
            "ftp://ftp.ensembl.org",
            "pub",
            "release-93",
            "gtf",
            "homo_sapiens",
            basename(gtf_file),
            sep = "/"
        ),
        destfile = gtf_file
    )
}

object <- CellRanger(
    dir = dir,
    organism = "Homo sapiens",
    gffFile = gtf_file
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
top_genes <- Matrix::rowSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
genes <- sort(names(top_genes))

top_cells <- Matrix::colSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 100L)
cells <- sort(names(top_cells))

## Subset the original object dataset to contain only top genes and cells.
object <- object[genes, cells]

## Report the size of each slot in bytes.
lapply(coerceS4ToList(object), object_size)
object_size(object)
stopifnot(object_size(object) < limit)
stopifnot(validObject(object))

pbmc_v3 <- object
usethis::use_data(pbmc_v3, compress = "xz", overwrite = TRUE)

## Example Cell Ranger v3 output ===============================================
input_dir <- file.path(
    outs_dir,
    "filtered_feature_bc_matrix"
)
stopifnot(dir.exists(input_dir))
output_dir <- file.path(
    "inst",
    "extdata",
    "cellranger_v3"
)
unlink(output_dir, recursive = TRUE)
sample_dir <- initDir(file.path(output_dir, "pbmc"))
outs_dir <- initDir(file.path(sample_dir, "outs"))
## Touch an empty file in the counter directory.
counter_dir <- file.path(sample_dir, "SC_RNA_COUNTER_CS")
initDir(counter_dir)
file.create(file.path(counter_dir, "empty"))
output_dir <- file.path(
    outs_dir,
    "filtered_feature_bc_matrix"
)
initDir(output_dir)

## Prepare the sparse matrix.
counts <- Matrix::readMM(file = file.path(input_dir, "matrix.mtx.gz"))
counts <- counts[seq_len(100), seq_len(100)]
Matrix::writeMM(counts, file = file.path(output_dir, "matrix.mtx.gz"))

features <- read_tsv(
    file = file.path(input_dir, "features.tsv.gz"),
    col_names = FALSE,
    n_max = 100
)
write_tsv(
    x = genes,
    path = file.path(output_dir, "features.tsv.gz"),
    col_names = FALSE
)

barcodes <- read_lines(
    file = file.path(input_dir, "barcodes.tsv.gz"),
    n_max = 100
)
write_lines(
    x = barcodes,
    path = file.path(output_dir, "barcodes.tsv.gz")
)
