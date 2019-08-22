## 10X Chromium Cell Ranger v2 example output.
## 4k PBMCs from a healthy donor.
## https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k
## Updated 2019-08-22.

library(usethis)
library(pryr)
library(readr)
library(Matrix)
library(basejump)

dataset_name <- "pbmc_v2"
data_raw_dir <- "data-raw"

## Restrict to 2 MB.
## Use `pryr::object_size()` instead of `utils::object.size()`.
limit <- structure(2e6, class = "object_size")

## Complete dataset =============================================================
## Create the example dataset directory structure.
dir <- initDir(file.path(data_raw_dir, dataset_name))
unlink(dir, recursive = TRUE)
sample_dir <- initDir(file.path(dir, "pbmc"))
outs_dir <- initDir(file.path(sample_dir, "outs"))
counter_dir <- file.path(sample_dir, "SC_RNA_COUNTER_CS")
initDir(counter_dir)
file.create(file.path(counter_dir, "empty"))

## Download the example files.
prefix <- "pbmc4k"
files <- paste0(
    prefix, "_",
    c(
        ## "filtered_gene_bc_matrices_h5.h5",  # 403 (2019-08-22)
        "filtered_gene_bc_matrices.tar.gz",
        "metrics_summary.csv",
        "molecule_info.h5",
        "raw_gene_bc_matrices_h5.h5",
        "raw_gene_bc_matrices.tar.gz"
    )
)
url <- pasteURL(
    "cf.10xgenomics.com",
    "samples",
    "cell-exp",
    "2.1.0",
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
tarfile <- file.path(dir, paste0(prefix, "_filtered_gene_bc_matrices.tar.gz"))
untar(tarfile = tarfile, exdir = outs_dir)
stopifnot(identical(dir(outs_dir), "filtered_gene_bc_matrices"))

## Using Ensembl 84 GTF annotations.
gtf_file <- file.path(data_raw_dir, "Homo_sapiens.GRCh38.84.gtf.gz")
if (!file.exists(gtf_file)) {
    download.file(
        url = paste(
            "ftp://ftp.ensembl.org",
            "pub",
            "release-84",
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

pbmc_v2 <- object
usethis::use_data(pbmc_v2, compress = "xz", overwrite = TRUE)

## Example Cell Ranger v2 output ===============================================
input_dir <- file.path(
    outs_dir,
    "filtered_gene_bc_matrices",
    "GRCh38"
)
stopifnot(dir.exists(input_dir))
output_dir <- file.path(
    "inst",
    "extdata",
    "cellranger_v2"
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
    "filtered_gene_bc_matrices",
    "GRCh38"
)
initDir(output_dir)

## Prepare the sparse matrix.
counts <- Matrix::readMM(file = file.path(input_dir, "matrix.mtx"))
counts <- counts[seq_len(100), seq_len(100)]
Matrix::writeMM(counts, file = file.path(output_dir, "matrix.mtx"))

genes <- read_tsv(
    file = file.path(input_dir, "genes.tsv"),
    col_names = FALSE,
    n_max = 100
)
write_tsv(
    x = genes,
    path = file.path(output_dir, "genes.tsv"),
    col_names = FALSE
)

barcodes <- read_lines(
    file = file.path(input_dir, "barcodes.tsv"),
    n_max = 100
)
write_lines(
    x = barcodes,
    path = file.path(output_dir, "barcodes.tsv")
)
