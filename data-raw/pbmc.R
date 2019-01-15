# FIXME Work on loading a v3 dataset and using as example.



# 10X Chromium Cell Ranger example output.
# 4k PBMCs from a healthy donor.
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k
# 2018-12-03

library(pryr)
library(tidyverse)
library(Matrix)

# Restrict to 2 MB.
# Use `pryr::object_size()` instead of `utils::object.size()`.
limit <- structure(2e6, class = "object_size")

# Complete dataset =============================================================
dir <- initDir("data-raw/cellranger")
# Example dataset contains a single sample ("pbmc4k").
outs_dir <- initDir(file.path(dir, "pbmc", "outs"))

# Directory structure:
# - pbmc
# - outs
# - filtered_gene_bc_matrices
# - GRCh38
tar_file <- "data-raw/pbmc.tar.gz"
if (!file.exists(tar_file)) {
    download.file(
        url = paste(
            "http://cf.10xgenomics.com",
            "samples",
            "cell-exp",
            "2.1.0",
            "pbmc4k",
            "pbmc4k_filtered_gene_bc_matrices.tar.gz",
            sep = "/"
        ),
        destfile = tar_file
    )
}
untar(tarfile = tar_file, exdir = outs_dir)
stopifnot(identical(dir(outs_dir), "filtered_gene_bc_matrices"))

# Using Ensembl 84 GTF annotations.
gtf_file <- file.path(dir, "genes.gtf")
if (!file.exists(gtf_file)) {
    download.file(
        url = paste(
            "ftp://ftp.ensembl.org",
            "pub",
            "release-84",
            "gtf",
            "homo_sapiens",
            "Homo_sapiens.GRCh38.84.gtf.gz",
            sep = "/"
        ),
        destfile = gtf_file
    )
}

# Note that this blows up in memory too much to run on RStudio AMI.
pbmc <- Chromium(dir = dir, organism = "Homo sapiens", gffFile = gtf_file)
# We're using a subset of this object for our working example (see below).
object_size(pbmc)
assignAndSaveData(name = "pbmc", object = pbmc, dir = "data-raw")

# Example object ===============================================================
counts <- counts(pbmc)

# Subset the matrix to include only the top genes and cells.
top_genes <- Matrix::rowSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
genes <- sort(names(top_genes))

top_cells <- Matrix::colSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 100L)
cells <- sort(names(top_cells))

# Subset the original pbmc dataset to contain only top genes and cells.
pbmc <- pbmc[genes, cells]

# Include only minimal metadata columns in rowRanges.
rowRanges(pbmc) <- rowRanges(pbmc) %>%
    .[, c("broadClass", "geneBiotype", "geneID", "geneName")] %>%
    relevelRowRanges()

# Report the size of each slot in bytes.
vapply(
    X = coerceS4ToList(pbmc),
    FUN = object_size,
    FUN.VALUE = numeric(1L)
)
object_size(pbmc)
stopifnot(object_size(pbmc) < limit)
stopifnot(validObject(pbmc))

usethis::use_data(pbmc, compress = "xz", overwrite = TRUE)



# Example Cell Ranger output ===================================================
input_dir <- file.path(
    dir,
    "pbmc",
    "outs",
    "filtered_gene_bc_matrices",
    "GRCh38"
)
stopifnot(dir.exists(input_dir))
output_dir <- file.path(
    "inst",
    "extdata",
    "cellranger",
    "pbmc",
    "outs",
    "filtered_gene_bc_matrices",
    "GRCh38"
)
unlink("inst/extdata/cellranger", recursive = TRUE)
dir.create(output_dir, recursive = TRUE)

# Prepare the sparse matrix.
counts <- Matrix::readMM(file = file.path(input_dir, "matrix.mtx"))
counts <- counts[seq_len(100), seq_len(100)]
Matrix::writeMM(counts, file = file.path(output_dir, "matrix.mtx"))

genes <- readr::read_tsv(
    file = file.path(input_dir, "genes.tsv"),
    col_names = c("geneID", "geneName"),
    n_max = 100
)
readr::write_tsv(
    x = genes,
    path = file.path(output_dir, "genes.tsv"),
    col_names = FALSE
)

barcodes <- readr::read_lines(
    file = file.path(input_dir, "barcodes.tsv"),
    n_max = 100
)
readr::write_lines(
    x = barcodes,
    path = file.path(output_dir, "barcodes.tsv")
)
