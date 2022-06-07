## 10X Chromium Cell Ranger v2 example output.
## 4k PBMCs from a healthy donor.
## https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k
## Updated 2022-06-07.

## nolint start
suppressPackageStartupMessages({
    library(devtools)
    library(usethis)
    library(AcidBase)
    library(Matrix)
    library(pipette)
})
## nolint end

load_all()
datasetName <- "pbmc_v2"
dataRawDir <- "data-raw"

## Restrict to 2 MB.
## Use `pryr::object_size()` instead of `utils::object.size()`.
limit <- structure(2e6, class = "object_size")

## Complete dataset =============================================================
## Create the example dataset directory structure.
dir <- initDir(file.path(dataRawDir, datasetName))
sampleDir <- initDir(file.path(dir, "pbmc"))
outsDir <- initDir(file.path(sampleDir, "outs"))
counterDir <- initDir(file.path(sampleDir, "SC_RNA_COUNTER_CS"))
file.create(file.path(counterDir, "empty"))
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
untar(tarfile = tarfile, exdir = outsDir)
stopifnot(identical(dir(outsDir), "filtered_gene_bc_matrices"))
## Copy the example outs files.
files <- c("metrics_summary.csv", "molecule_info.h5")
invisible(lapply(
    X = files,
    FUN = function(file) {
        file.copy(
            from = file.path(dir, paste0(prefix, "_", file)),
            to = file.path(outsDir, file),
            overwrite = TRUE
        )
    }
))
## Using Ensembl 84 GTF annotations.
## Note that ensembldb only supports back to 87.
gffFile <- file.path(dataRawDir, "Homo_sapiens.GRCh38.84.gtf.gz")
if (!file.exists(gffFile)) {
    download.file(
        url = paste(
            "ftp://ftp.ensembl.org",
            "pub",
            "release-84",
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
    name = datasetName,
    object = object,
    dir = dataRawDir
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
lapply(coerceToList(object), object_size)
object_size(object)
stopifnot(object_size(object) < limit)
stopifnot(validObject(object))
pbmc_v2 <- object
use_data(pbmc_v2, compress = "xz", overwrite = TRUE)

## Example Cell Ranger v2 output ===============================================
input_dir <- dir
input_sampleDir <- sampleDir
input_outsDir <- outsDir
input_matrix_dir <- file.path(
    input_outsDir,
    "filtered_gene_bc_matrices",
    "GRCh38"
)
output_dir <- file.path(
    "inst",
    "extdata",
    "cellranger_v2"
)
unlink(output_dir, recursive = TRUE)
output_sampleDir <- initDir(file.path(output_dir, "pbmc"))
output_outsDir <- initDir(file.path(output_sampleDir, "outs"))
output_counterDir <- file.path(output_sampleDir, "SC_RNA_COUNTER_CS")
initDir(output_counterDir)
file.create(file.path(output_counterDir, "empty"))
output_matrix_dir <- file.path(
    output_outsDir,
    "filtered_gene_bc_matrices",
    "GRCh38"
)
initDir(output_matrix_dir)
## Copy metrics summary.
file.copy(
    from = file.path(input_outsDir, "metrics_summary.csv"),
    to = file.path(output_outsDir, "metrics_summary.csv"),
    overwrite = TRUE
)
## Prepare the sparse matrix.
counts <- Matrix::readMM(file = file.path(input_matrix_dir, "matrix.mtx"))
counts <- counts[seq_len(100), seq_len(100)]
Matrix::writeMM(counts, file = file.path(output_matrix_dir, "matrix.mtx"))
## FIXME What's the n_max argument here?
genes <- import(
    file = file.path(input_matrix_dir, "genes.tsv"),
    colnames = FALSE,
    n_max = 100
)
export(
    x = genes,
    file = file.path(output_matrix_dir, "genes.tsv"),
    col_names = FALSE
)
barcodes <- read_lines(
    file = file.path(input_matrix_dir, "barcodes.tsv"),
    n_max = 100
)
write_lines(
    x = barcodes,
    file = file.path(output_matrix_dir, "barcodes.tsv")
)
