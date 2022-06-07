## 10X Chromium Cell Ranger v2 example output.
## 
## 4k PBMCs from a healthy donor.
##
## Updated 2022-06-07.
##
## https://support.10xgenomics.com/single-cell-gene-expression/datasets/
## 2.1.0/pbmc4k

## FIXME This file is now failing to import:
## /Users/mike/git/monorepo/r-packages/r-chromium/data-raw/data-raw/pbmc_v2/pbmc/outs/filtered_gene_bc_matrices/GRCh38/matrix.mtx

## nolint start
suppressPackageStartupMessages({
    library(devtools)
    library(usethis)
    library(AcidBase)
    library(AcidExperiment)
    library(pipette)
})
## nolint end

load_all()
datasetName <- "pbmc_v2"
## Restrict to 2 MB.
limit <- structure(2e6L, class = "object_size")

## Complete dataset ============================================================
## Create the example dataset directory structure.
dir <- initDir(datasetName)
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
gffFile <- "Homo_sapiens.GRCh38.84.gtf.gz"
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
assignAndSaveData(
    name = datasetName,
    object = object,
    dir = getwd()
)

## Example object ==============================================================
counts <- counts(object)
topGenes <-
    counts |>
    Matrix::rowSums() |>
    sort(decreasing = TRUE) |>
    head(n = 500L)
genes <- sort(names(topGenes))
topCells <-
    counts |>
    Matrix::colSums() |>
    sort(decreasing = TRUE) |>
    head(n = 100L)
cells <- sort(names(topCells))
object <- object[genes, cells]
stopifnot(
    object.size(object) < limit,
    validObject(object)
)
pbmc_v2 <- object # nolint
use_data(pbmc_v2, compress = "xz", overwrite = TRUE)

## Example Cell Ranger v2 output ===============================================
inputDir <- dir
inputSampleDir <- sampleDir
inputOutsDir <- outsDir
inputMatrixDir <- file.path(
    inputOutsDir,
    "filtered_gene_bc_matrices",
    "GRCh38"
)
outputDir <- file.path(
    "..",
    "inst",
    "extdata",
    "cellranger_v2"
)
if (dir.exists(outputDir)) {
    unlink2(outputDir)
}
outputSampleDir <- initDir(file.path(outputDir, "pbmc"))
outputOutsDir <- initDir(file.path(outputSampleDir, "outs"))
outputCounterDir <- initDir(file.path(outputSampleDir, "SC_RNA_COUNTER_CS"))
file.create(file.path(outputCounterDir, "empty"))
outputMatrixDir <- initDir(file.path(
    outputOutsDir,
    "filtered_gene_bc_matrices",
    "GRCh38"
))
file.copy(
    from = file.path(inputOutsDir, "metrics_summary.csv"),
    to = file.path(outputOutsDir, "metrics_summary.csv"),
    overwrite = TRUE
)
counts <- import(file = file.path(inputMatrixDir, "matrix.mtx"))
counts <- counts[seq_len(100L), seq_len(100L)]
export(
    object = counts,
    con = file.path(outputMatrixDir, "matrix.mtx")
)
genes <- import(
    file = file.path(inputMatrixDir, "genes.tsv"),
    colnames = FALSE,
    nMax = 100L
)
stopifnot(identical(nrow(genes), 100L))
export(
    object = genes,
    con = file.path(outputMatrixDir, "genes.tsv"),
    colnames = FALSE
)
barcodes <- import(
    file = file.path(inputMatrixDir, "barcodes.tsv"),
    nMax = 100L
)
stopifnot(identical(nrow(barcodes), 100L))
export(
    object = barcodes,
    file = file.path(outputMatrixDir, "barcodes.tsv"),
    colnames = FALSE
)
