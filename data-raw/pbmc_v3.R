## 10X Chromium Cell Ranger v3 example output.
##
## 5k Peripheral blood mononuclear cells (PBMCs) from a healthy donor with cell
## surface proteins (v3 chemistry).
##
## Updated 2022-06-07.
##
## https://support.10xgenomics.com/single-cell-gene-expression/datasets/
## 3.1.0/5k_pbmc_protein_v3

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
datasetName <- "pbmc_v3"
limit <- structure(2e6L, class = "object_size") # nolint

## Complete dataset ============================================================
## Create the example dataset directory structure.
dir <- initDir(datasetName)
if (dir.exists(dir)) {
    unlink2(dir)
}
sampleDir <- initDir(file.path(dir, "pbmc"))
outsDir <- initDir(file.path(sampleDir, "outs"))
counterDir <- initDir(file.path(sampleDir, "SC_RNA_COUNTER_CS"))
file.create(file.path(counterDir, "empty"))
urlStem <- pasteUrl(
    "cf.10xgenomics.com",
    "samples",
    "cell-exp",
    "3.1.0",
    "5k_pbmc_protein_v3",
    protocol = "http"
)
## nolint start
urls <- paste(
    urlStem,
    paste0(
        "5k_pbmc_protein_v3_",
        c(
            "filtered_feature_bc_matrix.h5",
            "filtered_feature_bc_matrix.tar.gz",
            "metrics_summary.csv",
            "molecule_info.h5",
            "raw_feature_bc_matrix.h5",
            "raw_feature_bc_matrix.tar.gz"
        )
    ),
    sep = "/"
)
## nolint end
files <- vapply(
    X = urls,
    pkg = .pkgName,
    FUN = cacheUrl,
    FUN.VALUE = character(1L),
    USE.NAMES = FALSE
)
names(files) <- gsub(
    pattern = "^5k_pbmc_protein_v3_",
    replacement = "",
    x = basename(urls)
)
untar(
    tarfile = files[["filtered_feature_bc_matrix.tar.gz"]],
    exdir = outsDir
)
file.copy(
    from = files[["filtered_feature_bc_matrix.h5"]],
    to = file.path(outsDir, "filtered_feature_bc_matrix.h5")
)
file.copy(
    from = files[["metrics_summary.csv"]],
    to = file.path(outsDir, "metrics_summary.csv")
)
file.copy(
    from = files[["molecule_info.h5"]],
    to = file.path(outsDir, "molecule_info.h5")
)
gffFile <- cacheUrl(
    url = pasteUrl(
        "ftp.ensembl.org",
        "pub",
        "release-93",
        "gtf",
        "homo_sapiens",
        "Homo_sapiens.GRCh38.93.gtf.gz",
        protocol = "ftp"
    ),
    pkg = .pkgName
)
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
## Subset the matrix to include only the top genes and cells.
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
## Subset the original object dataset to contain only top genes and cells.
object <- object[genes, cells]
## Report the size of each slot in bytes.
stopifnot(
    object.size(object) < limit,
    validObject(object)
)
pbmc_v3 <- object # nolint
use_data(pbmc_v3, compress = "xz", overwrite = TRUE)

## Example Cell Ranger v3 output ===============================================
inputDir <- dir
inputSampleDir <- sampleDir
inputOutsDir <- outsDir
inputMatrixDir <- file.path(
    inputOutsDir,
    "filtered_feature_bc_matrix"
)
outputDir <- file.path(
    "..",
    "inst",
    "extdata",
    "cellranger_v3"
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
    "filtered_feature_bc_matrix"
))
file.copy(
    from = file.path(inputOutsDir, "metrics_summary.csv"),
    to = file.path(outputOutsDir, "metrics_summary.csv"),
    overwrite = TRUE
)
## Prepare the sparse matrix.
counts <- import(file.path(inputMatrixDir, "matrix.mtx.gz"))
counts <- counts[seq_len(100L), seq_len(100L)]
export(
    object = counts,
    con = file.path(outputMatrixDir, "matrix.mtx.gz")
)
features <- import(
    con = file.path(inputMatrixDir, "features.tsv.gz"),
    colnames = FALSE,
    nMax = 100L
)
expect_identical(nrow(features), 100L)
export(
    object = features,
    con = file.path(outputMatrixDir, "features.tsv.gz"),
    colnames = FALSE
)
barcodes <- import(
    con = file.path(inputMatrixDir, "barcodes.tsv.gz"),
    nMax = 100L
)
expect_identical(nrow(barcodes), 100L)
export(
    object = barcodes,
    con = file.path(outputMatrixDir, "barcodes.tsv.gz"),
    colnames = FALSE
)
