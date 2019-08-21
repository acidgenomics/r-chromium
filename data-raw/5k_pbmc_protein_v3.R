## 5k Peripheral blood mononuclear cells (PBMCs) from a healthy donor with cell
## surface proteins (v3 chemistry)
## 
## Single Cell Gene Expression Dataset by Cell Ranger 3.1.0
##
## Peripheral blood mononuclear cells (PBMCs) from a healthy donor cells were
## stained with a panel of 31 TotalSeq™-B antibodies (BioLegend).
##
## PBMCs are primary cells with relatively small amounts of RNA (~1pg RNA/cell).
##
## Libraries were prepared following the Chromium Single Cell 3ʹ Reagent Kits v3
## with Feature Barcoding technology for Cell Surface Protein User Guide
## (CG000185 RevB).
##
## - 5,247 cells detected
## - Sequenced on Illumina NovaSeq with approximately 28,918 reads per cell
## - 28bp read1 (16bp Chromium barcode and 12bp UMI)
## - 91bp read2 (transcript)
## - 8bp I7 sample barcode
## - run with `--expect-cells=5000`
##
## Published on July 24, 2019
##
## This dataset is licensed under the Creative Commons Attribution license.
##
## https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.1.0/5k_pbmc_protein_v3

dataset <- "5k_pbmc_protein_v3"
initDir(dataset)
prefix <- pasteURL(
    "cf.10xgenomics.com",
    "samples",
    "cell-exp",
    "3.1.0",
    dataset,
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
