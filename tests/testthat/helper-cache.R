## 10X Genomics example datasets:
## - 4k PBMCs from a Healthy Donor
##   Single Cell Gene Expression Dataset by Cell Ranger 2.1.0
##   https://support.10xgenomics.com/single-cell-gene-expression/
##     datasets/2.1.0/pbmc4k
## - 5k Peripheral Blood Mononuclear Cells (PBMCs) from a Healthy Donor with a
##   Panel of TotalSeq™-B Antibodies (v3 chemistry) Single Cell Gene Expression
##   Dataset by Cell Ranger 3.1.0
##   https://www.10xgenomics.com/resources/datasets/5-k-peripheral-blood-
##     mononuclear-cells-pbm-cs-from-a-healthy-donor-with-cell-surface-proteins-
##     v-3-chemistry-3-1-standard-3-1-0
lst <- AcidDevTools::cacheTestFiles(
    pkg = .pkgName,
    files = c(
        ## pbmc_v2.
        "pbmc4k_raw_gene_bc_matrices_h5.h5",
        ## pbmc_v3.
        "5k_pbmc_protein_v3_filtered_feature_bc_matrix.h5"
    )
)
cacheDir <- lst[["cacheDir"]]
rm(lst)
