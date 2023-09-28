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
