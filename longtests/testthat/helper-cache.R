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
        "pbmc4k_filtered_gene_bc_matrices.tar.gz",
        "pbmc4k_molecule_info.h5",
        "pbmc4k_raw_gene_bc_matrices_h5.h5",
        "pbmc4k_raw_gene_bc_matrices.tar.gz",
        ## pbmc_v3.
        "5k_pbmc_protein_v3_filtered_feature_bc_matrix.h5",
        "5k_pbmc_protein_v3_filtered_feature_bc_matrix.tar.gz",
        "5k_pbmc_protein_v3_molecule_info.h5",
        "5k_pbmc_protein_v3_raw_feature_bc_matrix.h5",
        "5k_pbmc_protein_v3_raw_feature_bc_matrix.tar.gz"
    )
)
cacheDir <- lst[["cacheDir"]]
rm(lst)
dirs <- c(
    "pbmc_v2_hdf5" = file.path(cacheDir, "pbmc_v2_hdf5"),
    "pbmc_v2_mtx" = file.path(cacheDir, "pbmc_v2_mtx"),
    "pbmc_v3_hdf5" = file.path(cacheDir, "pbmc_v3_hdf5"),
    "pbmc_v3_mtx" = file.path(cacheDir, "pbmc_v3_mtx")
)
if (!goalie::isADir(dirs[["pbmc_v2_hdf5"]])) {
    AcidBase::initDir(dirs[["pbmc_v2_hdf5"]])
    file.copy(
        from = file.path(
            cacheDir,
            "pbmc4k_raw_gene_bc_matrices_h5.h5"
        ),
        to = file.path(
            dirs[["pbmc_v2_hdf5"]],
            "pbmc4k_raw_gene_bc_matrices_h5.h5"
        )
    )
}
if (!goalie::isADir(dirs[["pbmc_v2_mtx"]])) {
    AcidBase::initDir(dirs[["pbmc_v2_mtx"]])
    file.copy(
        from = file.path(
            cacheDir,
            "pbmc4k_filtered_gene_bc_matrices.tar.gz"
        ),
        to = file.path(
            dirs[["pbmc_v2_mtx"]],
            "pbmc4k_filtered_gene_bc_matrices.tar.gz"
        )
    )
    file.copy(
        from = file.path(
            cacheDir,
            "pbmc4k_raw_gene_bc_matrices.tar.gz"
        ),
        to = file.path(
            dirs[["pbmc_v2_mtx"]],
            "pbmc4k_raw_gene_bc_matrices.tar.gz"
        )
    )
    utils::untar(
        tarfile = file.path(
            dirs[["pbmc_v2_mtx"]],
            "pbmc4k_filtered_gene_bc_matrices.tar.gz"
        ),
        exdir = dirs[["pbmc_v2_mtx"]]
    )
    utils::untar(
        tarfile = file.path(
            dirs[["pbmc_v2_mtx"]],
            "pbmc4k_raw_gene_bc_matrices.tar.gz"
        ),
        exdir = dirs[["pbmc_v2_mtx"]]
    )
}
if (!goalie::isADir(dirs[["pbmc_v3_hdf5"]])) {
    AcidBase::initDir(dirs[["pbmc_v3_hdf5"]])
    file.copy(
        from = file.path(
            cacheDir,
            "5k_pbmc_protein_v3_filtered_feature_bc_matrix.h5"
        ),
        to = file.path(
            dirs[["pbmc_v3_hdf5"]],
            "5k_pbmc_protein_v3_filtered_feature_bc_matrix.h5"
        )
    )
    file.copy(
        from = file.path(
            cacheDir,
            "5k_pbmc_protein_v3_raw_feature_bc_matrix.h5"
        ),
        to = file.path(
            dirs[["pbmc_v3_hdf5"]],
            "5k_pbmc_protein_v3_raw_feature_bc_matrix.h5"
        )
    )
}
if (!goalie::isADir(dirs[["pbmc_v3_mtx"]])) {
    AcidBase::initDir(dirs[["pbmc_v3_mtx"]])
    file.copy(
        from = file.path(
            cacheDir,
            "5k_pbmc_protein_v3_filtered_feature_bc_matrix.tar.gz"
        ),
        to = file.path(
            dirs[["pbmc_v3_mtx"]],
            "5k_pbmc_protein_v3_filtered_feature_bc_matrix.tar.gz"
        )
    )
    file.copy(
        from = file.path(
            cacheDir,
            "5k_pbmc_protein_v3_raw_feature_bc_matrix.tar.gz"
        ),
        to = file.path(
            dirs[["pbmc_v3_mtx"]],
            "5k_pbmc_protein_v3_raw_feature_bc_matrix.tar.gz"
        )
    )
    utils::untar(
        tarfile = file.path(
            dirs[["pbmc_v3_mtx"]],
            "5k_pbmc_protein_v3_filtered_feature_bc_matrix.tar.gz"
        ),
        exdir = dirs[["pbmc_v3_mtx"]]
    )
    utils::untar(
        tarfile = file.path(
            dirs[["pbmc_v3_mtx"]],
            "5k_pbmc_protein_v3_raw_feature_bc_matrix.tar.gz"
        ),
        exdir = dirs[["pbmc_v3_mtx"]]
    )
}
