# Chromium 0.2.0 (2022-06-07)

## Major changes

- Now requiring R 4.2 / Bioconductor 3.15.
- Split out basejump imports in NAMESPACE into component packages.
- Updated lintr and testthat checks.

## Minor changes

- Resaved `pbmc_v2` and `pbmc_v3` example objects.

# Chromium 0.1.9 (2021-03-12)

## Minor changes

- Updated basejump dependencies and removed unnecessary stringr import.

# Chromium 0.1.8 (2021-03-04)

## Minor changes

- Reworked and simplified internal NAMESPACE using basejump v0.14.

# Chromium 0.1.7 (2020-10-12)

## Minor changes

- Updated minimum Acid Genomics dependency versions.

# Chromium 0.1.6 (2020-07-23)

## Minor changes

- Maintenance release, bumping the minimum R version to 4.0.

# Chromium 0.1.5 (2020-03-11)

## Minor changes

- NAMESPACE fix removing `isSpike`, which is now defunct in Bioconductor 3.11
  SingleCellExperiment update.

# Chromium 0.1.4 (2020-02-24)

## Minor changes

- `CellRanger`: Removed now defunct `spikeNames` argument. Refer to
  `makeSummarizedExperiment`, `makeSingleCellExperiment` documentation and
  release notes in basejump package for details.

# Chromium 0.1.3 (2020-01-27)

## Minor changes

- Updated basejump dependencies.
- Now using cli package for messages.

# Chromium 0.1.2 (2019-10-30)

## Minor changes

- Compatibility fixes for Bioconductor 3.10.

# Chromium 0.1.1 (2019-08-28)

## Minor changes

- NAMESPACE and basejump dependency updates.
- Now using `leftJoin` instead of `left_join` for `DataFrame`.

# Chromium 0.1.0 (2019-08-23)

Initial release.
