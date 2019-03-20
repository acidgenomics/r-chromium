# Chromium

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

Toolkit for 10X Genomics single-cell RNA-seq data.

🚨🚨🚨 NOT YET COMPATIBLE WITH CELLRANGER 3 🚨🚨🚨

## Installation

This is an [R][] package.

### [Bioconductor][] method

We recommend installing the package with [BiocManager][].

```r
if (!require("BiocManager")) {
    install.packages("BiocManager")
}
BiocManager::install("remotes")
BiocManager::install("steinbaugh/Chromium")
```

[BiocManager]: https://cran.r-project.org/package=BiocManager
[Bioconductor]: https://bioconductor.org
[R]: https://www.r-project.org
