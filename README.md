# Chromium

[![Repo status: active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Travis CI build status](https://travis-ci.com/acidgenomics/Chromium.svg?branch=master)](https://travis-ci.com/acidgenomics/Chromium)
[![AppVeyor CI build status](https://ci.appveyor.com/api/projects/status/kq9ecwl1nktap64f/branch/master?svg=true)](https://ci.appveyor.com/project/mjsteinbaugh/chromium/branch/master)

Toolkit for 10X Genomics Chromium single cell data.

## Installation

### [R][] method

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}
Sys.setenv(R_REMOTES_UPGRADE = "always")
# Set `GITHUB_PAT` in `~/.Renviron` if you get a rate limit error.
remotes::install_github("acidgenomics/Chromium")
remotes::update_packages()
```

[Bioconductor]: https://bioconductor.org/
[BiocManager]: https://cran.r-project.org/package=BiocManager
[R]: https://www.r-project.org/
