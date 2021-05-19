# Chromium

Toolkit for 10X Genomics Chromium single cell data.

## Installation

### [R][] method

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
install.packages(
    pkgs = "Chromium",
    repos = c(
        "https://r.acidgenomics.com",
        BiocManager::repositories()
    )
)
```

[r]: https://www.r-project.org/
