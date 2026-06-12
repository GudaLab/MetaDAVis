# MetaDAVis

## Introduction

MetaDAVis, interactive Metagenome Data Analysis and Visualization, is a browser-based R Shiny application for analyzing and visualizing 16S and whole-metagenome sequencing results from kingdom to species level. It is designed for researchers who want to run common metagenomics analyses without writing R code.

MetaDAVis includes:

1. Data summary and abundance distribution
2. Diversity analysis
3. Dimension reduction
4. Correlation analysis
5. Heatmap
6. Differential abundance for two-group and multi-group comparisons
7. Bulk download for all completed outputs

## Use MetaDAVis Online

MetaDAVis is deployed at:

https://www.gudalab-rtools.net/MetaDAVis

## Local Requirements

Recommended versions:

- R >= 4.4.2
- RStudio >= 2024.12.0
- Bioconductor >= 3.20
- Shiny >= 1.10.0

Recommended local build tools:

- Windows: install Rtools for your R version from https://cran.r-project.org/bin/windows/Rtools/
- macOS: install Xcode Command Line Tools with `xcode-select --install`
- Ubuntu/Debian Linux: install development libraries before installing R packages:

```bash
sudo apt-get update
sudo apt-get install -y build-essential gfortran libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libharfbuzz-dev libfribidi-dev
```

- Red Hat/CentOS/Fedora Linux: install development libraries before installing R packages:

```bash
sudo yum install -y gcc gcc-c++ gcc-gfortran make libcurl-devel openssl-devel libxml2-devel fontconfig-devel freetype-devel libpng-devel libtiff-devel libjpeg-turbo-devel harfbuzz-devel fribidi-devel
```

## Install MetaDAVis Locally

### Option 1: Run Directly From GitHub

Open R or RStudio and run:

```r
install.packages("shiny", repos = "https://cloud.r-project.org")
library(shiny)
shiny::runGitHub("MetaDAVis", "GudaLab")
```

### Option 2: Download and Run From a Local Folder

Download or clone the MetaDAVis repository:

```bash
git clone https://github.com/GudaLab/MetaDAVis.git
cd MetaDAVis
```

Then open R or RStudio in the MetaDAVis folder and run:

```r
library(shiny)
runApp(".", launch.browser = TRUE)
```

If you downloaded a ZIP file instead of using Git, unzip it, open R/RStudio in that folder, and run the same `runApp()` command.

## Install All Required R Packages

MetaDAVis checks for missing packages at startup, but installing all dependencies first is recommended for a smoother local launch.

Run this once in R or RStudio:

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_packages <- c(
  "shiny",
  "DT",
  "shinythemes",
  "shinyFiles",
  "shinyjs",
  "shinydashboard",
  "ggplot2",
  "ggpubr",
  "vegan",
  "ggfortify",
  "ggplotify",
  "reshape2",
  "tibble",
  "scales",
  "dunn.test",
  "tidyr",
  "dplyr",
  "devtools",
  "patchwork",
  "GGally",
  "plotly",
  "zip",
  "filelock",
  "shinycssloaders",
  "RColorBrewer",
  "circlize"
)

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  install.packages(missing_cran, dependencies = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "phyloseq",
  "microbiome",
  "ComplexHeatmap",
  "qvalue",
  "scater",
  "DESeq2",
  "limma",
  "edgeR",
  "metagenomeSeq",
  "bluster",
  "mia",
  "lefser"
)

missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

github_packages <- c(
  "microsud/microbiomeutilities",
  "biobakery/maaslin3"
)

for (pkg in github_packages) {
  package_name <- sub(".*/", "", pkg)
  if (!requireNamespace(package_name, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}
```

If a GitHub package fails to install through `BiocManager::install()`, install it with `devtools`:

```r
devtools::install_github("microsud/microbiomeutilities")
devtools::install_github("biobakery/maaslin3")
```

## Tool and Package Summary

| MetaDAVis tool | Main R packages used |
| --- | --- |
| Upload and preprocessing | `tidyr`, `dplyr`, `tibble` |
| Group and individual abundance distribution | `ggplot2`, `RColorBrewer` |
| Alpha diversity | `phyloseq`, `ggplot2`, `ggpubr`, `RColorBrewer` |
| Beta diversity | `phyloseq`, `vegan`, `reshape2`, `ggplot2`, `ggpubr` |
| PCA-2D | `ggfortify`, `ggplot2`, `RColorBrewer`, `dplyr` |
| PCA-3D | `plotly`, `RColorBrewer` |
| t-SNE and UMAP | `phyloseq`, `bluster`, `patchwork`, `scater`, `mia`, `RColorBrewer`, `dplyr` |
| Taxa and sample correlation | `ggpubr`, `GGally`, `ggplot2`, `RColorBrewer`, `dplyr` |
| Heatmap | `ComplexHeatmap`, `circlize`, `scales`, `ggplotify`, `RColorBrewer` |
| Wilcoxon and t-test | `ggplot2`, `tibble`, `qvalue`, `ComplexHeatmap`, `ggplotify`, `RColorBrewer`, `dplyr` |
| metagenomeSeq | `metagenomeSeq`, `ggplot2`, `tibble`, `ComplexHeatmap`, `ggplotify`, `dplyr` |
| DESeq2 | `DESeq2`, `ggplot2`, `tibble`, `qvalue`, `ComplexHeatmap`, `ggplotify`, `dplyr` |
| LEfSe | `phyloseq`, `lefser`, `mia`, `tibble`, `ggplot2`, `RColorBrewer`, `dplyr` |
| MaAsLin3 | `maaslin3`, `tibble`, `dplyr` |
| Limma-Voom and edgeR | `limma`, `edgeR`, `ggplot2`, `tibble`, `ComplexHeatmap`, `ggplotify`, `circlize`, `dplyr` |
| Kruskal-Wallis and ANOVA | `ggplot2`, `tibble`, `qvalue`, `dunn.test`, `ComplexHeatmap`, `ggplotify`, `RColorBrewer` |
| Bulk download | `zip`, `filelock`, base R |
| User interface | `shiny`, `DT`, `shinythemes`, `shinyFiles`, `shinyjs`, `shinydashboard`, `shinycssloaders` |

## Usage

Open the Manual tab inside the MetaDAVis application for the step-by-step tutorial, example-data instructions, analysis options, and Bulk Download instructions.

Example data are included in:

```text
www/example_data
```

## Tested Platforms

MetaDAVis has been tested on:

- Linux: Red Hat and Ubuntu
- Windows: 10 and 11

## Citation

Jagadesan S, Guda C (2025) MetaDAVis: An R shiny application for metagenomic data analysis and visualization. PLoS ONE 20(4): e0319949. https://doi.org/10.1371/journal.pone.0319949
