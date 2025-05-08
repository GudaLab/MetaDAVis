# MetaDAVis
## Introduction
interactive Metagenome Data Analysis and Visualization (MetaDAVis) application analyzes 16S and whole metagenome sequence results at various levels (kingdom to species). It is a browser-based and user-friendly R Shiny application for researchers to analyze and visualize without programming proficiency. It comprises six functional analyses.
1.	Data Summary and Distribution
2.	Diversity analysis
3.	Dimension reduction
4.	Correlation analysis
5.	Heatmap
6.	Differential abundance (Two and multiple groups)

## use MetaDAVis online
iMGDAViz is deployed at: https://www.gudalab-rtools.net/MetaDAVis

## Installation
In this tutorial, we will go through the installation and usage of each functional module using the example dataset. The MetaDAVis is publicly available at (https://github.com/GudaLab/MetaDAVis). The example dataset is provided on the GitHub page (https://github.com/GudaLab/MetaDAVis/tree/main/www/example_data).

How to start MetaDAVis locally

Download the MetaDAVis application locally from the GitHub page (https://github.com/GudaLab/MetaDAVis).

### Requirement

•	R (≥ 4.4.2), available at (https://www.r-project.org/)

•	RStudio (≥ 2024.12.0) available at (https://posit.co/download/rstudio-desktop/) 

•	Bioconductor (≥ 3.20) and 

•	Shiny (≥ 1.10.0)

Start an R session using RStudio and run the following commands to install the shiny package:
if Bioconductor version is less than 3.19. Bioconductor could be updated by :
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```
```
install.packages("shiny")
```
To run MetaDAVis by the following commands in R:
```
library(shiny)
shiny::runGitHub("MetaDAVis","GudaLab")
```
Or 
Alternatively, download the source code from GitHub and run the following command in the R session using RStudio:
```
library(shiny)
runApp('/path/to/the/MetaDAVis-master', launch.browser=TRUE)
```
## Usage
Tutorial for MetaDAVis https://github.com/GudaLab/MetaDAVis/blob/main/www/manual/MetaDAVis_manual.pdf

## Tested
This Application was tested in Linux (Red Hat and Ubuntu) and Windows (10 and 11)

## Citation
Jagadesan S, Guda C (2025) MetaDAVis: An R shiny application for metagenomic data analysis and visualization. PLoS ONE 20(4): e0319949. (https://doi.org/10.1371/journal.pone.0319949)

## Contributors
- [@s.jagadesan@unmc.edu](https://github.com/s-jagadesan) – Maintainer
- [@janedoe](https://github.com/janedoe) – Feature X
