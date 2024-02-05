# iMGDAVis
## Introduction
interactive Metagenome Data Analysis and Visualization (iMGDAVis) application analyzes 16S and whole metagenome sequence results at various levels (kingdom to species). It is a browser-based and user-friendly R Shiny application for researchers to analyze and visualize without programming proficiency. It comprises six functional analyses.
1.	Data Summary and Distribution
2.	Diversity analysis
3.	Dimension reduction
4.	Correlation analysis
5.	Heatmap
6.	Differential abundance (Two and multiple groups)

## Installation
In this tutorial, we will go through the installation and usage of each functional module using the example dataset. The iMGDAVis is publicly available at (). The example dataset is provided on the GitHub page ().

How to start iMGDAVis locally

Download the iMGDAVis application locally from the GitHub page ().

Requirement:

•	R (≥ 4.3.1), available at (https://www.r-project.org/)

•	RStudio (≥ 2023.09.0) available at (https://posit.co/download/rstudio-desktop/) 

•	Bioconductor (3.18) and 

•	Shiny (≥ 1.7.5)

This Application was tested in Linux (Red Hat) and Windows 10

Start an R session using RStudio and run the following commands to install the shiny package:
```
install.packages("shiny")
```
To run iMGDAVis by the following commands in R:
```
library(shiny)
shiny::runGitHub("iMGDAVis ","GudaLab")
```
Or 
Alternatively, download the source code from GitHub and run the following command in the R session using RStudio:
```
library(shiny)
runApp('/path/to/the/iMGDAVis-master', launch.browser=TRUE)
```
## Usage
Tutorial for iMDGAVis https://github.com/GudaLab/iMGDAVis/blob/main/www/manual/iMDGAVis_manual.docx
