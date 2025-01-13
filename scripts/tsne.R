library(phyloseq)
library(bluster)
library(patchwork)
library(scater)
library(mia)
library(RColorBrewer)
library(dplyr)

tsne_plot_table <- function(OTU_input, group_index, tax_index, method, dimension, select_tsne_color_palette) {
  # Prepare taxonomic data
  colnames(tax_index)[1] <- "Species"
  rownames(tax_index) <- tax_index[, 1]
  taxmat <- as.matrix(tax_index)
  
  # Prepare phyloseq object
  OTU <- otu_table(OTU_input, taxa_are_rows = TRUE)
  TAX <- tax_table(taxmat)
  sampledata <- sample_data(group_index)
  physeq1 <- phyloseq(OTU, TAX, sampledata)
  
  # Generate dynamic colors
  unique_conditions <- unique(sample_data(physeq1)$Condition)
  colors <- colorRampPalette(brewer.pal(9, select_tsne_color_palette))(length(unique_conditions))
  
  # Convert to TreeSummarizedExperiment
  tse <- makeTreeSummarizedExperimentFromPhyloseq(physeq1)
  
  # Ensure dimension is numeric
  dimension <- as.numeric(dimension)
  
  # Transform counts and run TSNE
  if (method == "counts") {
    tse <- runTSNE(tse, name = "TSNE", exprs_values = "counts", ncomponents = dimension)
  } else if (method %in% c("rclr", "hellinger", "pa", "rank", "relabundance")) {
    tse <- transformCounts(tse, method = method)
    tse <- runTSNE(tse, name = "TSNE", exprs_values = method, ncomponents = dimension)
  } else {
    stop("Invalid method. Please select a valid method: 'counts', 'rclr', 'hellinger', 'pa', 'rank', or 'relabundance'.")
  }
  
  # Plot TSNE
  plottsne <- plotReducedDim(
    tse,
    "TSNE",
    colour_by = "Condition",
    ncomponents = seq_len(dimension)
  ) + scale_color_manual(values = colors)
  
  # Round numeric values in the plot data
  plottsne$data <- plottsne$data %>% mutate(across(where(is.numeric), ~ round(., 2)))
  
  # Return the plot and data
  return(list(plot = plottsne, data = plottsne$data))
}
