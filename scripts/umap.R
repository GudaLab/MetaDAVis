library(phyloseq)
library(bluster)
library(patchwork)
library(scater)
library(mia)
library(RColorBrewer)
library(dplyr)

umap_plot_table <- function(OTU_input, group_index, tax_index, method, kvalue, select_umap_color_palette) {
  # Prepare taxonomic data
  colnames(tax_index)[1] <- "Species"
  rownames(tax_index) <- tax_index[, 1]
  taxmat <- as.matrix(tax_index)
  
  # Prepare phyloseq object
  OTU <- otu_table(OTU_input, taxa_are_rows = TRUE)
  TAX <- tax_table(taxmat)
  sampledata <- sample_data(group_index)
  physeq1 <- phyloseq(OTU, TAX, sampledata)
  
  # Generate dynamic colors with enough values for all unique conditions
  num_conditions <- length(unique(sample_data(physeq1)$Condition))
  if (num_conditions > 9) {
    colors <- colorRampPalette(brewer.pal(9, select_umap_color_palette))(num_conditions)
  } else {
    colors <- brewer.pal(num_conditions, select_umap_color_palette)
  }
  
  # Convert to TreeSummarizedExperiment
  tse <- makeTreeSummarizedExperimentFromPhyloseq(physeq1)
  
  # Ensure kvalue is numeric
  kvalue <- as.numeric(kvalue)
  
  # Validate the selected method
  valid_methods <- c("counts", "rclr", "hellinger", "pa", "rank", "relabundance")
  if (!method %in% valid_methods) {
    stop("Invalid method. Please select one of: ", paste(valid_methods, collapse = ", "))
  }
  
  # Apply the transformation and run UMAP
  if (method != "counts") {
    tse <- transformCounts(tse, method = method)
  }
  tse <- runUMAP(tse, name = "UMAP", exprs_values = method)
  
  # Cluster rows based on kvalue
  graph_clusters <- clusterRows(t(assays(tse)[[method]]), NNGraphParam(k = kvalue))
  
  # Create UMAP plots
  plots <- plotUMAP(tse, colour_by = "Condition") +
    labs(title = method) +
    scale_color_manual(values = colors)
  
  plots1 <- plotUMAP(tse, colour_by = I(graph_clusters)) +
    labs(title = paste0("k = ", kvalue))
  
  # Combine plots
  umapplot <- plots + plots1
  
  # Round numeric values in plot data
  plots$data <- plots$data %>% mutate(across(where(is.numeric), ~ round(., 2)))
  plots1$data <- plots1$data %>% mutate(across(where(is.numeric), ~ round(., 2)))
  
  # Return results
  return(list(plot = umapplot, data = plots$data, data1 = plots1$data))
}
