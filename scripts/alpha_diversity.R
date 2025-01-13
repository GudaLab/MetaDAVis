alpha_diversity_boxplot <- function(OTU_input, group_index, tax_index, index_method = "1", plot_method = "1", pvalue = "No", select_alpha_color_palette = "Set3") {
  # Ensure tax_index has appropriate row and column names
  colnames(tax_index)[1] <- "Species"
  rownames(tax_index) <- tax_index[, 1]
  taxmat <- as.matrix(tax_index)
  
  # Create phyloseq object
  OTU <- otu_table(OTU_input, taxa_are_rows = TRUE)
  TAX <- tax_table(taxmat)
  sampledata <- sample_data(group_index)
  physeq1 <- phyloseq(OTU, TAX, sampledata)
  
  # Validate color palette
  if (!select_alpha_color_palette %in% rownames(brewer.pal.info)) {
    stop("Invalid color palette. Please choose a valid palette from RColorBrewer.")
  }
  
  # Create dynamic colors
  colors <- colorRampPalette(brewer.pal(9, select_alpha_color_palette))(length(unique(sample_data(physeq1)$Condition)))
  
  # Define richness measures based on index_method
  richness_measures <- switch(
    index_method,
    "1" = "Observed",
    "2" = "Chao1",
    "3" = "ACE",
    "4" = "Shannon",
    "5" = "Simpson",
    "6" = "InvSimpson",
    "7" = "Fisher",
    "8" = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"),
    stop("Invalid index_method. Please choose a value between 1 and 8.")
  )
  
  # Define plot type
  plot_type <- ifelse(plot_method == "1", "boxplot", "violin")
  
  # Create richness plot
  create_richness_plot <- function(measures, plot_type) {
    p <- plot_richness(
      physeq1, 
      x = "Condition", 
      measures = measures, 
      color = "Condition"
    ) +
      theme_bw() +
      scale_color_manual(values = colors)
    
    if (plot_type == "boxplot") {
      p <- p + geom_boxplot()
    } else if (plot_type == "violin") {
      p <- p + geom_violin()
    }
    return(p)
  }
  
  # Generate the richness plot
  richness_plot <- create_richness_plot(richness_measures, plot_type)
  
  # Add statistical comparisons if specified
  if (pvalue %in% c("Yes", "star")) {
    comps <- combn(unique(sample_data(physeq1)$Condition), 2, simplify = FALSE)
    stat_layer <- stat_compare_means(
      comparisons = comps, 
      label = "p.format", 
      tip.length = 0.01, 
      method = "wilcox.test", 
      symnum.args = if (pvalue == "star") {
        list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
      } else NULL
    )
    richness_plot <- richness_plot + stat_layer
  }
  
  # Return the final plot and additional data
  return(list(
    plot = richness_plot, 
    data1 = estimate_richness(physeq1, measures = richness_measures), 
    data2 = taxmat
  ))
}
