# Requires: phyloseq, ggplot2, ggpubr, RColorBrewer

alpha_diversity_boxplot <- function(
    OTU_input,
    group_index,
    tax_index,
    index_method = "1",
    plot_method  = "1",
    pvalue = "No",
    select_alpha_color_palette = "Set3",
    show_head = TRUE,
    head_n = 5,
    show_colors = TRUE
){
  log_head <- function(x, label){
    if (show_head) {
      msg <- paste0(
        "\n--- ", label, " (head ", head_n, ") ---\n",
        paste(utils::capture.output(utils::head(x, head_n)), collapse = "\n"),
        "\n"
      )
      cat(msg, file = stderr())
      flush.console()
    }
  }
  log_head_vec <- function(x, label){
    if (show_head) {
      x <- as.vector(x)
      msg <- paste0(
        "\n--- ", label, " (first ", head_n, ") ---\n",
        paste(utils::capture.output(utils::head(x, head_n)), collapse = "\n"),
        "\n"
      )
      cat(msg, file = stderr())
      flush.console()
    }
  }
  
  # Ensure tax_index has appropriate row and column names
  colnames(tax_index)[1] <- "Species"
  rownames(tax_index) <- tax_index[, 1, drop = TRUE]
  log_head(tax_index, "tax_index (with Species col & rownames)")
  
  taxmat <- as.matrix(tax_index)
  
  # Create phyloseq object
  OTU <- otu_table(OTU_input, taxa_are_rows = TRUE)
  sampledata <- sample_data(group_index)
  TAX <- tax_table(taxmat)
  physeq1 <- phyloseq(OTU, TAX, sampledata)
  
  log_head(OTU_input,   "OTU_input (raw)")
  log_head(group_index, "group_index (sample metadata)")
  log_head(tax_index,   "tax_index (taxonomy matrix)")
  
  # Validate color palette
  if (!select_alpha_color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    stop("Invalid color palette. Please choose a valid palette from RColorBrewer.")
  }
  
  # Create dynamic colors
  conds <- unique(sample_data(physeq1)$Condition)
  n_cond <- length(conds)
  max_n <- RColorBrewer::brewer.pal.info[select_alpha_color_palette, "maxcolors"]
  base_cols <- RColorBrewer::brewer.pal(min(max_n, max(3, min(max_n, n_cond))), select_alpha_color_palette)
  colors <- colorRampPalette(base_cols)(n_cond)
  if (show_head && show_colors) {
    msg <- paste0(
      "\n--- colors (first few) ---\n",
      paste(utils::capture.output(utils::head(colors, head_n)), collapse = "\n"),
      "\n"
    )
    cat(msg, file = stderr()); flush.console()
  }
  
  
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
  log_head_vec(richness_measures, "richness_measures")
  
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
  
  # Compute and preview richness data
  richness_df <- estimate_richness(physeq1, measures = richness_measures)
  log_head(richness_df, "estimate_richness (output)")
  
  # Add statistical comparisons if specified
  if (pvalue %in% c("Yes", "star")) {
    comps <- combn(conds, 2, simplify = FALSE)
    log_head(comps, "pairwise comparisons (Condition)")
    
    stat_layer <- ggpubr::stat_compare_means(
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
    plot  = richness_plot,
    data1 = richness_df,
    data2 = taxmat
  ))
}
