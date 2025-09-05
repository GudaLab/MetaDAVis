library(phyloseq)
library(bluster)
library(patchwork)
library(scater)
library(mia)
library(RColorBrewer)
library(dplyr)

tsne_plot_table <- function(
    OTU_input,
    group_index,
    tax_index,
    method,
    dimension,
    select_tsne_color_palette,
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
  
  # Prepare taxonomic data
  colnames(tax_index)[1] <- "Species"
  rownames(tax_index) <- tax_index[, 1]
  taxmat <- as.matrix(tax_index)
  log_head(tax_index, "tax_index (with Species col & rownames)")
  
  # Prepare phyloseq object
  OTU <- otu_table(OTU_input, taxa_are_rows = TRUE)
  TAX <- tax_table(taxmat)
  sampledata <- sample_data(group_index)
  physeq1 <- phyloseq(OTU, TAX, sampledata)
  
  log_head(OTU_input,   "OTU_input (raw)")
  log_head(group_index, "group_index (sample metadata)")
  
  # Validate palette
  if (!select_tsne_color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    stop("Invalid color palette. Choose a valid palette from RColorBrewer.")
  }
  
  # Generate dynamic colors
  unique_conditions <- unique(sample_data(physeq1)$Condition)
  max_n <- RColorBrewer::brewer.pal.info[select_tsne_color_palette, "maxcolors"]
  base_cols <- RColorBrewer::brewer.pal(min(max_n, max(3, min(max_n, length(unique_conditions)))), select_tsne_color_palette)
  colors <- colorRampPalette(base_cols)(length(unique_conditions))
  if (show_head && show_colors) {
    cat("\n--- colors (first few) ---\n"); print(utils::head(colors, head_n))
  }
  log_head_vec(unique_conditions, "Condition levels")
  
  # Convert to TreeSummarizedExperiment
  tse <- makeTreeSummarizedExperimentFromPhyloseq(physeq1)
  
  # Ensure dimension is numeric
  dimension <- as.numeric(dimension)
  if (is.na(dimension) || dimension < 2) stop("`dimension` must be a number >= 2.")
  
  # Transform counts and run TSNE
  if (method == "counts") {
    tse <- runTSNE(tse, name = "TSNE", exprs_values = "counts", ncomponents = dimension)
  } else if (method %in% c("rclr", "hellinger", "pa", "rank", "relabundance")) {
    tse <- transformCounts(tse, method = method)
    tse <- runTSNE(tse, name = "TSNE", exprs_values = method, ncomponents = dimension)
  } else {
    stop("Invalid method. Choose: 'counts', 'rclr', 'hellinger', 'pa', 'rank', or 'relabundance'.")
  }
  
  # Peek at reduced coordinates
  if ("TSNE" %in% reducedDimNames(tse)) {
    log_head(reducedDim(tse, "TSNE"), "TSNE coordinates (from reducedDim)")
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
  log_head(plottsne$data, "plottsne$data (rounded)")
  
  # Return the plot and data
  return(list(plot = plottsne, data = plottsne$data))
}
