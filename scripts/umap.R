library(phyloseq)
library(bluster)
library(patchwork)
library(scater)
library(mia)
library(RColorBrewer)
library(dplyr)

umap_plot_table <- function(
    OTU_input,
    group_index,
    tax_index,
    method,
    kvalue,
    select_umap_color_palette,
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
  
  # Validate palette name
  if (!select_umap_color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    stop("Invalid color palette. Choose a valid name from RColorBrewer.")
  }
  
  # Generate dynamic colors (handle > maxcolors)
  conds <- unique(sample_data(physeq1)$Condition)
  n_cond <- length(conds)
  max_n <- RColorBrewer::brewer.pal.info[select_umap_color_palette, "maxcolors"]
  base_cols <- RColorBrewer::brewer.pal(min(max_n, max(3, min(max_n, n_cond))), select_umap_color_palette)
  colors <- colorRampPalette(base_cols)(n_cond)
  if (show_head && show_colors) {
    cat("\n--- colors (first few) ---\n"); print(utils::head(colors, head_n))
  }
  log_head_vec(conds, "Condition levels")
  
  # Convert to TreeSummarizedExperiment
  tse <- makeTreeSummarizedExperimentFromPhyloseq(physeq1)
  
  # Ensure kvalue is numeric and valid
  kvalue <- as.numeric(kvalue)
  if (is.na(kvalue) || kvalue < 2) stop("`kvalue` must be a number >= 2.")
  
  # Validate the selected method
  valid_methods <- c("counts", "rclr", "hellinger", "pa", "rank", "relabundance")
  if (!method %in% valid_methods) {
    stop("Invalid method. Please select one of: ", paste(valid_methods, collapse = ", "))
  }
  
  # Apply the transformation and run UMAP
  if (method != "counts") {
    tse <- transformCounts(tse, method = method)
  }
  expr_slot <- if (method == "counts") "counts" else method
  tse <- runUMAP(tse, name = "UMAP", exprs_values = expr_slot)
  
  # Peek at UMAP coordinates
  if ("UMAP" %in% reducedDimNames(tse)) {
    log_head(reducedDim(tse, "UMAP"), "UMAP coordinates (from reducedDim)")
  }
  
  # Cluster rows based on kvalue
  emb_mat <- assays(tse)[[expr_slot]]
  log_head(t(emb_mat), paste0("assay '", expr_slot, "' (samples x features, preview)"))
  graph_clusters <- clusterRows(t(assays(tse)[[expr_slot]]), NNGraphParam(k = kvalue))
  
  # Create UMAP plots
  plots <- plotUMAP(tse, colour_by = "Condition") +
    labs(title = expr_slot) +
    scale_color_manual(values = colors)
  
  plots1 <- plotUMAP(tse, colour_by = I(graph_clusters)) +
    labs(title = paste0("k = ", kvalue))
  
  # Combine plots
  umapplot <- plots + plots1
  
  # Round numeric values in plot data
  if (!is.null(plots$data))  { plots$data  <- plots$data  %>% mutate(across(where(is.numeric), ~ round(., 2))) }
  if (!is.null(plots1$data)) { plots1$data <- plots1$data %>% mutate(across(where(is.numeric), ~ round(., 2))) }
  log_head(plots$data,  "plots$data (rounded)")
  log_head(plots1$data, "plots1$data (rounded)")
  
  # Return results
  return(list(plot = umapplot, data = plots$data, data1 = plots1$data))
}
