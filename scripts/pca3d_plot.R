library(plotly)

pca3d_plot <- function(
    OTU_input,
    group_index,
    select_pca3d_color_palette,
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
  
  # 1) Keep only samples present in both the count matrix and metadata.
  shared_samples <- intersect(rownames(group_index), colnames(OTU_input))
  if (length(shared_samples) == 0) {
    stop("No overlapping samples between OTU_input (columns) and group_index (rownames).")
  }

  OTU_input <- OTU_input[, shared_samples, drop = FALSE]
  group_index <- group_index[shared_samples, , drop = FALSE]
  OTU_input <- OTU_input[rowSums(OTU_input) > 0, , drop = FALSE]
  log_head_vec(shared_samples, "shared_samples")
  log_head(OTU_input, "OTU_input (synced and filtered)")

  # 2) PCA with up to 3 components; pad missing PCs so small datasets still render.
  pca_input <- as.data.frame(t(OTU_input))
  rank_components <- min(3, nrow(pca_input), ncol(pca_input))
  if (rank_components < 1) {
    stop("PCA-3D needs at least one non-empty sample and one non-empty taxon.")
  }

  prin_comp <- prcomp(pca_input, rank. = rank_components)
  components <- data.frame(prin_comp$x)
  for (pc_name in c("PC1", "PC2", "PC3")) {
    if (!pc_name %in% names(components)) {
      components[[pc_name]] <- 0
    }
  }
  components <- components[, c("PC1", "PC2", "PC3"), drop = FALSE]
  log_head(components, "PCA scores (PC1-3)")
  
  # 3) Invert PC2 and PC3 for better visualization
  components$PC2 <- -components$PC2
  components$PC3 <- -components$PC3
  
  # 4) Add metadata without creating unmatched rows.
  components$Samples <- rownames(components)
  components <- cbind(components, group_index[rownames(components), , drop = FALSE])
  components <- components[, c("Samples", setdiff(names(components), "Samples")), drop = FALSE]
  log_head(components, "components (with metadata)")
  
  # 5) Round numeric columns
  components <- components %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(., 2)))
  log_head(components, "components (rounded)")
  
  # 6) Generate colors for unique conditions
  unique_conditions <- unique(components$Condition)
  colors <- colorRampPalette(RColorBrewer::brewer.pal(9, select_pca3d_color_palette))(length(unique_conditions))
  if (show_head && show_colors) {
    cat("\n--- colors (first few) ---\n"); print(utils::head(colors, head_n))
  }
  log_head_vec(unique_conditions, "Condition levels")
  
  # 7) 3D PCA scatter
  fig <- plot_ly(
    components,
    x = ~PC1,
    y = ~PC2,
    z = ~PC3,
    text = ~Samples,
    color = ~Condition,
    colors = colors,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 6, opacity = 0.9)
  ) %>%
    layout(
      legend = list(title = list(text = "Condition")),
      scene = list(
        xaxis = list(title = "PC1"),
        yaxis = list(title = "PC2"),
        zaxis = list(title = "PC3"),
        bgcolor = "white"
      )
    )
  
  return(list(plot = fig, data = components))
}
