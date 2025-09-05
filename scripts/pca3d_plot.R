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
  
  # 1) Transpose and convert input data
  OTU_input <- as.data.frame(t(OTU_input))
  log_head(OTU_input, "OTU_input (transposed)")
  
  # 2) PCA with 3 components
  prin_comp <- prcomp(OTU_input, rank. = 3)
  components <- data.frame(prin_comp$x)
  log_head(components, "PCA scores (PC1-3)")
  
  # 3) Invert PC2 and PC3 for better visualization
  components$PC2 <- -components$PC2
  components$PC3 <- -components$PC3
  
  # 4) Merge with metadata
  components <- merge(components, group_index, by = "row.names", all = TRUE)
  colnames(components)[1] <- "Samples"
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
    colors = colors
  ) %>%
    add_markers(size = 12) %>%
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
