library(ggfortify)

pca_plot <- function(OTU_input, group_index, pca_label1, pca_label_size, pca_frame, select_pca_color_palette) {
  # Remove rows with zero sums
  OTU_input <- OTU_input[rowSums(OTU_input) > 0, ]
  
  # Perform PCA
  pcDat <- prcomp(t(OTU_input))
  
  # Define theme for the plot
  theme <- theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "line")
  )
  
  # Create dynamic colors
  unique_conditions <- unique(group_index$Condition)
  colors <- colorRampPalette(brewer.pal(9, select_pca_color_palette))(length(unique_conditions))
  
  # Generate PCA plot based on label settings
  if (pca_label1 == TRUE) {
    plots1 <- autoplot(
      pcDat,
      data = group_index,
      colour = "Condition",
      size = 3,  # Adjusted size to be consistent
      label = TRUE,
      label.size = pca_label_size,
      frame = pca_frame,
      frame.type = 'norm'
    ) +
      coord_equal(ratio = 1) +
      theme +
      scale_color_manual(values = colors)
  } else {
    plots1 <- autoplot(
      pcDat,
      data = group_index,
      colour = "Condition",
      size = 3,  # Adjusted size to be consistent
      label = FALSE,
      label.size = pca_label_size,
      frame = pca_frame,
      frame.type = 'norm'
    ) +
      coord_equal(ratio = 1) +
      theme +
      scale_color_manual(values = colors)
  }
  
  # Extract PCA data
  pca_data <- as.data.frame(plots1$data[, c(1:3)])
  pca_data <- pca_data %>% mutate(across(where(is.numeric), ~ round(., 2)))
  
  # Return the plot and PCA data
  return(list(plot = plots1, data = pca_data))
}
