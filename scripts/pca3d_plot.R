library("plotly")
pca3d_plot <- function(OTU_input, group_index, select_pca3d_color_palette) {
  # Transpose and convert input data
  OTU_input <- as.data.frame(t(OTU_input))
  
  # Perform PCA with 3 components
  prin_comp <- prcomp(OTU_input, rank. = 3)
  components <- prin_comp[["x"]]
  components <- data.frame(components)
  
  # Invert PC2 and PC3 for better visualization
  components$PC2 <- -components$PC2
  components$PC3 <- -components$PC3
  
  # Merge with group index
  components <- merge(components, group_index, by = 'row.names', all = TRUE)
  colnames(components)[1] <- "Samples"
  
  # Round numeric columns for cleaner data
  components <- components %>% mutate(across(where(is.numeric), ~ round(., 2)))
  
  # Generate colors for unique conditions
  unique_conditions <- unique(components$Condition)
  colors <- colorRampPalette(brewer.pal(9, select_pca3d_color_palette))(length(unique_conditions))
  
  # Create a plotly 3D scatter plot
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
      legend = list(title = list(text = 'Condition')),
      scene = list(
        xaxis = list(title = "PC1"),
        yaxis = list(title = "PC2"),
        zaxis = list(title = "PC3"),
        bgcolor = "white"
      )
    )
  
  # Return the plot and data
  return(list(plot = fig, data = components))
}