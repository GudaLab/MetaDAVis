library(ggpubr)
library(GGally)
library(ggplot2)
library(dplyr)

samples_based_correlation_plot_table <- function(OTU_input, group_index, method, labe_size, select_sample_geom_shape) {
  # Generate the correlation plot
  p <- ggcorr(
    OTU_input,
    label = FALSE,
    label_alpha = TRUE,
    size = labe_size,
    geom = select_sample_geom_shape,
    method = c("pairwise", method)
  )
  
   
  p <- p +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12, color = "black"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  # Round numeric values in the plot data
  p$data <- p$data %>% mutate(across(where(is.numeric), ~ round(., 2)))
  
  # Return the plot and data
  return(list(plot = p, data = p$data))
}
