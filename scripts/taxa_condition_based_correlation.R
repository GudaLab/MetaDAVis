library(ggpubr)
library(GGally)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

taxa_condition_based_correlation_plot_table <- function(OTU_input, group_index, group1, method, labe_size, select_taxa_condition_geom_shape) {
  # Transpose OTU data for taxa-based correlation

  
  group_index <- subset(group_index, (group_index$Condition %in% group1))
  
  OTU_input <- as.data.frame(t(OTU_input))
  #convert row name to 1st column
  OTU_input <- tibble::rownames_to_column(OTU_input, "Samples")
  #mapping column based on metadata
  OTU_input <- OTU_input[OTU_input$Samples %in% group_index$Samples, ]
  #convert 1st column to row name
  OTU_input <- OTU_input %>% remove_rownames %>% column_to_rownames(var="Samples")
  
  OTU_input <- as.data.frame((OTU_input))
  
  # Generate correlation plot
  p <- ggcorr(
    OTU_input,
    label = FALSE,
    label_alpha = TRUE,
    size = labe_size,
    geom = select_taxa_condition_geom_shape,
    method = c("pairwise", method)
  )
  
  # Add color scale and improve aesthetics
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
