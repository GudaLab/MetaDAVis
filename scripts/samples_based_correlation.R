library(ggpubr)
library(GGally)
library(ggplot2)
library(dplyr)

samples_based_correlation_plot_table <- function(
    OTU_input,
    group_index,
    method,
    labe_size,
    group1,
    select_sample_geom_shape,
    show_head = TRUE,
    head_n = 5
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

  
  # 1) subset metadata to selected groups
  group_index <- subset(group_index, (group_index$Condition %in% group1))
  log_head(group_index, "group_index (filtered to group1)")
  
  # 2) transpose OTU -> samples as rows
  OTU_input <- as.data.frame(t(OTU_input))
  log_head(OTU_input, "OTU_input (transposed; samples as rows)")
  
  # 3) rownames -> first column
  OTU_input <- tibble::rownames_to_column(OTU_input, "Samples")
  log_head(OTU_input, "OTU_input (with 'Samples' column)")
  
  # 4) keep only samples present in metadata
  OTU_input <- OTU_input[OTU_input$Samples %in% group_index$Samples, ]
  log_head(OTU_input, "OTU_input (matched to metadata Samples)")
  
  # 5) first column back to rownames
  OTU_input <- OTU_input %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var = "Samples")
  
  # 6) transpose back so taxa as rows (GGally::ggcorr expects variables as columns)
  OTU_input <- as.data.frame(t(OTU_input))
  log_head(OTU_input, "OTU_input (final for ggcorr)")
  
  # 7) correlation plot
  p <- GGally::ggcorr(
    OTU_input,
    label       = FALSE,
    label_alpha = TRUE,
    size        = labe_size,
    geom        = select_sample_geom_shape,
    method      = c("pairwise", method)
  )
  
  p <- p +
    theme_minimal() +
    theme(
      axis.text      = element_text(size = 12, color = "black"),
      legend.title   = element_text(size = 12),
      legend.text    = element_text(size = 10)
    )
  
  # 8) round numeric values in plot data (if available)
  if (!is.null(p$data)) {
    p$data <- p$data %>% mutate(across(where(is.numeric), ~ round(., 2)))
    log_head(p$data, "p$data (rounded)")
  }
  
  return(list(plot = p, data = p$data))
}
