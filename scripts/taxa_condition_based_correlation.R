# library(ggpubr)
# library(GGally)
# library(ggplot2)
# library(dplyr)

taxa_condition_based_correlation_plot_table <- function(
    OTU_input,
    group_index,
    group1,
    method,
    labe_size,
    select_taxa_condition_geom_shape,
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

  
  # 1) subset metadata to selected conditions
  group_index <- subset(group_index, (group_index$Condition %in% group1))
  log_head(group_index, "group_index (filtered to selected Condition(s))")
  
  # 2) transpose OTU so samples become rows
  OTU_input <- as.data.frame(t(OTU_input))
  log_head(OTU_input, "OTU_input (transposed; samples as rows)")
  
  # 3) move rownames -> first column "Samples"
  OTU_input <- tibble::rownames_to_column(OTU_input, "Samples")
  log_head(OTU_input, "OTU_input (+ 'Samples' column)")
  
  # 4) keep only samples present in filtered metadata
  OTU_input <- OTU_input[OTU_input$Samples %in% group_index$Samples, ]
  log_head(OTU_input, "OTU_input (matched to metadata Samples)")
  
  # 5) first column back to rownames
  OTU_input <- OTU_input %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "Samples")
  
  # 6) ensure it's a plain data.frame (variables = taxa as columns)
  OTU_input <- as.data.frame(OTU_input)
  log_head(OTU_input, "OTU_input (final for ggcorr)")
  
  # 7) correlation plot
  p <- GGally::ggcorr(
    OTU_input,
    label       = FALSE,
    label_alpha = TRUE,
    size        = labe_size,
    geom        = select_taxa_condition_geom_shape,
    method      = c("pairwise", method)
  ) +
    theme_minimal() +
    theme(
      axis.text     = element_text(size = 12, color = "black"),
      legend.title  = element_text(size = 12),
      legend.text   = element_text(size = 10)
    )
  
  # 8) round numeric values in plot data (if present)
  if (!is.null(p$data)) {
    p$data <- p$data %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(., 2)))
    log_head(p$data, "p$data (rounded)")
  }
  
  return(list(plot = p, data = p$data))
}
