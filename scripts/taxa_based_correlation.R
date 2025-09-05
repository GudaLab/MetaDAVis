library(ggpubr)
library(GGally)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

taxa_based_correlation_plot_table <- function(
    OTU_input,
    group_index,
    method,
    labe_size,
    select_taxa_geom_shape,
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

  
  # 1) transpose so taxa are columns (variables) and samples are rows
  OTU_input <- as.data.frame(t(OTU_input))
  log_head(OTU_input, "OTU_input (transposed: samples x taxa)")
  
  # 2) correlation plot
  p <- GGally::ggcorr(
    OTU_input,
    label       = FALSE,
    label_alpha = TRUE,
    size        = labe_size,
    geom        = select_taxa_geom_shape,
    method      = c("pairwise", method)
  )
  
  p <- p +
    theme_minimal() +
    theme(
      axis.text     = element_text(size = 12, color = "black"),
      legend.title  = element_text(size = 12),
      legend.text   = element_text(size = 10)
    )
  
  # 3) round numeric values in plot data (if present)
  if (!is.null(p$data)) {
    p$data <- p$data %>% mutate(across(where(is.numeric), ~ round(., 2)))
    log_head(p$data, "p$data (rounded)")
  }
  
  return(list(plot = p, data = p$data))
}
