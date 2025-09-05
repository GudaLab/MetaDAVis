# Requires: ggplot2, ggfortify, RColorBrewer, dplyr
library(ggfortify)
pca_plot <- function(
    OTU_input,
    group_index,
    pca_label1,
    pca_label_size,
    pca_frame,
    select_pca_color_palette,
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
  
  # 0) sync samples by intersection of sample IDs (rownames in group_index vs colnames in OTU_input)
  if (is.null(rownames(group_index))) {
    stop("group_index must have rownames that are sample IDs.")
  }
  shared_samples <- intersect(rownames(group_index), colnames(OTU_input))
  if (length(shared_samples) == 0) stop("No overlapping samples between OTU_input (columns) and group_index (rownames).")
  OTU_input  <- OTU_input[, shared_samples, drop = FALSE]
  group_index <- group_index[shared_samples, , drop = FALSE]
  log_head_vec(shared_samples, "shared_samples")
  
  # 1) remove rows (taxa) with zero sums
  OTU_input <- OTU_input[rowSums(OTU_input) > 0, , drop = FALSE]
  log_head(OTU_input, "OTU_input (taxa filtered: row sums > 0)")
  
  # 2) PCA on samples (transpose so samples are rows)
  pcDat <- prcomp(t(OTU_input))  # center=TRUE by default
  log_head(pcDat$x, "PCA scores (samples x PCs)")
  
  # 3) theme (rename to avoid shadowing ggplot2::theme)
  theme_pca <- ggplot2::theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks  = element_line(colour = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "line")
  )
  
  # 4) colors
  if (!select_pca_color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    stop("Invalid color palette. Choose a valid name from RColorBrewer.")
  }
  unique_conditions <- unique(group_index$Condition)
  max_n <- RColorBrewer::brewer.pal.info[select_pca_color_palette, "maxcolors"]
  base_cols <- RColorBrewer::brewer.pal(min(max_n, max(3, min(max_n, length(unique_conditions)))), select_pca_color_palette)
  colors <- colorRampPalette(base_cols)(length(unique_conditions))
  if (show_head && show_colors) {
    cat("\n--- colors (first few) ---\n"); print(utils::head(colors, head_n))
  }
  log_head_vec(unique_conditions, "Condition levels")
  
  # 5) Plot
  plots1 <- autoplot(
    pcDat,
    data  = group_index,
    colour = "Condition",
    size   = 3,
    label  = isTRUE(pca_label1),
    label.size = pca_label_size,
    frame  = pca_frame,
    frame.type = "norm"
  ) +
    coord_equal(ratio = 1) +
    theme_pca +
    scale_color_manual(values = colors)
  
  # 6) Extract plotted data (PC1/PC2 + Condition)
  pdat <- as.data.frame(plots1$data[, 1:3])
  colnames(pdat)[1:3] <- c("PC1","PC2","Condition")  # enforce clear names if needed
  suppressWarnings({
    pdat <- dplyr::mutate(pdat, dplyr::across(where(is.numeric), ~ round(., 2)))
  })
  log_head(pdat, "pca_data (rounded)")
  
  return(list(plot = plots1, data = pdat))
}
