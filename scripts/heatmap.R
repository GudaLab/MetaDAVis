library("ComplexHeatmap")
library("scales")
library("ggplotify")
library("circlize")
library("RColorBrewer")

set.seed(123)
heatmap_abundance <- function(
    OTU_input,
    group_index,
    heatmap_clustering_method_rows,
    heatmap_clustering_method_columns,
    heatmap_color_palette,
    heatmap_normalization,
    heatmap_row_names,
    heatmap_row_names_size,
    heatmap_column_names,
    heatmap_column_names_size,
    heatmap_row_dend,
    heatmap_column_dend,
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
  
  heatmap_row_dend    <- as.logical(heatmap_row_dend)
  heatmap_column_dend <- as.logical(heatmap_column_dend)
  
  # filter out empty taxa
  OTU_input <- OTU_input[rowSums(OTU_input[]) > 0, , drop = FALSE]
  log_head(OTU_input, "OTU_input (filtered: row sums > 0)")
  
  # ----- normalization -----
  if (heatmap_normalization == "scale") {
    df1 <- suppressWarnings(scale(OTU_input))  # Z-score across columns
  } else if (heatmap_normalization == "minmax") {
    # safe min-max per column (avoid divide-by-zero)
    df1 <- apply(OTU_input, 2, function(x) {
      rng <- max(x) - min(x)
      if (rng == 0) return(rep(0, length(x)))
      (x - min(x)) / rng
    })
  } else if (heatmap_normalization == "log") {
    df1 <- log(OTU_input + 1)
  } else if (heatmap_normalization == "row") {
    df1 <- t(apply(OTU_input, 1, function(x) if (sum(x) == 0) x else x / sum(x)))
  } else if (heatmap_normalization == "column") {
    df1 <- apply(OTU_input, 2, function(x) if (sum(x) == 0) x else x / sum(x))
  } else {
    df1 <- OTU_input
  }
  df1 <- as.matrix(df1)
  log_head(df1, paste0("df1 (after '", heatmap_normalization, "' normalization)"))
  
  # ----- palette checks & construction -----
  if (!heatmap_color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    stop("Invalid heatmap_color_palette. Choose a valid palette from RColorBrewer.")
  }
  # Breakpoints for color mapping; adapt to data when not scaled
  if (heatmap_normalization %in% c("scale")) {
    breaks <- c(-4, 0, 4)
  } else {
    rng <- range(df1, finite = TRUE)
    mid <- mean(rng)
    breaks <- c(rng[1], mid, rng[2])
    if (!is.finite(mid)) breaks <- c(-4, 0, 4)  # fallback
  }
  base_cols <- RColorBrewer::brewer.pal(3, heatmap_color_palette)
  color_scheme <- circlize::colorRamp2(breaks, base_cols)
  if (show_head && show_colors) {
    cat("\n--- heatmap color breakpoints ---\n"); print(breaks)
    cat("\n--- base colors ---\n"); print(base_cols)
  }
  
  # ----- sample annotation -----
  conds <- unique(group_index$Condition)
  n_cond <- length(conds)
  max_n  <- RColorBrewer::brewer.pal.info[heatmap_color_palette, "maxcolors"]
  anno_base <- RColorBrewer::brewer.pal(min(max_n, max(3, min(max_n, n_cond))), heatmap_color_palette)
  anno_cols <- colorRampPalette(anno_base)(n_cond)
  names(anno_cols) <- conds
  if (show_head && show_colors) {
    msg <- paste0(
      "\n--- colors (first few) ---\n",
      paste(utils::capture.output(utils::head(colors, head_n)), collapse = "\n"),
      "\n"
    )
    cat(msg, file = stderr()); flush.console()
  }
  log_head_vec(conds, "Condition levels (annotation)")
  
  ha <- ComplexHeatmap::HeatmapAnnotation(
    Group = group_index$Condition,
    annotation_height = unit(4, "mm"),
    col = list(Group = anno_cols)
  )
  
  # ----- build heatmap -----
  p <- ComplexHeatmap::Heatmap(
    df1, name = ifelse(heatmap_normalization == "scale", "Z", "Value"),
    show_row_names = heatmap_row_names,
    row_names_gp = grid::gpar(fontsize = heatmap_row_names_size),
    show_column_names = heatmap_column_names,
    column_names_gp = grid::gpar(fontsize = heatmap_column_names_size),
    show_row_dend = heatmap_row_dend,
    show_column_dend = heatmap_column_dend,
    top_annotation = ha,
    col = color_scheme,
    clustering_method_columns = heatmap_clustering_method_columns,
    clustering_method_rows = heatmap_clustering_method_rows
  )
  
  q <- ggplotify::as.ggplot(p)
  
  return(list(plot = p, plot1 = q))
}
