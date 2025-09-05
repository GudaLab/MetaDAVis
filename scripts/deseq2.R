library("DESeq2")
library("ggplot2")
library("tibble")
library("qvalue")
library("RColorBrewer")
library("ComplexHeatmap")
library("ggplotify")
library("dplyr")

deseq2_summary <- function(
    OTU_input,
    group_index,
    index_pvalue,
    group1,
    group2,
    plot_method,
    alpha,
    deseq2_color_palette,
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
  
  # early check
  if (group1 == group2) {
    text_deseq_summary <- "Please select different group in condition1 or condition2"
    message(text_deseq_summary)
    return(list(deseq2_result_summary = text_deseq_summary))
  }
  
  # keep only target groups
  group_index <- subset(group_index, (group_index$Condition %in% c(group1, group2)))
  log_head(group_index, "group_index (filtered to two groups)")
  
  # align OTU_input columns to metadata (order-preserving)
  OTU_input <- as.data.frame(t(OTU_input))
  OTU_input <- tibble::rownames_to_column(OTU_input, "Samples")
  OTU_input <- OTU_input[OTU_input$Samples %in% group_index$Samples, ]
  OTU_input <- OTU_input[match(group_index$Samples, OTU_input$Samples), , drop = FALSE]
  OTU_input <- OTU_input %>% tibble::remove_rownames() %>% tibble::column_to_rownames("Samples")
  OTU_input <- as.data.frame(t(OTU_input))
  OTU_input <- OTU_input[rowSums(OTU_input[]) > 0, , drop = FALSE]
  log_head(OTU_input, "OTU_input (aligned & filtered taxa)")
  
  # totals
  parent_seq_count <- as.data.frame(colSums(OTU_input))
  parent_seq_count <- tibble::rownames_to_column(parent_seq_count, "Samples")
  colnames(parent_seq_count)[2] <- "Total_counts"
  log_head(parent_seq_count, "parent_seq_count (per-sample totals)")
  
  # design data
  setupinfo <- data.frame(sample_ID = group_index$Samples, Condition = group_index$Condition)
  setupinfo$Condition <- factor(setupinfo$Condition, levels = c(group1, group2))
  rownames(setupinfo) <- setupinfo$sample_ID
  log_head(setupinfo, "setupinfo (design)")
  
  # ensure DESeq2 matrix matches colData ordering
  OTU_input <- OTU_input[, rownames(setupinfo), drop = FALSE]
  
  # DESeq2 pipeline
  dds <- DESeqDataSetFromMatrix(countData = round(OTU_input),
                                colData = setupinfo,
                                design = ~ Condition)
  dds <- dds[, dds$Condition %in% c(group1, group2)]
  dds <- DESeq(dds, parallel = FALSE)
  
  res <- DESeq2::results(dds, alpha = alpha)
  res_tax <- as.data.frame(res)
  res_tax <- cbind(res_tax, OTU = rownames(res_tax))
  res_tax$foldChange <- 2 ^ abs(res_tax$log2FoldChange)
  res_tax <- res_tax[, c("OTU", "baseMean", "foldChange", "log2FoldChange", "pvalue", "padj")]
  colnames(res_tax) <- c("OTU", "Basemean", "foldChange", "log2FC", "PValue", "FDR")
  rownames(res_tax) <- res_tax$OTU
  res_tax <- res_tax[order(res_tax$FDR, res_tax$PValue), ]
  log_head(res_tax, "res_tax (DESeq2 results)")
  
  normalized.counts <- as.data.frame(counts(dds, normalized = TRUE))
  log_head(normalized.counts, "normalized.counts")
  
  # significance filter
  fold <- 1  # |log2FC| > 1 (~2-fold)
  if (identical(index_pvalue, "padj")) {
    res_tax_sig <- subset(res_tax, FDR < alpha & fold < abs(log2FC))
  } else if (identical(index_pvalue, "pvalue")) {
    res_tax_sig <- subset(res_tax, PValue < alpha & ((log2FC >= 1) | (log2FC <= -1)))
  } else {
    stop("`index_pvalue` must be 'padj' or 'pvalue'.")
  }
  log_head(res_tax_sig, "res_tax_sig (significant)")
  
  # mark significance
  res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig), "Yes", "No")
  
  # plotting helpers
  # log-relative normalized table for box/violin plots
  abund_table <- t(OTU_input)
  data_lr <- log((abund_table + 1) / (rowSums(abund_table) + ncol(abund_table)))
  data_lr <- as.data.frame(data_lr)
  
  df <- NULL
  df1 <- NULL
  if (nrow(res_tax_sig) >= 1) {
    # heatmap matrix (z-scored counts of significant taxa)
    OTU_input_sig <- OTU_input[rownames(OTU_input) %in% rownames(res_tax_sig), , drop = FALSE]
    df1 <- scale(OTU_input_sig)
    
    groups <- group_index$Condition
    for (i in res_tax[rownames(res_tax_sig), "OTU"]) {
      label_col <- if (identical(index_pvalue, "padj")) "FDR" else "PValue"
      tmp <- data.frame(
        Value = data_lr[, i],
        Condition = groups,
        Taxa1 = i,
        Taxa = paste(i, label_col, "=", round(res_tax[i, label_col], 4)),
        check.names = FALSE
      )
      df <- if (is.null(df)) tmp else rbind(df, tmp)
    }
    log_head(df, "df (plotting table for significant taxa)")
  }
  
  # colors
  if (!deseq2_color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    stop("Invalid `deseq2_color_palette`. Choose a valid RColorBrewer palette.")
  }
  unique_conditions <- unique(group_index$Condition)
  n_cond <- length(unique_conditions)
  max_n  <- RColorBrewer::brewer.pal.info[deseq2_color_palette, "maxcolors"]
  base_cols <- RColorBrewer::brewer.pal(min(max_n, max(3, min(max_n, n_cond))), deseq2_color_palette)
  colors <- colorRampPalette(base_cols)(n_cond)
  if (show_head && show_colors) {
    msg <- paste0(
      "\n--- colors (first few) ---\n",
      paste(utils::capture.output(utils::head(colors, head_n)), collapse = "\n"),
      "\n"
    )
    cat(msg, file = stderr()); flush.console()
  }
  log_head_vec(unique_conditions, "Condition levels")
  
  # plots
  p1 <- NULL
  if (nrow(res_tax_sig) >= 1) {
    if (plot_method == "1") {
      p1 <- ggplot(df, aes(x = Taxa1, y = Value, color = Condition)) +
        geom_boxplot() +
        ggtitle("Significant Taxa") +
        xlab("Taxa") + ylab("Log-relative normalized") +
        theme_bw() +
        theme(legend.position = "right",
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_color_manual(values = colors)
    } else if (plot_method == "2") {
      p1 <- ggplot(df, aes(Condition, Value, colour = Condition)) +
        ylab("Log-relative normalized") +
        geom_boxplot() + geom_jitter(width = 0.2, alpha = 0.7) + theme_bw() +
        facet_wrap(~ Taxa, scales = "free", ncol = 3) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_color_manual(values = colors)
    } else if (plot_method == "3") {
      p1 <- ggplot(res_tax, aes(x = Basemean, y = log2FC, color = Significant)) +
        geom_point(size = 2) +
        scale_x_log10() +
        scale_color_manual(values = c("black", "red")) +
        labs(x = "log10(Mean abundance)", y = "Log2FC") +
        theme_bw() +
        geom_text(data = subset(res_tax, Significant == "Yes"),
                  aes(label = OTU), size = 3.5, vjust = 1)
    } else if (plot_method == "4") {
      # Align annotation vector to the heatmap columns
      annot_vec <- setNames(group_index$Condition, group_index$Samples)[colnames(df1)]
      conds     <- unique(annot_vec)
      
      # Build a named color map: names must be the condition levels
      max_n   <- RColorBrewer::brewer.pal.info[deseq2_color_palette, "maxcolors"]
      base_cs <- RColorBrewer::brewer.pal(min(max_n, max(3, length(conds))), deseq2_color_palette)
      ann_cols <- setNames(colorRampPalette(base_cs)(length(conds)), conds)
      
      ha <- ComplexHeatmap::HeatmapAnnotation(
        Group = annot_vec,
        col   = list(Group = ann_cols)  # <-- named vector!
      )
      
      # Continuous color mapping for the matrix
      hm_cols <- circlize::colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick"))
      
      p4 <- ComplexHeatmap::Heatmap(
        df1, name = "Z",
        show_row_names = TRUE,  row_names_gp = grid::gpar(fontsize = 7),
        show_column_names = TRUE, column_names_gp = grid::gpar(fontsize = 7),
        show_row_dend = TRUE, show_column_dend = FALSE, cluster_columns = FALSE,
        top_annotation = ha, col = hm_cols
      )
      p1 <- ggplotify::as.ggplot(p4)
    }
    
  }
  
  # summary text
  text_deseq_summary <- if (exists("res_tax_sig") && nrow(res_tax_sig) >= 1) {
    paste("Total of ", nrow(res_tax_sig), " taxa were identified as significant.", sep = "")
  } else if (identical(index_pvalue, "padj")) {
    "No significant taxa were identified. Adjust the FDR value or change the test correction method."
  } else {
    "No significant taxa were identified. Adjust the P-value value or change the test correction method."
  }
  
  return(list(
    deseq2_result_summary   = text_deseq_summary,
    deseq2_result_significant = if (exists("res_tax_sig")) res_tax_sig else NULL,
    deseq2_result_table     = res_tax,
    deseq2_normalized       = normalized.counts,
    total_counts            = parent_seq_count,
    plot                    = p1
  ))
}
