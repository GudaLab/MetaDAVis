library(ggplot2)
library(tibble)
library(edgeR)
library(dplyr)

edger_summary <- function(
    OTU_input,
    group_index,
    index_pvalue,
    group1,
    group2,
    plot_method,
    alpha,
    edger_color_palette,
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
  
  if (group1 == group2) {
    text_edger_summary <- "Please select different group in condition1 or condition2"
    message(text_edger_summary)
    return(list(edger_result_summary = text_edger_summary))
  }
  
  # Validate palette
  if (!edger_color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    stop("Invalid `edger_color_palette`. Choose a valid RColorBrewer palette.")
  }
  
  # Keep only target groups & align order
  group_index <- subset(group_index, group_index$Condition %in% c(group1, group2))
  log_head(group_index, "group_index (filtered to two groups)")
  
  # Transpose -> samples as rows; align to metadata; back to taxa x samples
  OTU_input <- as.data.frame(t(OTU_input))
  OTU_input <- tibble::rownames_to_column(OTU_input, "Samples")
  OTU_input <- OTU_input[OTU_input$Samples %in% group_index$Samples, ]
  OTU_input <- OTU_input[match(group_index$Samples, OTU_input$Samples), , drop = FALSE]
  OTU_input <- OTU_input %>% tibble::remove_rownames() %>% tibble::column_to_rownames("Samples")
  OTU_input <- as.data.frame(t(OTU_input))
  OTU_input <- OTU_input[rowSums(OTU_input[]) > 0, , drop = FALSE]
  log_head(OTU_input, "OTU_input (aligned & nonzero taxa)")
  
  # Per-sample totals
  parent_seq_count <- as.data.frame(colSums(OTU_input))
  parent_seq_count <- tibble::rownames_to_column(parent_seq_count, "Samples")
  colnames(parent_seq_count)[2] <- "Total_counts"
  log_head(parent_seq_count, "parent_seq_count (per-sample totals)")
  
  # Design
  setupinfo <- data.frame(sample_ID = group_index$Samples, Condition = group_index$Condition)
  setupinfo$Condition <- factor(setupinfo$Condition, levels = c(group1, group2))
  rownames(setupinfo) <- setupinfo$sample_ID
  log_head(setupinfo, "setupinfo (design)")
  
  # Condition colors
  unique_conditions <- unique(setupinfo$Condition)
  n_cond <- length(unique_conditions)
  max_n  <- RColorBrewer::brewer.pal.info[edger_color_palette, "maxcolors"]
  base_cols <- RColorBrewer::brewer.pal(min(max_n, max(3, n_cond)), edger_color_palette)
  cond_colors <- colorRampPalette(base_cols)(n_cond)
  if (show_head && show_colors) {
    msg <- paste0(
      "\n--- colors (first few) ---\n",
      paste(utils::capture.output(utils::head(colors, head_n)), collapse = "\n"),
      "\n"
    )
    cat(msg, file = stderr()); flush.console()
  }
  log_head_vec(unique_conditions, "Condition levels")
  
  # edgeR (robust GLM)
  design <- model.matrix(~ factor(setupinfo$Condition, levels = c(group1, group2)))
  dge <- edgeR::DGEList(counts = OTU_input, group = setupinfo$Condition)
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateGLMRobustDisp(dge, design = design)
  fit <- edgeR::glmFit(dge, design = design)
  lrt <- edgeR::glmLRT(fit, coef = 2)
  
  results <- edgeR::topTags(lrt, n = Inf)$table
  results$foldChange <- 2 ^ results$logFC
  res_tax <- cbind(as.data.frame(results), OTU = rownames(results))
  res_tax <- res_tax[, c("OTU", "foldChange", "logFC", "PValue", "FDR")]
  colnames(res_tax) <- c("OTU", "foldChange", "log2FC", "PValue", "FDR")
  rownames(res_tax) <- res_tax$OTU
  res_tax <- res_tax[order(res_tax$FDR, res_tax$PValue), ]
  log_head(res_tax, "res_tax (edgeR results)")
  
  # Significance filter
  fold <- 1  # |log2FC| > 1 (~2-fold)
  if (identical(index_pvalue, "padj")) {
    res_tax_sig <- subset(res_tax, FDR < alpha & fold < abs(log2FC))
  } else if (identical(index_pvalue, "pvalue")) {
    res_tax_sig <- subset(res_tax, PValue < alpha & fold < abs(log2FC))
  } else {
    stop("`index_pvalue` must be 'padj' or 'pvalue'.")
  }
  log_head(res_tax_sig, "res_tax_sig (significant)")
  
  # Mark significance
  res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig), "Yes", "No")
  
  # Build plotting data (if significant)
  df <- NULL
  df1 <- NULL
  if (nrow(res_tax_sig) >= 1) {
    # Heatmap matrix (z-scored counts of significant taxa)
    OTU_input1 <- OTU_input[rownames(OTU_input) %in% rownames(res_tax_sig), , drop = FALSE]
    df1 <- scale(OTU_input1)
    
    # Log-relative normalization for box/violin plots
    abund_table <- t(OTU_input) # samples x taxa
    data_lr <- log((abund_table + 1) / (rowSums(abund_table) + ncol(abund_table)))
    data_lr <- as.data.frame(data_lr)
    
    groups <- setupinfo$Condition
    for (i in res_tax[rownames(res_tax_sig), "OTU"]) {
      label_col <- if (identical(index_pvalue, "padj")) "FDR" else "PValue"
      tmp <- data.frame(
        Value     = data_lr[, i],
        Condition = groups,
        Taxa1     = i,
        Taxa      = paste(i, label_col, "=", round(res_tax[i, label_col], 4)),
        check.names = FALSE
      )
      df <- if (is.null(df)) tmp else rbind(df, tmp)
    }
    log_head(df, "df (plotting table for significant taxa)")
  }
  
  # Summary text
  text_edger_summary <- if (exists("res_tax_sig") && nrow(res_tax_sig) >= 1) {
    paste0("Total of ", nrow(res_tax_sig), " taxa were identified as significant.")
  } else if (identical(index_pvalue, "padj")) {
    "No significant taxa were identified. Adjust the FDR value or change the test correction method."
  } else {
    "No significant taxa were identified. Adjust the P-value value or change the test correction method."
  }
  
  # Plots
  p1 <- NULL
  if (exists("res_tax_sig") && nrow(res_tax_sig) >= 1) {
    if (plot_method == "1") {
      p1 <- ggplot(df, aes(x = Taxa1, y = Value, color = Condition)) +
        geom_boxplot() +
        ggtitle("Significant Taxa") +
        xlab("Taxa") + ylab("Log-relative normalized") +
        theme_bw() +
        theme(legend.position = "right",
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_color_manual(values = cond_colors)
    } else if (plot_method == "2") {
      p1 <- ggplot(df, aes(Condition, Value, colour = Condition)) +
        ylab("Log-relative normalized") +
        geom_boxplot() + geom_jitter(width = 0.2, alpha = 0.7) + theme_bw() +
        facet_wrap(~ Taxa, scales = "free", ncol = 3) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_color_manual(values = cond_colors)
    } else if (plot_method == "3") {
      # Heatmap with properly named annotation colors aligned to columns
      annot_vec <- setNames(setupinfo$Condition, rownames(setupinfo))[colnames(df1)]
      conds     <- unique(annot_vec)
      max_n   <- RColorBrewer::brewer.pal.info[edger_color_palette, "maxcolors"]
      base_cs <- RColorBrewer::brewer.pal(min(max_n, max(3, length(conds))), edger_color_palette)
      ann_cols <- setNames(colorRampPalette(base_cs)(length(conds)), conds)
      
      ha <- ComplexHeatmap::HeatmapAnnotation(
        Group = annot_vec,
        col   = list(Group = ann_cols)
      )
      hm_cols <- circlize::colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick"))
      
      p3 <- ComplexHeatmap::Heatmap(
        df1, name = "Z",
        show_row_names = TRUE,  row_names_gp = grid::gpar(fontsize = 7),
        show_column_names = TRUE, column_names_gp = grid::gpar(fontsize = 7),
        show_row_dend = TRUE, show_column_dend = FALSE, cluster_columns = FALSE,
        top_annotation = ha, col = hm_cols
      )
      p1 <- ggplotify::as.ggplot(p3)
    }
  }
  
  return(list(
    edger_result_summary     = text_edger_summary,
    edger_result_significant = if (exists("res_tax_sig")) res_tax_sig else NULL,
    edger_result_table       = res_tax,
    total_counts             = parent_seq_count,
    plot                     = p1
  ))
}
