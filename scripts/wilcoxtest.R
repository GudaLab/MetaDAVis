library("ggplot2")
library("tibble")
library("qvalue")
library("ComplexHeatmap")
library("scales")
library("ggplotify")
library("RColorBrewer")
library("dplyr")

wilcoxtest_summary <- function(
    OTU_input,
    group_index,
    index_pvalue,
    group1,
    group2,
    plot_method,
    alpha,
    wilcox_color_palette,
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
  
  # 0) quick checks
  if (group1 == group2) {
    text_wilcoxtest_summary <- "Please select different group in condition1 or condition2"
    message(text_wilcoxtest_summary)
    return(list(wilcoxtest_result_summary = text_wilcoxtest_summary))
  }
  
  # 1) keep only the two groups of interest
  group_index <- subset(group_index, (group_index$Condition %in% c(group1, group2)))
  log_head(group_index, "group_index (filtered to two groups)")
  
  # 2) map OTU_input columns to samples in metadata (and keep order by metadata)
  OTU_input <- as.data.frame(t(OTU_input))
  OTU_input <- tibble::rownames_to_column(OTU_input, "Samples")
  OTU_input <- OTU_input[OTU_input$Samples %in% group_index$Samples, ]
  OTU_input <- OTU_input[match(group_index$Samples, OTU_input$Samples), , drop = FALSE]  # align order
  OTU_input <- OTU_input %>% tibble::remove_rownames() %>% tibble::column_to_rownames("Samples")
  OTU_input <- as.data.frame(t(OTU_input))
  OTU_input  <- OTU_input[rowSums(OTU_input[]) > 0, , drop = FALSE]
  OTU_input1 <- OTU_input
  log_head(OTU_input, "OTU_input (aligned & filtered taxa)")
  
  # 3) per-sample totals
  parent_seq_count <- as.data.frame(colSums(OTU_input))
  parent_seq_count <- tibble::rownames_to_column(parent_seq_count, "Samples")
  colnames(parent_seq_count)[2] <- "Total_counts"
  log_head(parent_seq_count, "parent_seq_count (per-sample totals)")
  
  # 4) relative frequencies (%)
  rel_freq <- sweep(OTU_input, 2, colSums(OTU_input), "/") * 100
  rel_freq[!is.finite(as.matrix(rel_freq))] <- 0
  OTU_input <- as.data.frame(rel_freq)
  log_head(OTU_input, "OTU_input (relative frequency %)")
  
  # 5) prep groups as a vector aligned to sample order
  group_index1 <- group_index
  group_vec <- setNames(group_index1$Condition, group_index1$Samples)
  group_vec <- group_vec[colnames(OTU_input)]
  log_head_vec(group_vec, "group_vec (aligned to samples)")
  
  # 6) wilcoxon tests (per taxon)
  g_names <- names(table(group_vec))
  test_summary <- matrix(NA, nrow = nrow(OTU_input), ncol = 11)
  colnames(test_summary) <- c(
    "OTU",
    paste("Present_in_no_of", g_names, sep = "_"),
    paste("Mean_relative_frequency", g_names, sep = "_"),
    "All_mean_relative_frequency", "Difference_between_means",
    "fold_change", "log2FC", "PValue", "q_value"
  )
  
  for (i in 1:nrow(OTU_input)) {
    OTU_tmp <- OTU_input[i, ]
    g1_vals <- as.numeric(OTU_tmp[group_vec == g_names[1]])
    g2_vals <- as.numeric(OTU_tmp[group_vec == g_names[2]])
    wt <- suppressWarnings(wilcox.test(g1_vals, g2_vals))
    
    mean1 <- mean(g1_vals, na.rm = TRUE)
    mean2 <- mean(g2_vals, na.rm = TRUE)
    
    test_summary[i, paste("Mean_relative_frequency", g_names[1], sep = "_")] <- signif(mean1, 4)
    test_summary[i, paste("Mean_relative_frequency", g_names[2], sep = "_")] <- signif(mean2, 4)
    test_summary[i, "All_mean_relative_frequency"] <- signif((mean1 + mean2) / 2, 4)
    test_summary[i, "Difference_between_means"] <- signif(mean1 - mean2, 4)
    # ratio mean2/mean1; protect divide-by-zero
    ratio <- if (mean1 == 0) Inf else mean2 / mean1
    test_summary[i, "fold_change"] <- signif(ratio, 4)
    test_summary[i, "log2FC"] <- signif(log2(if (is.finite(ratio) && ratio > 0) ratio else NA), 4)
    test_summary[i, paste("Present_in_no_of", g_names[1], sep = "_")] <- sum(g1_vals != 0, na.rm = TRUE)
    test_summary[i, paste("Present_in_no_of", g_names[2], sep = "_")] <- sum(g2_vals != 0, na.rm = TRUE)
    test_summary[i, "OTU"] <- rownames(OTU_input)[i]
    test_summary[i, "PValue"] <- signif(wt$p.value, 4)
  }
  
  test_summary2 <- test_summary[!is.na(test_summary[, "PValue"]), , drop = FALSE]
  test_summary2[, "q_value"] <- signif(p.adjust(as.numeric(test_summary2[, "PValue"]), "BH"), 4)
  
  res_tax <- as.data.frame(test_summary2, stringsAsFactors = FALSE)
  rownames(res_tax) <- res_tax[, "OTU"]
  # handle non-finite log2FC
  suppressWarnings({
    res_tax$log2FC <- as.numeric(res_tax$log2FC)
    res_tax$PValue <- as.numeric(res_tax$PValue)
    res_tax$q_value <- as.numeric(res_tax$q_value)
  })
  res_tax$log2FC[!is.finite(res_tax$log2FC)] <- 0
  log_head(res_tax, "res_tax (test summary)")
  
  # 7) significance filtering
  fold <- 1  # threshold on |log2FC|; 1 ≈ 2-fold
  if (identical(index_pvalue, "padj")) {
    res_tax_sig <- subset(res_tax, q_value < alpha & fold < abs(log2FC))
  } else if (identical(index_pvalue, "pvalue")) {
    res_tax_sig <- subset(res_tax, PValue < alpha & fold < abs(log2FC))
  } else {
    stop("`index_pvalue` must be 'padj' or 'pvalue'.")
  }
  log_head(res_tax_sig, "res_tax_sig (significant)")
  
  # mark significance
  res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig), "Yes", "No")
  
  # 8) palette for conditions
  if (!wilcox_color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    stop("Invalid `wilcox_color_palette`. Choose a valid RColorBrewer palette.")
  }
  unique_conditions <- unique(group_index1$Condition)
  n_cond <- length(unique_conditions)
  max_n  <- RColorBrewer::brewer.pal.info[wilcox_color_palette, "maxcolors"]
  base_cols <- RColorBrewer::brewer.pal(min(max_n, max(3, min(max_n, n_cond))), wilcox_color_palette)
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
  
  # 9) build long df for plots 1/2 if significant
  df <- NULL
  if (nrow(res_tax_sig) >= 1) {
    res_tax$q_value <- as.numeric(res_tax$q_value)
    data_long <- as.data.frame(t(OTU_input))  # samples x taxa
    # align condition vector to rows of data_long
    cond_vec <- setNames(group_index1$Condition, group_index1$Samples)[rownames(data_long)]
    for (i in res_tax[rownames(res_tax_sig), "OTU"]) {
      tmp <- data.frame(
        Value = data_long[, i],
        Condition = cond_vec,
        Taxa1 = i,
        Taxa = if (identical(index_pvalue, "padj")) {
          paste0(i, " q_value = ", round(res_tax[i, "q_value"], 4))
        } else {
          paste0(i, " PValue = ", round(res_tax[i, "PValue"], 4))
        },
        row.names = NULL,
        check.names = FALSE
      )
      df <- if (is.null(df)) tmp else rbind(df, tmp)
    }
    log_head(df, "df (plotting table for significant taxa)")
  }
  
  # 10) text summary + plots
  if (nrow(res_tax_sig) >= 1) {
    text_wilcoxtest_summary <- paste("Total of ", nrow(res_tax_sig), " taxa were identified as significant.", sep = "")
  } else {
    text_wilcoxtest_summary <- if (identical(index_pvalue, "padj")) {
      "No significant taxa were identified. Adjust FDR (alpha) or change the correction method."
    } else {
      "No significant taxa were identified. Adjust P-value (alpha) or change the correction method."
    }
  }
  
  p1 <- NULL
  if (nrow(res_tax_sig) >= 1) {
    if (plot_method == "1") {
      p1 <- ggplot(df, aes(x = Taxa1, y = Value, color = Condition)) +
        geom_boxplot() +
        scale_y_log10() +
        ggtitle("Significant Taxa") +
        xlab("Taxa") + ylab("log10(Relative frequency)") +
        theme_bw() +
        theme(legend.position = "right",
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_color_manual(values = colors)
    } else if (plot_method == "2") {
      p1 <- ggplot(df, aes(Condition, Value, colour = Condition)) +
        ylab("Relative frequency") +
        geom_boxplot() + geom_jitter(width = 0.2, alpha = 0.7) + theme_bw() +
        facet_wrap(~ Taxa, scales = "free", ncol = 3) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_color_manual(values = colors)
    } else if (plot_method == "3") {
      p1 <- ggplot(res_tax, aes(x = as.numeric(All_mean_relative_frequency), y = as.numeric(log2FC), color = Significant)) +
        geom_point(size = 2) +
        scale_x_log10() +
        scale_color_manual(values = c("black", "red")) +
        labs(x = "log10(Mean relative abundance)", y = "Log2FC") +
        theme_bw() +
        geom_text(
          data = subset(res_tax, Significant == "Yes"),
          aes(label = OTU), size = 3.5, vjust = 1
        )
    } else if (plot_method == "4") {
      # Heatmap for significant taxa
      OTU_input2 <- OTU_input1[rownames(OTU_input1) %in% rownames(res_tax_sig), , drop = FALSE]
      df1 <- scale(OTU_input2)
      # annotation aligned to columns of df1
      annot_vec <- setNames(group_index1$Condition, group_index1$Samples)[colnames(df1)]
      color_mapping <- setNames(colors, unique_conditions)
      ha <- ComplexHeatmap::HeatmapAnnotation(Group = annot_vec, col = list(Group = color_mapping))
      # simple continuous palette for matrix
      hm_cols <- colorRampPalette(c("navy", "white", "firebrick"))(101)
      p4 <- ComplexHeatmap::Heatmap(
        df1, name = "Z", show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = 7),
        show_column_names = TRUE, column_names_gp = grid::gpar(fontsize = 7),
        show_row_dend = TRUE, show_column_dend = FALSE, cluster_columns = FALSE,
        top_annotation = ha, col = hm_cols
      )
      p1 <- ggplotify::as.ggplot(p4)
    } else {
      message("Unknown plot_method; returning results without plot.")
    }
  }
  
  return(list(
    wilcoxtest_result_summary    = text_wilcoxtest_summary,
    wilcoxtest_result_significant= if (exists("res_tax_sig")) res_tax_sig else NULL,
    wilcoxtest_result_table      = res_tax,
    wilcoxtest_relative_frequency= rel_freq,
    total_counts                 = parent_seq_count,
    plot                         = p1
  ))
}
