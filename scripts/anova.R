library(ggplot2)
library(tibble)
library(qvalue)

anova_summary <- function(
    OTU_input, group_index, index_pvalue, ad_hoc, plot_method, alpha, anova_color_palette,
    show_head = TRUE, head_n = 5, show_colors = TRUE
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

  
  # Validate palette and prep condition colors
  if (!anova_color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    stop("Invalid `anova_color_palette`. Choose a valid RColorBrewer palette.")
  }
  unique_conditions <- unique(group_index$Condition)
  n_cond <- length(unique_conditions)
  max_n  <- RColorBrewer::brewer.pal.info[anova_color_palette, "maxcolors"]
  base_cols <- RColorBrewer::brewer.pal(min(max_n, max(3, n_cond)), anova_color_palette)
  cond_cols <- colorRampPalette(base_cols)(n_cond)
  if (show_head && show_colors) {
    msg <- paste0(
      "\n--- colors (first few) ---\n",
      paste(utils::capture.output(utils::head(colors, head_n)), collapse = "\n"),
      "\n"
    )
    cat(msg, file = stderr()); flush.console()
  }
  
  # Keep nonzero taxa
  OTU_input  <- OTU_input[rowSums(OTU_input[]) > 0, ]
  OTU_input1 <- OTU_input
  
  # Per-sample totals
  parent_seq_count <- as.data.frame(colSums(OTU_input))
  parent_seq_count <- tibble::rownames_to_column(parent_seq_count, "Samples")
  colnames(parent_seq_count)[2] <- "Total_counts"
  log_head(parent_seq_count, "parent_seq_count")
  
  # Relative frequency (%)
  rel_freq  <- sweep(OTU_input, 2, colSums(OTU_input), "/") * 100
  OTU_input <- as.data.frame(rel_freq)
  log_head(OTU_input, "OTU_input (% rel. freq.)")
  
  # Keep a full copy for heatmap annotation later
  group_index1 <- group_index
  
  # Use only Conditions vector for tests
  group_index <- group_index[, 2]
  groups <- names(table(group_index))
  
  if (length(groups) <= 2) {
    text_anova_summary <- "If you have ≤2 groups, use a two-group method (e.g., t-test/Wilcoxon)."
    return(list(anova_result_summary = text_anova_summary))
  }
  
  # ---------------- One-way ANOVA across all groups ----------------
  K <- length(groups)
  test_summary <- matrix(NA, nrow = nrow(OTU_input), ncol = 3 + 2*K)
  colnames(test_summary) <- c(
    "OTU",
    paste("Present_in_no_of", groups, sep = "_"),
    paste("Mean_relative_frequency", groups, sep = "_"),
    "PValue", "q_value"
  )
  
  for (i in seq_len(nrow(OTU_input))) {
    OTU_tmp <- OTU_input[i, ]
    # Simple guard; ANOVA itself handles equal variances across groups
    test_tmp <- aov(as.numeric(OTU_tmp) ~ as.factor(group_index))
    for (j in seq_len(K)) {
      gname <- groups[j]
      gvals <- as.numeric(OTU_tmp[which(group_index == gname)])
      test_summary[i, paste("Mean_relative_frequency", gname, sep = "_")] <- signif(mean(gvals, na.rm = TRUE), 4)
      test_summary[i, paste("Present_in_no_of", gname, sep = "_")]        <- sum(gvals != 0, na.rm = TRUE)
    }
    test_summary[i, "OTU"]    <- rownames(OTU_input)[i]
    test_summary[i, "PValue"] <- signif(summary(test_tmp)[[1]][["Pr(>F)"]][1], 4)
  }
  
  test_summary2 <- test_summary[!is.na(test_summary[, "PValue"]), , drop = FALSE]
  test_summary2[, "q_value"] <- signif(p.adjust(as.numeric(test_summary2[, "PValue"]), "BH"), 4)
  
  # ---------------- Optional Tukey HSD post-hoc ----------------
  if (identical(ad_hoc, "Yes")) {
    npairs <- choose(K, 2)
    pairwise_p_matrix <- matrix(NA, nrow = nrow(test_summary2), ncol = npairs)
    pair_names <- NULL
    
    for (i in seq_len(nrow(test_summary2))) {
      OTU_tmp <- OTU_input[match(test_summary2[i, "OTU"], rownames(OTU_input)), ]
      aov_fit  <- aov(as.numeric(OTU_tmp) ~ as.factor(group_index))
      thsd     <- TukeyHSD(aov_fit, conf.level = 0.95)
      padj_vec <- signif(thsd$`as.factor(group_index)`[, "p adj"], 4)
      if (is.null(pair_names)) pair_names <- rownames(thsd$`as.factor(group_index)`)
      pairwise_p_matrix[i, ] <- padj_vec
    }
    colnames(pairwise_p_matrix) <- pair_names
    test_summary2 <- cbind(test_summary2, pairwise_p_matrix)
  }
  
  res_tax <- as.data.frame(test_summary2)
  rownames(res_tax) <- res_tax[, 1]
  log_head(res_tax, "res_tax (ANOVA results)")
  
  # ---------------- Significance filter ----------------
  if (identical(index_pvalue, "padj")) {
    res_tax_sig <- subset(res_tax, q_value < alpha)
  } else if (identical(index_pvalue, "pvalue")) {
    res_tax_sig <- subset(res_tax, PValue < alpha)
  } else {
    stop("`index_pvalue` must be 'padj' or 'pvalue'.")
  }
  log_head(res_tax_sig, "res_tax_sig (significant)")
  
  res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig), "Yes", "No")
  
  # ---------------- Plot data (if any significant) ----------------
  df  <- NULL
  df1 <- NULL
  if (nrow(res_tax_sig) >= 1) {
    # z-scored raw counts (not % rel freq) for heatmap
    OTU_input2 <- OTU_input1[rownames(OTU_input1) %in% rownames(res_tax_sig), , drop = FALSE]
    df1 <- scale(OTU_input2)
    
    # long format for box/violin
    data_long <- as.data.frame(t(OTU_input))  # samples x taxa (in %)
    label_col <- if (identical(index_pvalue, "padj")) "q_value" else "PValue"
    df <- do.call(rbind, lapply(res_tax[rownames(res_tax_sig), "OTU"], function(i) {
      data.frame(
        Value     = data_long[, i],
        Condition = group_index,
        Taxa1     = i,
        Taxa      = paste(i,
                          if (label_col == "q_value") " q_value = " else " PValue = ",
                          round(as.numeric(res_tax[i, label_col]), 4), sep = ""),
        check.names = FALSE
      )
    }))
    log_head(df, "df (plotting table)")
  }
  
  text_anova_summary <- if (nrow(res_tax_sig) >= 1) {
    paste0("Total of ", nrow(res_tax_sig), " taxa were identified as significant.")
  } else if (identical(index_pvalue, "padj")) {
    "No significant taxa were identified. Adjust the FDR value or change the test correction method."
  } else {
    "No significant taxa were identified. Adjust the P-value value or change the test correction method."
  }
  
  # ---------------- Plots ----------------
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
        scale_color_manual(values = cond_cols)
      
    } else if (plot_method == "2") {
      p1 <- ggplot(df, aes(Condition, Value, colour = Condition)) +
        ylab("Relative frequency") +
        geom_boxplot() + geom_jitter(width = 0.2, alpha = 0.7) + theme_bw() +
        facet_wrap(~ Taxa, scales = "free", ncol = 3) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_color_manual(values = cond_cols)
      
    } else if (plot_method == "3") {
      # Heatmap: align annotation to column order and use named color mapping
      gi <- if ("Samples" %in% colnames(group_index1)) {
        group_index1 %>% tibble::remove_rownames() %>% tibble::column_to_rownames("Samples")
      } else {
        group_index1
      }
      annot_vec <- setNames(gi$Condition, rownames(gi))[colnames(df1)]
      conds     <- unique(annot_vec)
      ann_cols  <- setNames(colorRampPalette(base_cols)(length(conds)), conds)
      
      ha <- ComplexHeatmap::HeatmapAnnotation(
        Group = annot_vec,
        col   = list(Group = ann_cols)
      )
      # continuous colormap for z-scored values
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
    anova_result_summary        = text_anova_summary,
    anova_result_significant    = if (exists("res_tax_sig")) res_tax_sig else NULL,
    anova_result_table          = res_tax,
    anova_relative_frequency    = rel_freq,
    total_counts                = parent_seq_count,
    plot                        = p1
  ))
}
