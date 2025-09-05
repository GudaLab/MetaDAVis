library("ggplot2")
library("tibble")
library("qvalue")
library("dunn.test")
library("RColorBrewer")
library("ComplexHeatmap")
library("ggplotify")
library("circlize")

kruskal_wallis_test_summary <- function(
    OTU_input, group_index, index_pvalue, ad_hoc, plot_method, alpha,
    kruskal_wallis_test_color_palette,
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

  
  # Palette validation + condition colors
  if (!kruskal_wallis_test_color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    stop("Invalid `kruskal_wallis_test_color_palette`. Choose a valid RColorBrewer palette.")
  }
  unique_conditions <- unique(group_index$Condition)
  n_cond <- length(unique_conditions)
  max_n  <- RColorBrewer::brewer.pal.info[kruskal_wallis_test_color_palette, "maxcolors"]
  base_cols <- RColorBrewer::brewer.pal(min(max_n, max(3, n_cond)), kruskal_wallis_test_color_palette)
  cond_cols <- colorRampPalette(base_cols)(n_cond)
  if (show_head && show_colors) {
    msg <- paste0(
      "\n--- colors (first few) ---\n",
      paste(utils::capture.output(utils::head(colors, head_n)), collapse = "\n"),
      "\n"
    )
    cat(msg, file = stderr()); flush.console()
  }
  
  OTU_input <- OTU_input[rowSums(OTU_input[]) > 0, ]
  OTU_input1 <- OTU_input
  
  # totals
  parent_seq_count <- as.data.frame(colSums(OTU_input))
  parent_seq_count <- tibble::rownames_to_column(parent_seq_count, "Samples")
  colnames(parent_seq_count)[2] <- "Total_counts"
  
  # relative frequency (%)
  rel_freq <- sweep(OTU_input, 2, colSums(OTU_input), "/") * 100
  OTU_input <- as.data.frame(rel_freq)
  
  group_index1 <- group_index                       # keep full metadata for later (heatmap)
  group_index   <- group_index[, 2]                 # just the Condition vector used for tests
  
  if (length(names(table(group_index))) <= 2) {
    text_kruskal_wallis_test_summary <- "If you have ≤2 groups, use a two-group method (e.g., Wilcoxon/t-test)."
    return(list(kruskal_wallis_test_result_summary = text_kruskal_wallis_test_summary))
  }
  
  # ---------------- Kruskal–Wallis across all groups ----------------
  K <- length(names(table(group_index)))
  test_summary <- matrix(NA, nrow = nrow(OTU_input), ncol = 3 + 2*K)
  colnames(test_summary) <-
    c("OTU",
      paste("Present_in_no_of", names(table(group_index)), sep = "_"),
      paste("Mean_relative_frequency", names(table(group_index)), sep = "_"),
      "PValue", "q_value")
  
  for (i in seq_len(nrow(OTU_input))) {
    OTU_tmp <- OTU_input[i, ]
    # guard against all-NA/constant-by-group pathologies
    if (all(sapply(seq_len(K), function(k)
      sd(OTU_tmp[which(group_index == names(table(group_index))[k])], na.rm = TRUE)) >= 0)) {
      
      test_tmp <- kruskal.test(as.numeric(OTU_tmp) ~ as.factor(group_index))
      
      for (j in seq_len(K)) {
        gname <- names(table(group_index))[j]
        gvals <- as.numeric(OTU_tmp[which(group_index == gname)])
        test_summary[i, paste("Mean_relative_frequency", gname, sep = "_")] <- signif(mean(gvals, na.rm = TRUE), 4)
        test_summary[i, paste("Present_in_no_of", gname, sep = "_")] <- sum(gvals != 0, na.rm = TRUE)
      }
      test_summary[i, "OTU"]    <- rownames(OTU_input)[i]
      test_summary[i, "PValue"] <- signif(test_tmp$p.value, 4)
    }
  }
  
  test_summary2 <- test_summary[!is.na(test_summary[, "PValue"]), , drop = FALSE]
  test_summary2[, "q_value"] <- signif(p.adjust(as.numeric(test_summary2[, "PValue"]), "BH"), 4)
  
  # ---------------- Optional Dunn’s post-hoc ----------------
  if (identical(ad_hoc, "Yes")) {
    # number of pairwise comps = K choose 2
    npairs <- K * (K - 1) / 2
    pairwise_p_matrix <- matrix(NA, nrow = nrow(test_summary2), ncol = npairs)
    
    pair_names <- NULL
    for (i in seq_len(nrow(test_summary2))) {
      OTU_tmp <- OTU_input[match(test_summary2[i, "OTU"], rownames(OTU_input)), ]
      dt <- dunn.test::dunn.test(as.numeric(OTU_tmp), group_index, method = "bh", kw = FALSE, list = TRUE)
      pvec <- signif(dt$P, 4)
      if (is.null(pair_names)) pair_names <- dt$comparisons
      pairwise_p_matrix[i, ] <- pvec
    }
    colnames(pairwise_p_matrix) <- pair_names
    test_summary2 <- cbind(test_summary2, pairwise_p_matrix)
  }
  
  res_tax <- as.data.frame(test_summary2)
  rownames(res_tax) <- res_tax[, 1]
  
  # ---------------- Significance filter ----------------
  if (identical(index_pvalue, "padj")) {
    res_tax_sig <- subset(res_tax, q_value < alpha)
  } else if (identical(index_pvalue, "pvalue")) {
    res_tax_sig <- subset(res_tax, PValue < alpha)
  } else {
    stop("`index_pvalue` must be 'padj' or 'pvalue'.")
  }
  
  # build plotting tables if any significant
  df <- NULL
  df1 <- NULL
  if (nrow(res_tax_sig) >= 1) {
    res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig), "Yes", "No")
    
    # z-scored matrix for heatmap
    OTU_input2 <- OTU_input1[rownames(OTU_input1) %in% rownames(res_tax_sig), , drop = FALSE]
    df1 <- scale(OTU_input2)
    
    # long frame for box/violin
    data_long <- as.data.frame(t(OTU_input))
    df <- do.call(rbind, lapply(res_tax[rownames(res_tax_sig), "OTU"], function(i) {
      data.frame(
        Value = data_long[, i],
        Condition = group_index,
        Taxa1 = i,
        Taxa = paste0(i, if (identical(index_pvalue, "padj")) " q_value = " else " PValue = ",
                      round(as.numeric(res_tax[i, if (identical(index_pvalue, "padj")) "q_value" else "PValue"]), 4)),
        check.names = FALSE
      )
    }))
  }
  
  text_kruskal_wallis_test_summary <- if (nrow(res_tax_sig) >= 1) {
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
      # Heatmap annotation: align to df1 columns and use NAMED vector for colors
      # Ensure rownames of metadata are sample IDs
      group_index1_rown <- group_index1 %>%
        tibble::remove_rownames() %>%
        tibble::column_to_rownames("Samples")
      annot_vec <- setNames(group_index1_rown$Condition, rownames(group_index1_rown))[colnames(df1)]
      conds     <- unique(annot_vec)
      ann_cols  <- setNames(colorRampPalette(base_cols)(length(conds)), conds)
      
      ha <- ComplexHeatmap::HeatmapAnnotation(
        Group = annot_vec,
        col   = list(Group = ann_cols)
      )
      # continuous colormap for matrix values
      hm_cols <- circlize::colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick"))
      
      p3 <- ComplexHeatmap::Heatmap(
        df1, name = "Z",
        show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = 7),
        show_column_names = TRUE, column_names_gp = grid::gpar(fontsize = 7),
        show_row_dend = TRUE, show_column_dend = FALSE, cluster_columns = FALSE,
        top_annotation = ha, col = hm_cols
      )
      p1 <- ggplotify::as.ggplot(p3)
    }
  }
  
  return(list(
    kruskal_wallis_test_result_summary       = text_kruskal_wallis_test_summary,
    kruskal_wallis_test_result_significant   = if (exists("res_tax_sig")) res_tax_sig else NULL,
    kruskal_wallis_test_result_table         = res_tax,
    kruskal_wallis_test_relative_frequency   = rel_freq,
    total_counts                             = parent_seq_count,
    plot                                     = p1
  ))
}
