differential_individual_boxplots <- function(df, colors, y_label, facets_per_page = 25, ncol = 5, nrow = 5) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df) || !"Taxa" %in% names(df)) {
    return(list())
  }

  taxa_levels <- unique(as.character(df$Taxa))
  taxa_levels <- taxa_levels[!is.na(taxa_levels) & nzchar(taxa_levels)]

  if (!length(taxa_levels)) {
    return(list())
  }

  condition_levels <- unique(as.character(df$Condition))
  condition_levels <- condition_levels[!is.na(condition_levels) & nzchar(condition_levels)]

  if (!is.null(colors) && length(colors) >= length(condition_levels) && is.null(names(colors))) {
    colors <- setNames(colors[seq_along(condition_levels)], condition_levels)
  }

  taxa_pages <- split(taxa_levels, ceiling(seq_along(taxa_levels) / facets_per_page))
  plots <- lapply(seq_along(taxa_pages), function(page_index) {
    taxa_page <- taxa_pages[[page_index]]
    plot_df <- df[df$Taxa %in% taxa_page, , drop = FALSE]
    plot_df$Taxa <- factor(plot_df$Taxa, levels = taxa_page)

    ggplot(plot_df, aes(Condition, Value, colour = Condition)) +
      ylab(y_label) +
      geom_boxplot() +
      geom_jitter(width = 0.2, alpha = 0.7) +
      theme_bw() +
      facet_wrap(~ Taxa, scales = "free", ncol = ncol, nrow = nrow) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_color_manual(values = colors)
  })

  names(plots) <- vapply(seq_along(taxa_pages), function(page_index) {
    start_index <- ((page_index - 1) * facets_per_page) + 1
    end_index <- min(page_index * facets_per_page, length(taxa_levels))
    paste0("Plot ", page_index, " (taxa ", start_index, "-", end_index, ")")
  }, character(1))

  plots
}

differential_plot_top_n <- function(plot_top_n, default = 25, minimum = 2, maximum = 100) {
  plot_top_n <- suppressWarnings(as.integer(plot_top_n))

  if (!length(plot_top_n) || is.na(plot_top_n[1])) {
    plot_top_n <- default
  }

  min(maximum, max(minimum, plot_top_n[1]))
}

differential_plot_taxa <- function(res_tax_sig, plot_top_n, plot_taxa_mode = "all") {
  if (is.null(res_tax_sig) || !is.data.frame(res_tax_sig) || !nrow(res_tax_sig)) {
    return(character(0))
  }

  taxa <- rownames(res_tax_sig)

  if (is.null(taxa) || any(!nzchar(taxa))) {
    taxa <- as.character(res_tax_sig[[1]])
  }

  taxa <- taxa[!is.na(taxa) & nzchar(taxa)]
  plot_taxa_mode <- tolower(as.character(plot_taxa_mode))
  if (!length(plot_taxa_mode) || is.na(plot_taxa_mode[1]) || !nzchar(plot_taxa_mode[1])) {
    plot_taxa_mode <- "all"
  } else {
    plot_taxa_mode <- plot_taxa_mode[1]
  }

  if (plot_taxa_mode %in% c("selected", "selected_number", "top_n", "limited")) {
    taxa[seq_len(min(length(taxa), differential_plot_top_n(plot_top_n)))]
  } else {
    taxa
  }
}
