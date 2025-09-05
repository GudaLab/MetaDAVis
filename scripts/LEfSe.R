library(phyloseq)
library(lefser)
library(mia)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

LEfSe_summary <- function(
    OTU_input,
    group_index,
    group1,
    group2,
    select_LEfSe_method      = "BH",   # p-adjust: "BH","holm","bonferroni","BY","none"
    select_LEfSe_pvalue      = 0.05,   # Kruskal–Wallis threshold
    select_LEfSe_threshold   = 2,      # LDA score threshold
    select_LEfSe_color_palette = "Set2",
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
  
  # --- early checks ---
  if (group1 == group2) {
    txt <- "Please select different group in condition1 or condition2"
    message(txt)
    return(list(LEfSe_result_summary = txt))
  }
  if (!select_LEfSe_color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    stop("Invalid `select_LEfSe_color_palette`. Choose a valid RColorBrewer palette.")
  }
  
  # --- keep only the two groups & align sample order ---
  group_index <- subset(group_index, (group_index$Condition %in% c(group1, group2)))
  log_head(group_index, "group_index (filtered to two groups)")
  
  # transpose -> samples as rows
  OTU_input <- as.data.frame(t(OTU_input))
  OTU_input <- tibble::rownames_to_column(OTU_input, "Samples")
  # keep & order to metadata
  OTU_input <- OTU_input[OTU_input$Samples %in% group_index$Samples, ]
  OTU_input <- OTU_input[match(group_index$Samples, OTU_input$Samples), , drop = FALSE]
  # back to matrix with samples as rownames
  OTU_input <- OTU_input %>% tibble::remove_rownames() %>% tibble::column_to_rownames("Samples")
  # back to taxa x samples
  OTU_input <- as.data.frame(t(OTU_input))
  # drop empty taxa
  OTU_input <- OTU_input[rowSums(OTU_input[]) > 0, , drop = FALSE]
  log_head(OTU_input, "OTU_input (aligned & nonzero taxa)")
  
  # totals
  parent_seq_count <- as.data.frame(colSums(OTU_input))
  parent_seq_count <- tibble::rownames_to_column(parent_seq_count, "Samples")
  colnames(parent_seq_count)[2] <- "Total_counts"
  log_head(parent_seq_count, "parent_seq_count")
  
  # --- build phyloseq ---
  group_index1 <- group_index %>% tibble::remove_rownames() %>% tibble::column_to_rownames("Samples")
  tax_index <- data.frame(Species = rownames(OTU_input), row.names = rownames(OTU_input), check.names = FALSE)
  physeq1 <- phyloseq(
    otu_table(OTU_input, taxa_are_rows = TRUE),
    tax_table(as.matrix(tax_index)),
    sample_data(group_index1)
  )
  
  # --- colors (2 conditions) ---
  conds <- unique(sample_data(physeq1)$Condition)
  max_n <- RColorBrewer::brewer.pal.info[select_LEfSe_color_palette, "maxcolors"]
  base_cols <- RColorBrewer::brewer.pal(min(max_n, max(3, length(conds))), select_LEfSe_color_palette)
  colors <- colorRampPalette(base_cols)(length(conds))
  if (show_head && show_colors) {
    msg <- paste0(
      "\n--- colors (first few) ---\n",
      paste(utils::capture.output(utils::head(colors, head_n)), collapse = "\n"),
      "\n"
    )
    cat(msg, file = stderr()); flush.console()
  }
  log_head_vec(conds, "Condition levels")
  
  # --- convert and run LEfSe ---
  lefse_input <- mia::convertFromPhyloseq(physeq1)
  
  # sanitize method
  valid_methods <- c("BH","holm","bonferroni","BY","none")
  if (!select_LEfSe_method %in% valid_methods) {
    warning(sprintf("Unknown method '%s'; defaulting to 'BH'.", select_LEfSe_method))
    select_LEfSe_method <- "BH"
  }
  
  res <- lefser::lefser(
    lefse_input,
    classCol          = "Condition",
    kruskal.threshold = as.numeric(select_LEfSe_pvalue),
    lda.threshold     = as.numeric(select_LEfSe_threshold),
    method            = select_LEfSe_method
  )
  
  # --- plot & summarize ---
  plots1 <- NULL
  data_tbl <- NULL
  txt <- NULL
  
  # Try to plot; handle case with no significant taxa
  plots1 <- tryCatch(
    {
      p <- lefser::lefserPlot(res, colors = colors)
      # round numeric columns in the attached data
      if (!is.null(p$data)) {
        data_tbl <- p$data %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(., 3)))
        log_head(data_tbl, "lefserPlot$data (rounded)")
        n_sig <- nrow(data_tbl)
        txt <- paste("Total of ", n_sig, " taxa were identified as significant.", sep = "")
      } else {
        txt <- "LEfSe completed, but no plotting data were produced (no significant taxa?)."
      }
      p
    },
    error = function(e){
      txt <<- paste0("LEfSe completed but plotting failed: ", e$message)
      NULL
    }
  )
  
  if (is.null(txt)) {
    txt <- "LEfSe completed."
  }
  
  return(list(
    LEfSe_result_summary = txt,
    data1 = data_tbl,
    data2 = parent_seq_count,
    plot  = plots1
  ))
}
