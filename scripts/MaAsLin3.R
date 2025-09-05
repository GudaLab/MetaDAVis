library("maaslin3")
library("tibble")
library("dplyr")

MaAsLin3_summary <- function(
    OTU_input,
    group_index,
    group1,
    group2,
    select_MaAsLin3_normalization,
    select_MaAsLin3_transformation,
    select_MaAsLin3_correction,
    select_MaAsLin3_pvalue,
    output_dir = file.path(getwd(), "/www/hmp2_output"),
    show_head = TRUE,
    head_n = 5,
    show_files = TRUE
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

  
  if (group1 == group2) {
    txt <- "Please select different group in condition1 or condition2"
    message(txt)
    return(list(MaAsLin3_result_summary = txt))
  }
  
  # 1) keep only the two groups
  group_index <- subset(group_index, group_index$Condition %in% c(group1, group2))
  log_head(group_index, "group_index (filtered to two groups)")
  
  # 2) transpose -> samples as rows, align to metadata order, back to taxa x samples
  OTU_input <- as.data.frame(t(OTU_input))
  OTU_input <- tibble::rownames_to_column(OTU_input, "Samples")
  OTU_input <- OTU_input[OTU_input$Samples %in% group_index$Samples, ]
  OTU_input <- OTU_input[match(group_index$Samples, OTU_input$Samples), , drop = FALSE]
  OTU_input <- OTU_input %>% tibble::remove_rownames() %>% tibble::column_to_rownames("Samples")
  OTU_input <- as.data.frame(t(OTU_input))
  OTU_input <- OTU_input[rowSums(OTU_input[]) > 0, , drop = FALSE]
  log_head(OTU_input, "OTU_input (aligned & nonzero taxa)")
  
  # 3) totals
  parent_seq_count <- as.data.frame(colSums(OTU_input))
  parent_seq_count <- tibble::rownames_to_column(parent_seq_count, "Samples")
  colnames(parent_seq_count)[2] <- "Total_counts"
  log_head(parent_seq_count, "parent_seq_count (per-sample totals)")
  
  # 4) metadata with rownames = Samples, ensure same order as columns of OTU_input
  group_index1 <- group_index %>% tibble::remove_rownames() %>% tibble::column_to_rownames("Samples")
  group_index1 <- group_index1[colnames(OTU_input), , drop = FALSE]
  log_head(group_index1, "group_index1 (aligned to OTU_input cols)")
  
  # 5) ensure output directory exists (clean old run optionally)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 6) run MaAsLin3
  set.seed(1)
  fit_out <- maaslin3(
    input_data                    = OTU_input,
    input_metadata                = group_index1,
    output                        = output_dir,
    formula                       = 'Condition',
    normalization                 = select_MaAsLin3_normalization,   # e.g., "TSS","RA","CLR","NONE"
    transform                     = select_MaAsLin3_transformation,  # e.g., "LOG","NONE","AST","CSS"
    correction                    = select_MaAsLin3_correction,      # e.g., "BH","holm","bonferroni","BY","none"
    augment                       = TRUE,
    standardize                   = TRUE,
    plot_summary_plot             = TRUE,
    plot_associations             = TRUE,
    max_significance              = select_MaAsLin3_pvalue,
    median_comparison_abundance   = TRUE,
    median_comparison_prevalence  = FALSE,
    summary_plot_first_n          = 25,
    max_pngs                      = 30,
    cores                         = 1
  )
  
  # 7) zip outputs
  if (!requireNamespace("zip", quietly = TRUE)) {
    warning("Package 'zip' is not installed; skipping ZIP creation.")
    output_zip <- NULL
  } else {
    output_zip <- file.path(dirname(output_dir), paste0(basename(output_dir), ".zip"))
    if (file.exists(output_zip)) unlink(output_zip)
    zip::zipr(zipfile = output_zip, files = list.files(output_dir, full.names = TRUE, recursive = TRUE))
  }
  
  if (isTRUE(show_head) && isTRUE(show_files)) {
    if (dir.exists(output_dir)) {
      files <- list.files(output_dir, recursive = TRUE, include.dirs = FALSE)
      n <- min(length(files), head_n)
      
      cat(sprintf(
        "\n--- output files (first %d of %d) in `%s` ---\n",
        n, length(files), normalizePath(output_dir, mustWork = FALSE)
      ), file = stderr())
      
      if (n > 0) {
        cat(paste(utils::head(files, n), collapse = "\n"), "\n", file = stderr())
      } else {
        cat("(no files found)\n", file = stderr())
      }
      flush.console()
    } else {
      cat(sprintf("\n--- output files ---\nDirectory does not exist: %s\n",
                  output_dir), file = stderr())
      flush.console()
    }
  }
  
  text_maaslin3_summary <- if (!is.null(output_zip)) {
    "Done. Use the 'Download Output as ZIP' button to get the results."
  } else {
    "Done. Results were written to the output directory."
  }
  
  return(list(
    MaAsLin3_result_summary = text_maaslin3_summary,
    output_dir              = output_dir,
    output_zip              = if (!is.null(output_zip)) output_zip else NA_character_,
    total_counts            = parent_seq_count,
    fit                     = fit_out
  ))
}
