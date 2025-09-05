#library("RColorBrewer")   # needed for brewer.pal()
#library("ggplot2")        # needed for ggplot()

Bar_Plot_Individual_OTU <- function(
    OTU_input,
    group_index,
    color_palette_individual,
    n_top = 15,
    names = "OTU",
    plot_type,
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
  
  # 1) remove empty rows
  OTU_input <- OTU_input[rowSums(OTU_input[]) > 0, ]
  log_head(OTU_input, "OTU_input (after removing empty rows)")
  
  # 2) relative frequency if requested
  if (plot_type == "3"){
    rel_freq <- (OTU_input) / (colSums(OTU_input)) * 100
    OTU_input <- as.data.frame(rel_freq)
    log_head(OTU_input, "OTU_input (relative frequency)")
  }
  
  # 3) set up grouping
  group_index1 <- group_index[, 2]
  if (n_top > nrow(OTU_input)){
    print("Number of selected top OTU number is greater than the total OTU numbers")
  } else {
    if (ncol(OTU_input) != length(group_index1)){
      print("The length of the group index does not match the number of columns in OTU table")
    } else {
      Overall_mean <- apply(OTU_input, 1, mean, na.rm = TRUE)
      Top_index <- order(Overall_mean, decreasing = TRUE)[1:n_top]
      Top_OTU_names <- c(rownames(OTU_input)[Top_index], "OTHERS")
      log_head_vec(Top_OTU_names, "Top_OTU_names (selected + OTHERS)")
      
      group_tmp <- names(table(group_index1))
      Samples <- OTU <- Abundance <- c()
      
      for (i in 1:length(group_tmp)){
        matrix_tmp <- OTU_input[, which(group_index1 == group_tmp[i])]
        for (j in 1:ncol(matrix_tmp)){
          top_OTU_abd_tmp <- matrix_tmp[Top_index, j]
          others_OTU_abd_tmp <- 1 - sum(top_OTU_abd_tmp)
          if (others_OTU_abd_tmp < 0) others_OTU_abd_tmp <- 0
          Samples <- c(Samples, rep(colnames(matrix_tmp)[j], n_top + 1))
          OTU     <- c(OTU, Top_OTU_names)
          Abundance <- c(Abundance, c(top_OTU_abd_tmp, others_OTU_abd_tmp))
        }
      }
      
      # 4) ggplot data
      data_gg <- data.frame(Samples, OTU, Abundance)
      level_order <- unique(as.character(data_gg$Samples))
      level_OTU <- Top_OTU_names[length(Top_OTU_names):1]
      log_head(data_gg, "data_gg (ggplot input)")
      log_head_vec(level_OTU, "level_OTU (OTU order)")
      
      # 5) merge with metadata
      result <- merge(data_gg, group_index, by = "Samples")
      log_head(result, "result (merged with group_index)")
      
      # 6) colors
      colors <- colorRampPalette(brewer.pal(min(12, length(level_OTU)), color_palette_individual))(length(level_OTU))
      if (show_head && show_colors) {
        msg <- paste0(
          "\n--- colors (first few) ---\n",
          paste(utils::capture.output(utils::head(colors, head_n)), collapse = "\n"),
          "\n"
        )
        cat(msg, file = stderr()); flush.console()
      }
      
      # 7) plot
      if (plot_type == "1"){
        ggplot(result, aes(fill = factor(OTU, level = level_OTU), y = Abundance, x = Samples)) +
          geom_bar(stat = "identity", position = "fill") +
          facet_wrap(. ~ Condition, scales = "free") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
          labs(y = "Abundance(%)", fill = names) +
          scale_fill_manual(values = colors) +
          scale_y_continuous(labels = scales::percent_format())
      }
      else if (plot_type == "2"){
        ggplot(result, aes(fill = factor(OTU, level = level_OTU), y = Abundance, x = Samples)) +
          geom_bar(stat = "identity", position = "stack") +
          facet_wrap(. ~ Condition, scales = "free") +
          scale_fill_manual(values = colors) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
          labs(y = "Abundance value", fill = names)
      }
      else if (plot_type == "3"){
        ggplot(result, aes(fill = factor(OTU, level = level_OTU), y = Abundance, x = Samples)) +
          geom_bar(stat = "identity", position = "fill") +
          facet_wrap(. ~ Condition, scales = "free") +
          scale_fill_manual(values = colors) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
          labs(y = "Relative frequency", fill = names) +
          scale_y_continuous(labels = scales::percent_format())
      }
      else {
        print("select plot type")
      }
    }
  }
}
