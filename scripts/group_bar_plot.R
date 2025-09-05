#library("RColorBrewer")

Bar_Plot_Group_OTU <- function(OTU_input, group_index, color_palette_group, n_top = 15, names = "OTU", plot_type, show_head = TRUE, head_n = 5){
  
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
  
  # remove empty rows
  OTU_input <- OTU_input[rowSums(OTU_input[]) > 0, ]
  log_head(OTU_input, "OTU_input (after removing empty rows)")
  
  # convert to relative freq if requested
  if (plot_type == "3"){
    rel_freq <- (OTU_input) / (colSums(OTU_input)) * 100
    OTU_input <- as.data.frame(rel_freq)
    log_head(OTU_input, "OTU_input (relative frequency)")
  }
  
  if (n_top > nrow(OTU_input)){
    print("Number of selected top OTU number is greater than the total OTU numbers")
  } else {
    if (ncol(OTU_input) != length(group_index)){
      print("The length of the group index does not match the number of columns in OTU table")
    } else {
      Overall_mean <- apply(OTU_input, 1, mean, na.rm = TRUE)
      Top_index <- order(Overall_mean, decreasing = TRUE)[1:n_top]
      Top_OTU_names <- c(rownames(OTU_input)[Top_index], "OTHERS")
      group_tmp <- names(table(group_index))
      group <- c()
      OTU <- c()
      Abundance <- c()
      
      for (i in 1:length(group_tmp)){
        matrix_tmp <- OTU_input[, which(group_index == group_tmp[i])]
        top_means <- apply(matrix_tmp[Top_index, ], 1, mean, na.rm = TRUE)
        other_mean <- mean(unlist(matrix_tmp[setdiff(c(1:nrow(OTU_input)), Top_index), ]))
        group <- c(group, rep(group_tmp[i], n_top + 1))
        OTU <- c(OTU, Top_OTU_names)
        Abundance <- c(Abundance, c(top_means, other_mean))
      }
      
      ## data structure for ggplot ##
      data_gg <- data.frame(group, OTU, Abundance)
      level_OTU <- Top_OTU_names[length(Top_OTU_names):1]
      
      log_head(data_gg, "data_gg (final ggplot input)")
      log_head(level_OTU, "level_OTU (OTU order)")
    }
  }
  
  colors <- colorRampPalette(brewer.pal(min(12, length(level_OTU)), color_palette_group))(length(level_OTU))
  
  if (plot_type == "1"){
    ggplot(data_gg, aes(fill = factor(OTU, level = level_OTU), y = Abundance, x = group)) + 
      geom_bar(stat = "identity", position = "fill") + 
      labs(y = "Abundance(%)", fill = names) +
      scale_fill_manual(values = colors) +
      scale_y_continuous(labels = scales::percent_format()) +
      theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20),
            legend.title = element_text(size = 15), legend.text = element_text(size = 10)) +
      theme_bw()
  }
  else if (plot_type == "2"){
    ggplot(data_gg, aes(fill = factor(OTU, level = level_OTU), y = Abundance, x = group)) + 
      geom_bar(stat = "identity", position = "stack") + 
      labs(y = "Abundance value", fill = names) +
      scale_fill_manual(values = colors) +
      theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20),
            legend.title = element_text(size = 15), legend.text = element_text(size = 10)) +
      theme_bw()
  }
  else if (plot_type == "3"){
    ggplot(data_gg, aes(fill = factor(OTU, level = level_OTU), y = Abundance, x = group)) + 
      geom_bar(stat = "identity", position = "fill") + 
      labs(y = "Relative frequency", fill = names) +
      scale_y_continuous(labels = scales::percent_format()) +
      scale_fill_manual(values = colors) +
      theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20),
            legend.title = element_text(size = 15), legend.text = element_text(size = 10)) +
      theme_bw()
  }
  else {
    print("select plot type")
  }
}
