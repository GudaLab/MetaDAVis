library("phyloseq")
library("vegan")
library("reshape2")
library("ggplot2")
library("RColorBrewer")
library("ggpubr")

beta_diversity_boxplot <- function(
    OTU_input,
    group_index,
    tax_index,
    select_beta_color_palette,
    index_method,
    plot_method,
    adonis_permutations,
    adonis_dissimilarities,   # renamed for clarity (was adonis_dissimilarities)
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
  
  ## Taxonomy prep
  colnames(tax_index)[1] <- "Species"
  rownames(tax_index) <- tax_index[, 1]
  log_head(tax_index, "tax_index (with Species col & rownames)")
  taxmat <- as.matrix(tax_index)
  
  ## Phyloseq object
  OTU <- otu_table(OTU_input, taxa_are_rows = TRUE)
  TAX <- tax_table(taxmat)
  sampledata <- sample_data(group_index)
  physeq1 <- phyloseq(OTU, TAX, sampledata)
  
  log_head(OTU_input,   "OTU_input (raw)")
  log_head(group_index, "group_index (sample metadata)")
  
  ## Abundance matrix (taxa as rows)
  abundance_data <- as.data.frame(otu_table(physeq1))
  if (!taxa_are_rows(physeq1)) {
    abundance_data <- t(abundance_data)
  }
  log_head(abundance_data, "abundance_data (before filtering zeros)")
  
  ## Remove empty rows
  empty_rows <- rowSums(abundance_data) == 0
  if (any(empty_rows)) {
    abundance_data <- abundance_data[!empty_rows, , drop = FALSE]
  }
  log_head(abundance_data, "abundance_data (after removing zero-sum taxa)")
  
  ## Sync samples between metadata and abundance
  shared_samples <- intersect(rownames(group_index), colnames(abundance_data))
  if (length(shared_samples) == 0) {
    stop("Error: No overlapping samples between abundance_data and group_index.")
  }
  group_index <- group_index[shared_samples, , drop = FALSE]
  abundance_data <- abundance_data[, shared_samples, drop = FALSE]
  log_head_vec(shared_samples, "shared_samples")
  
  ## Distance (e.g., 'bray', 'jaccard', etc. passed via index_method)
  bray_dist <- vegdist(t(abundance_data), method = index_method)
  log_head(as.matrix(bray_dist), paste0(index_method, " distance matrix (preview)"))
  
  ## PERMANOVA (adonis2)
  adonis_result <- adonis2(bray_dist ~ Condition,
                           data = group_index,
                           permutations = adonis_permutations,
                           sqrt.dist = adonis_dissimilarities)
  table_adonis <- as.data.frame(adonis_result)
  log_head(table_adonis, "adonis2 result")
  
  ## Ordination
  ordination_df <- NULL
  if (plot_method == "PCoA") {
    pcoa_results <- cmdscale(bray_dist, eig = TRUE, k = 2)
    ordination_df <- as.data.frame(pcoa_results$points)
    colnames(ordination_df) <- c("Axis1", "Axis2")
    ordination_title <- "PCoA Ordination"
  } else if (plot_method == "NMDS") {
    nmds_results <- metaMDS(bray_dist, k = 2)
    ordination_df <- as.data.frame(nmds_results$points)
    colnames(ordination_df) <- c("Axis1", "Axis2")
    ordination_title <- "NMDS Ordination"
  } else if (plot_method == "DCA") {
    dca_results <- decorana(t(abundance_data))
    ordination_df <- as.data.frame(scores(dca_results, display = "sites"))
    colnames(ordination_df) <- c("Axis1", "Axis2")
    ordination_title <- "DCA Ordination"
  } else if (plot_method == "RDA") {
    rda_results <- rda(t(abundance_data) ~ Condition, data = group_index)
    ordination_df <- as.data.frame(scores(rda_results, display = "sites"))
    colnames(ordination_df) <- c("Axis1", "Axis2")
    ordination_title <- "RDA Ordination"
  } else if (plot_method == "MDS") {
    mds_results <- MASS::isoMDS(as.matrix(bray_dist))
    ordination_df <- as.data.frame(mds_results$points)
    colnames(ordination_df) <- c("Axis1", "Axis2")
    ordination_title <- "MDS Ordination"
  } else {
    stop("Invalid ordination method selected. Choose from PCoA, NMDS, DCA, RDA, or MDS.")
  }
  
  ordination_df$Condition <- group_index$Condition
  log_head(ordination_df, paste0(ordination_title, " (scores preview)"))
  
  ## Validate and build colors
  if (!select_beta_color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    stop("Invalid color palette. Please choose a valid palette from RColorBrewer.")
  }
  conds <- unique(sample_data(physeq1)$Condition)
  n_cond <- length(conds)
  max_n <- RColorBrewer::brewer.pal.info[select_beta_color_palette, "maxcolors"]
  base_cols <- RColorBrewer::brewer.pal(min(max_n, max(3, min(max_n, n_cond))), select_beta_color_palette)
  colors <- colorRampPalette(base_cols)(n_cond)
  if (show_head && show_colors) {
    msg <- paste0(
      "\n--- colors (first few) ---\n",
      paste(utils::capture.output(utils::head(colors, head_n)), collapse = "\n"),
      "\n"
    )
    cat(msg, file = stderr()); flush.console()
  }
  log_head_vec(conds, "Condition levels")
  
  ## Boxplot of within-group distances
  sample_data(physeq1)$Group <- factor(sample_data(physeq1)$Condition, levels = conds)
  relab_genera <- transform_sample_counts(physeq1, function(x) x / sum(x) * 100)
  diveristy_distance <- phyloseq::distance(relab_genera, method = index_method)
  diveristy_distance <- as.matrix(diveristy_distance)
  log_head(diveristy_distance, paste0(index_method, " distances (relative, preview)"))
  
  sub_dist <- list()
  groups_all <- sample_data(relab_genera)$Group
  for (group in levels(groups_all)) {
    row_group <- which(groups_all == group)
    sample_group <- sample_names(relab_genera)[row_group]
    sub_dist[[group]] <- diveristy_distance[sample_group, sample_group, drop = FALSE]
    sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
  }
  
  braygroups <- melt(sub_dist)
  df.bray <- braygroups[complete.cases(braygroups), ]
  colnames(df.bray)[4] <- "Conditions"
  df.bray$Conditions <- factor(df.bray$Conditions, levels = names(sub_dist))
  log_head(df.bray, "df.bray (within-group distances)")
  
  plot1 <- ggplot(df.bray, aes(x = Conditions, y = value, colour = Conditions)) +
    geom_jitter() +
    geom_boxplot(alpha = 0.6) +
    ylab(paste(index_method, "_diversity")) +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
          axis.text.y  = element_text(size = 12)) +
    theme_bw() +
    scale_color_manual(values = colors)
  
  ## Ordination scatter
  plot2 <- ggplot(ordination_df, aes(x = Axis1, y = Axis2, color = Condition)) +
    geom_point(size = 3) +
    theme_minimal() +
    ggtitle(ordination_title) + theme_bw() +
    stat_ellipse(aes(group = Condition)) +
    scale_color_manual(values = colors)
  
  plots1 <- ggarrange(plot1, plot2, ncol = 2, nrow = 1)
  
  return(list(
    plot  = plots1,
    data1 = plot2$data,
    data2 = table_adonis
  ))
}
