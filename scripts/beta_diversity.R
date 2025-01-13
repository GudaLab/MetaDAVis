library("phyloseq")
library("vegan")
library("reshape2")
library("ggplot2")
library("RColorBrewer")
library("ggpubr")

beta_diversity_boxplot <- function(OTU_input, group_index, tax_index, select_beta_color_palette, index_method, plot_method, adonis_permutations, adonis_dissimilarities) {
#data_tmp <- read.delim("~/subu/Metagenome/PRJNA385949/Comparison-Taxonomy-all_level.tsv", header = TRUE, sep = "\t")
#data_species <- data_tmp[, c(7, 8:ncol(data_tmp))]
#data_species_sum <- aggregate(. ~ Species, data_species, sum)
#OTU_input <- data_species_sum[, 2:ncol(data_species_sum)]
#rownames(OTU_input) <- data_species_sum[, 1]
#tax_index <- data_species_sum[1]
colnames(tax_index)[1] <- "Species"
rownames(tax_index) <- tax_index[, 1]

#group_index <- read.delim("~/subu/Metagenome/PRJNA385949/sample_metadata1.tsv", row.names = 1, header = TRUE, sep = "\t")
# Create phyloseq object
taxmat <- as.matrix(tax_index)
OTU <- otu_table(OTU_input, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
sampledata <- sample_data(group_index)
physeq1 <- phyloseq(OTU, TAX, sampledata)

# Extract abundance data
abundance_data <- as.data.frame(otu_table(physeq1))

# Ensure taxa are rows
if (!taxa_are_rows(physeq1)) {
  abundance_data <- t(abundance_data)
}

# Remove empty rows
empty_rows <- rowSums(abundance_data) == 0
if (any(empty_rows)) {
  abundance_data <- abundance_data[!empty_rows, ]
}

# Synchronize metadata and abundance data
shared_samples <- intersect(rownames(group_index), colnames(abundance_data))
if (length(shared_samples) == 0) {
  stop("Error: No overlapping samples between abundance_data and group_index.")
}

group_index <- group_index[shared_samples, , drop = FALSE]
abundance_data <- abundance_data[, shared_samples, drop = FALSE]

# Calculate Bray-Curtis distance matrix
bray_dist <- vegdist(t(abundance_data), method = index_method)

adonis_result <- adonis2(bray_dist ~ Condition, data = group_index, permutations = adonis_permutations, sqrt.dist = adonis_dissimilarities)

# Output results
#print(adonis_result)
table_adonis <- as.data.frame(adonis_result)

# Select method
#ordination_method <- "PCoA"  # Change to "NMDS", "DCA", "RDA", or "MDS" as needed

# Perform ordination based on selected method
ordination_result <- NULL
ordination_df <- NULL

if (plot_method == "PCoA") {
  pcoa_results <- cmdscale(bray_dist, eig = TRUE, k = 2)
  ordination_df <- as.data.frame(pcoa_results$points)
  colnames(ordination_df) <- c("Axis1", "Axis2")
  ordination_df$Condition <- group_index$Condition
  ordination_title <- "PCoA Ordination"
} else if (plot_method == "NMDS") {
  nmds_results <- metaMDS(bray_dist, k = 2)
  ordination_df <- as.data.frame(nmds_results$points)
  colnames(ordination_df) <- c("Axis1", "Axis2")
  ordination_df$Condition <- group_index$Condition
  ordination_title <- "NMDS Ordination"
} else if (plot_method == "DCA") {
  dca_results <- decorana(t(abundance_data))
  ordination_df <- as.data.frame(scores(dca_results, display = "sites"))
  colnames(ordination_df) <- c("Axis1", "Axis2")
  ordination_df$Condition <- group_index$Condition
  ordination_title <- "DCA Ordination"
} else if (plot_method == "RDA") {
  rda_results <- rda(t(abundance_data) ~ Condition, data = group_index)
  ordination_df <- as.data.frame(scores(rda_results, display = "sites"))
  colnames(ordination_df) <- c("Axis1", "Axis2")
  ordination_df$Condition <- group_index$Condition
  ordination_title <- "RDA Ordination"
} else if (plot_method == "MDS") {
  mds_results <- isoMDS(as.matrix(bray_dist))
  ordination_df <- as.data.frame(mds_results$points)
  colnames(ordination_df) <- c("Axis1", "Axis2")
  ordination_df$Condition <- group_index$Condition
  ordination_title <- "MDS Ordination"
} else {
  stop("Invalid ordination method selected. Choose from PCoA, NMDS, DCA, RDA, or MDS.")
}

colors <- colorRampPalette(brewer.pal(9, select_beta_color_palette))(length(unique(sample_data(physeq1)$Condition)))

# Boxplot
abundance_long <- as.data.frame(as.table(as.matrix(abundance_data)))
colnames(abundance_long) <- c("Taxa", "Sample", "Abundance")
abundance_long$Condition <- group_index$Condition[match(abundance_long$Sample, rownames(group_index))]

sample_data(physeq1)$Group <- factor((sample_data(physeq1)$Condition), levels=unique(sampledata$Condition))
relab_genera = transform_sample_counts(physeq1, function(x) x / sum(x) * 100) 
diveristy_distance <- phyloseq::distance(relab_genera, method = index_method)
diveristy_distance <- as.matrix(diveristy_distance)
sub_dist <- list()
groups_all <- sample_data(relab_genera)$Group
for (group in levels(groups_all)) { 
  row_group <- which(groups_all == group)
  sample_group <- sample_names(relab_genera)[row_group]
  sub_dist[[group]] <- diveristy_distance[ sample_group, sample_group]
  sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
}

braygroups<- melt(sub_dist)
df.bray <- braygroups[complete.cases(braygroups), ]
colnames(df.bray)[4]  <- "Conditions"

df.bray$Conditions <- factor(df.bray$Conditions, levels=names(sub_dist))
plot1 <- ggplot(df.bray, aes(x=Conditions, y=value, colour=Conditions)) + geom_jitter() + geom_boxplot(alpha=0.6) + ylab(paste(index_method,"_diversity")) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12), axis.text.y=element_text(size=12)) + theme_bw()+
  scale_color_manual(values = colors)
  
#plot1<-ggplot(abundance_long, aes(x = Condition, y = Abundance)) +
#  geom_boxplot() +
#  theme_minimal() +
#  ggtitle("Abundance by Condition (Boxplot)") +
#  theme_bw() +
# scale_color_manual(values = colors)

# Plot ordination
plot2<-ggplot(ordination_df, aes(x = Axis1, y = Axis2, color = Condition)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle(ordination_title)+theme_bw() +
  stat_ellipse(aes(group=Condition)) +
  scale_color_manual(values = colors)

plots1<-ggarrange(plot1, plot2, ncol = 2, nrow = 1)

return(list(plot = plots1, data1=plot2$data, data2=table_adonis))
}
