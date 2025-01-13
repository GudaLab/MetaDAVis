library("lefser")
library("tibble")

LEfSe_summary <- function(OTU_input, group_index, group1, group2, select_LEfSe_method, select_LEfSe_pvalue, select_LEfSe_threshold, select_LEfSe_color_palette){

#data_tmp <- read.delim("~/subu/Metagenome/PRJNA385949/Comparison-Taxonomy-all_level.tsv", header = T, sep = "\t")
#data_species <- data_tmp[,c(7,8: ncol(data_tmp))]
#data_species_sum = aggregate(. ~ Species, data_species, sum)
#OTU_input <- data_species_sum[,2: ncol(data_species_sum)]
#rownames(OTU_input) <- data_species_sum[, 1] 

#group_index <- read.delim("~/subu/Metagenome/PRJNA385949/sample_metadata.tsv", header = T, sep = "\t")
#group1 ="CD"
#group2 = "HC"
unique_conditions <- unique(group_index$Condition)
colors <- colorRampPalette(brewer.pal(9, select_LEfSe_color_palette))(length(unique_conditions))

{
if(group1==group2)
{
  print("Please select different group  in condition1 or condition2")
  text_deseq_summary<- paste("Please select different group in condition1 or condition2")
  return(list(LEfSe_result_summary=text_deseq_summary))
}
else if(group1!=group2)
{
group_index <- subset(group_index, ((group_index$Condition == group1) | (group_index$Condition == group2)))

OTU_input <- as.data.frame(t(OTU_input))


#convert row name to 1st column
OTU_input <- tibble::rownames_to_column(OTU_input, "Samples")
#mapping column based on metadata
OTU_input <- OTU_input[OTU_input$Samples %in% group_index$Samples, ]
#convert 1st column to row name
OTU_input <- OTU_input %>% remove_rownames %>% column_to_rownames(var="Samples")

OTU_input <- as.data.frame(t(OTU_input))
OTU_input <- OTU_input[rowSums(OTU_input[])>0,]

#count total
parent_seq_count <- as.data.frame(colSums(OTU_input))
parent_seq_count <- tibble::rownames_to_column(parent_seq_count, "Samples")
colnames(parent_seq_count)[2] <- "Total_counts"


group_index1 <- group_index
group_index1<- group_index1 %>% remove_rownames %>% column_to_rownames(var="Samples")


# Create phyloseq object
tax_index <- as.data.frame(rownames(OTU_input))
colnames(tax_index)[1] <- "Species"
rownames(tax_index) <- tax_index[, 1]
taxmat <- as.matrix(tax_index)
OTU <- otu_table(OTU_input, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
sampledata <- sample_data(group_index1)
physeq1 <- phyloseq(OTU, TAX, sampledata)


counts <- unclass(otu_table(physeq1))
coldata <- as(sample_data(physeq1), "data.frame")
SummarizedExperiment(assays = list(counts = counts), colData = coldata)
lefse_input<-mia::convertFromPhyloseq(physeq1)

res <- lefser(lefse_input, # relative abundance only with terminal nodes
              classCol = "Condition",
              kruskal.threshold = 0.05,
              lda.threshold = 2,
              method="BH"
)
plots1<-lefserPlot(res, colors=colors)+
theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "line")
  )
text_deseq_summary<- paste("Total of ", nrow(plots1$data), " taxa were identified as significant.", sep="")

return(list(LEfSe_result_summary=text_deseq_summary, data1=plots1$data, data2=parent_seq_count, plot=plots1))
}
}
}
