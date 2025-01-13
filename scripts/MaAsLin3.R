library("maaslin3")
library("tibble")

MaAsLin3_summary <- function(OTU_input, group_index, group1, group2, select_MaAsLin3_normalization, select_MaAsLin3_transformation, select_MaAsLin3_correction, select_MaAsLin3_pvalue){

#data_tmp <- read.delim("~/subu/Metagenome/PRJNA385949/Comparison-Taxonomy-all_level.tsv", header = T, sep = "\t")
#data_species <- data_tmp[,c(7,8: ncol(data_tmp))]
#data_species_sum = aggregate(. ~ Species, data_species, sum)
#OTU_input <- data_species_sum[,2: ncol(data_species_sum)]
#rownames(OTU_input) <- data_species_sum[, 1] 

#group_index <- read.delim("~/subu/Metagenome/PRJNA385949/sample_metadata.tsv", header = T, sep = "\t")
#group1 ="CD"
#group2 = "HC"
#unique_conditions <- unique(group_index$Condition)
#colors <- colorRampPalette(brewer.pal(9, select_MaAsLin3_color_palette))(length(unique_conditions))

{
if(group1==group2)
{
  print("Please select different group  in condition1 or condition2")
  text_maaslin3_summary<- paste("Please select different group in condition1 or condition2")
  return(list(MaAsLin3_result_summary=text_maaslin3_summary))
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
    
    set.seed(1)
    tempdir <-getwd()
    output_dir <- paste0(getwd(),"/www/hmp2_output")
    fit_out <- maaslin3(input_data = OTU_input,
                        input_metadata = group_index1,
                        output = output_dir,
                        formula = 'Condition',
                        normalization = select_MaAsLin3_normalization,
                        transform = select_MaAsLin3_transformation,
                        correction = select_MaAsLin3_correction,
                        augment = TRUE,
                        standardize = TRUE,
                        plot_summary_plot = TRUE,
                        plot_associations = TRUE,
                        max_significance = select_MaAsLin3_pvalue,
                        median_comparison_abundance = TRUE,
                        median_comparison_prevalence = FALSE,
                        summary_plot_first_n = 25,
                        max_pngs = 30,
                        cores = 1)
    
    folder_to_zip <- paste0(getwd(),"/www/hmp2_output")
    output_zip <- paste0(getwd(),"/www/hmp2_output.zip")
    library(zip)
    zipr(output_zip, folder_to_zip)
text_maaslin3_summary<- paste("Click the 'Download Output as ZIP' button to download the results")
	
}
return(list(text_summary = text_maaslin3_summary, text_summary = text_maaslin3_summary))
}
}
