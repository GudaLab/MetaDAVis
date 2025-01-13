library("ggplot2")
library("tibble")
library("qvalue")
library("ComplexHeatmap")
library("scales")
library("ggplotify")
wilcoxtest_summary <- function(OTU_input, group_index, index_pvalue, group1, group2, plot_method, alpha, wilcox_color_palette){
#data_tmp <- read.delim("~/subu/Metagenome/PRJNA385949/Comparison-Taxonomy-all_level.tsv", header = T, sep = "\t")
#data_species <- data_tmp[,c(7,8: ncol(data_tmp))]
#data_species_sum = aggregate(. ~ Species, data_species, sum)
#OTU_input <- data_species_sum[,2: ncol(data_species_sum)]
#rownames(OTU_input) <- data_species_sum[, 1] 
#group_index <- read.delim("~/subu/Metagenome/PRJNA385949/sample_metadata.tsv", header = T, sep = "\t")
#group1 ="CD"
#group2 = "HC"

 
{
if(group1==group2)
{
  print("Please select different group  in condition1 or condition2")
  text_wilcoxtest_summary<- paste("Please select different group in condition1 or condition2")
  return(list(wilcoxtest_result_summary=text_wilcoxtest_summary))
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
OTU_input1 <- OTU_input[rowSums(OTU_input[])>0,]

#count total
parent_seq_count <- as.data.frame(colSums(OTU_input))
parent_seq_count <- tibble::rownames_to_column(parent_seq_count, "Samples")
colnames(parent_seq_count)[2] <- "Total_counts"

#calculate relative frequency
rel_freq <- (OTU_input)/(colSums(OTU_input))*100
OTU_input<-as.data.frame(rel_freq)

# adding suffix to column names 
colnames(rel_freq) <- paste(colnames(rel_freq),"rel_freq",sep="_")
group_index1 <- group_index
group_index <- group_index[,2]
test_summary <- matrix(NA, nrow = nrow (OTU_input), ncol = 11)
colnames(test_summary) <- c("OTU", paste("Present_in_no_of", names(table(group_index)), sep = "_"), paste("Mean_relative_frequency", names(table(group_index)), sep = "_"), "All_mean_relative_frequency","Difference_between_means", "fold_change", "log2FC", "PValue", "q_value")
for(i in 1: nrow(OTU_input)){
  OTU_tmp <- OTU_input[i, ]
  test_tmp <- wilcox.test(as.numeric(OTU_tmp[which(group_index == names(table(group_index))[1])]),  as.numeric(OTU_tmp[which(group_index == names(table(group_index))[2])]))
  test_summary[i, paste("Mean_relative_frequency", names(table(group_index))[1], sep = "_")] <- signif(mean(as.numeric(OTU_tmp[which(group_index == names(table(group_index))[1])])), 4)
  test_summary[i, paste("Mean_relative_frequency", names(table(group_index))[2], sep = "_")] <- signif(mean(as.numeric(OTU_tmp[which(group_index == names(table(group_index))[2])])), 4)
  test_summary[i, paste("All_mean_relative_frequency")] <- ((signif(mean(as.numeric(OTU_tmp[which(group_index == names(table(group_index))[1])]))))+(signif(mean(as.numeric(OTU_tmp[which(group_index == names(table(group_index))[2])]))))/2)
  test_summary[i, paste("Difference_between_means")] <- (signif(mean(as.numeric(OTU_tmp[which(group_index == names(table(group_index))[1])]))))-(signif(mean(as.numeric(OTU_tmp[which(group_index == names(table(group_index))[2])]))))
  test_summary[i, paste("fold_change")] <- (signif(mean(as.numeric(OTU_tmp[which(group_index == names(table(group_index))[2])]))))/(signif(mean(as.numeric(OTU_tmp[which(group_index == names(table(group_index))[1])]))))
  test_summary[i, paste("log2FC")] <- log2((signif(mean(as.numeric(OTU_tmp[which(group_index == names(table(group_index))[2])]))))/(signif(mean(as.numeric(OTU_tmp[which(group_index == names(table(group_index))[1])])))))
  test_summary[i, paste("Present_in_no_of", names(table(group_index))[1], sep = "_")] <- length(which(as.numeric(OTU_tmp[which(group_index == names(table(group_index))[1])]) != 0))
  test_summary[i, paste("Present_in_no_of", names(table(group_index))[2], sep = "_")] <- length(which(as.numeric(OTU_tmp[which(group_index == names(table(group_index))[2])]) != 0))
  test_summary[i, "OTU"] <- rownames(OTU_input)[i]
  test_summary[i, "PValue"] <- signif(test_tmp$p.value, 4)
}
test_summary2 <- test_summary[!is.na(test_summary[,"PValue"]), ]
test_summary2[ , "q_value"] <- signif(p.adjust(as.numeric(test_summary2[ ,"PValue"]), "BH"), 4)
res_tax <- as.data.frame(test_summary2)
rownames(res_tax) <- res_tax[, 1] 
res_tax[res_tax=="-Inf"]<-0
res_tax[res_tax=="Inf"]<-0
#alpha = 0.05
fold = 1

  {
    if(index_pvalue == "padj"){
res_tax_sig = subset(res_tax, q_value < alpha & fold < abs(as.numeric(log2FC)))
{
if (dim(res_tax_sig)[1] >= 1) {
#mark significant
res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig) , "Yes", "No")
OTU_input2 <- OTU_input1[rownames(OTU_input1) %in% rownames(res_tax_sig),]
df1 <- scale(OTU_input2)

#Now we plot taxa significantly different between the categories
df<-NULL

res_tax$q_value<-as.numeric(res_tax$q_value)
data <- as.data.frame(t(OTU_input))

for(i in res_tax[rownames(res_tax_sig),"OTU"]){
  tmp<-data.frame(data[,i],group_index,rep(paste(i)),rep(paste(i," q_value = ",round(res_tax[i,"q_value"],4),sep=""),dim(data)[1]))
  if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)} 
}
colnames(df)<-c("Value","Condition","Taxa1","Taxa")  
text_wilcoxtest_summary<- paste("Total of ", nrow(res_tax_sig), " taxa were identified as significant.", sep="")
}
else if(dim(res_tax_sig)[1]== 0)
  {
  text_wilcoxtest_summary<- paste("No. significant were identified. Adjust the FDR value or change the test correction method")
  #print("No. significant were identified")
}
}
}

  else if(index_pvalue == "pvalue"){
   
res_tax_sig = subset(res_tax, PValue < alpha & fold < abs(as.numeric(log2FC)))
{
if (dim(res_tax_sig)[1] >= 1) {
#mark significant
res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig) , "Yes", "No")
OTU_input2 <- OTU_input1[rownames(OTU_input1) %in% rownames(res_tax_sig),]
df1 <- scale(OTU_input2)

#Now we plot taxa significantly different between the categories
df<-NULL

res_tax$PValue<-as.numeric(res_tax$PValue)
data <- as.data.frame(t(OTU_input))

for(i in res_tax[rownames(res_tax_sig),"OTU"]){
  tmp<-data.frame(data[,i],group_index,rep(paste(i)),rep(paste(i," PValue = ",round(res_tax[i,"PValue"],4),sep=""),dim(data)[1]))
  if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)} 
}
colnames(df)<-c("Value","Condition","Taxa1","Taxa") 
text_wilcoxtest_summary<- paste("Total of ", nrow(res_tax_sig), " taxa were identified as significant.", sep="")
}
else if(dim(res_tax_sig)[1]== 0)
  {
  text_wilcoxtest_summary<- paste("No. significant were identified. Adjust the Pvalue value or change the test correction method")
  #print("No. significant were identified")
}
}
}
else
      {
        print("error")
      }
  }
# }
 unique_conditions <- unique(group_index)
colors <- colorRampPalette(brewer.pal(9, wilcox_color_palette))(length(unique_conditions))
 
   if(plot_method == "1"){
    p1 <- ggplot(data = df, aes(x=Taxa1, y=Value, color=Condition)) + 
    geom_boxplot() +
    scale_y_log10() +
    ggtitle("Significant Taxa's") + 
    xlab("Taxa's") +ylab("log10(Relative frequency)") +
    theme_bw() +
    theme(legend.position = "right")+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
      scale_color_manual(values = colors)
    #p1
   }
   else if(plot_method == "2")
   {
   p1<-ggplot(df,aes(Condition,Value,colour=Condition))+
   ylab("Relative frequency")+
   geom_boxplot()+geom_jitter()+theme_bw()+
   facet_wrap( ~ Taxa , scales="free", ncol=3)+
   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
      scale_color_manual(values = colors)
   #p2
   }
   else if(plot_method == "3")
   {
  p1<-ggplot(data = res_tax, aes(x = as.numeric(All_mean_relative_frequency), y = as.numeric(log2FC), color = Significant)) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_color_manual(values=c("black", "red")) +
  labs(x = "log10(Mean relative abundance)", y = "Log2FC")+theme_bw()+  
  geom_text(data = subset(res_tax, Significant == "Yes"), aes(label = OTU), size = 4, vjust = 1)
  #p3
   }
   else if(plot_method == "4")
   {
   #ha = HeatmapAnnotation(Group = group_index1$Condition, annotation_height = unit(4, "mm"))
   color_mapping <- setNames(colors, unique_conditions)
   ha <- HeatmapAnnotation(
  Group = group_index,
  col = list(Group = color_mapping)
)
   p4<-Heatmap(df1, name = "Scale", show_row_names = TRUE, row_names_gp = gpar(fontsize = 7), show_column_names = TRUE, column_names_gp = gpar(fontsize = 7), show_row_dend = TRUE, show_column_dend = FALSE, cluster_columns = FALSE, top_annotation = ha, col = colors)
   p1 <- as.ggplot(p4)
   #q
   }
 return(list(wilcoxtest_result_summary=text_wilcoxtest_summary, wilcoxtest_result_significant=res_tax_sig, wilcoxtest_result_table=res_tax, wilcoxtest_relative_frequency=rel_freq, total_counts=parent_seq_count, plot=p1))
 }
 }
}