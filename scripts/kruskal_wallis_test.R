library("ggplot2")
library("tibble")
library("qvalue")
library("dunn.test")
kruskal_wallis_test_summary <- function(OTU_input, group_index, index_pvalue, ad_hoc, plot_method, alpha, kruskal_wallis_test_color_palette){
unique_conditions <- unique(group_index$Condition)
colors <- colorRampPalette(brewer.pal(9, kruskal_wallis_test_color_palette))(length(unique_conditions))

OTU_input <- OTU_input[rowSums(OTU_input[])>0,]
OTU_input1 <- OTU_input[rowSums(OTU_input[])>0,]

#count total
parent_seq_count <- as.data.frame(colSums(OTU_input))
parent_seq_count <- tibble::rownames_to_column(parent_seq_count, "Samples")
colnames(parent_seq_count)[2] <- "Total_counts"

#calculate relative frequency
rel_freq <- (OTU_input)/(colSums(OTU_input))*100

#write.csv(OTU_input, file = "OTU_input1.csv")

OTU_input<-as.data.frame(rel_freq)
group_index1 <- group_index
group_index <- group_index[,2]
{
if(length(names(table(group_index))) > 2){


test_summary <- matrix(NA, nrow = nrow (OTU_input), ncol = 3 + (length(names(table(group_index)))*2) )
colnames(test_summary) <- c("OTU", paste("Present_in_no_of", names(table(group_index)), sep = "_"), paste("Mean_relative_frequency", names(table(group_index)), sep = "_"), "PValue", "q_value")
for(i in 1: nrow(OTU_input)){
  OTU_tmp <- OTU_input[i, ]
   if(all(sapply(1: length(names(table(group_index))), function(k) sd(OTU_tmp[which(group_index == names(table(group_index))[k])])) >=0)){
    test_tmp <- kruskal.test(as.numeric(OTU_tmp)~as.factor(group_index))
    for(j in 1:length(names(table(group_index)))){
      test_summary[i, paste("Mean_relative_frequency", names(table(group_index))[j], sep = "_")] <- signif(mean(as.numeric(OTU_tmp[which(group_index == names(table(group_index))[j])])), 4)
      test_summary[i, paste("Present_in_no_of", names(table(group_index))[j], sep = "_")] <- length(which(as.numeric(OTU_tmp[which(group_index == names(table(group_index))[j])]) != 0))
    }
    test_summary[i, "OTU"] <- rownames(OTU_input)[i]
    test_summary[i, "PValue"] <- signif(test_tmp$p.value, 4)
  }
}

test_summary2 <- test_summary[!is.na(test_summary[,"PValue"]), ]
test_summary2[ , "q_value"] <- signif(p.adjust(as.numeric(test_summary2[ ,"PValue"]), "BH"), 4)

#Dunn.test
pairwise_p_matrix <- matrix(NA, nrow = nrow(test_summary2), ncol = cumsum(1: (length(names(table(group_index))) - 1))[(length(names(table(group_index))) - 1)])
for(i in 1: nrow(test_summary2)){
  OTU_tmp <- OTU_input[match(test_summary2[i, "OTU"], rownames(OTU_input)), ]
  p_matrix_tmp <- dunn.test(as.numeric(OTU_tmp), group_index, method = "bh")
  pairwise_p_tmp <- signif(p_matrix_tmp$P, 4)
  names(pairwise_p_tmp) <- p_matrix_tmp$comparisons
  pairwise_p_matrix[i, ] <- pairwise_p_tmp
}
colnames(pairwise_p_matrix) <- names(pairwise_p_tmp)

if(ad_hoc == "Yes"){
test_summary2 <- cbind(test_summary2, pairwise_p_matrix)
res_tax <- as.data.frame(test_summary2)
rownames(res_tax) <- res_tax[, 1] 
}
else if(ad_hoc == "No"){
res_tax <- as.data.frame(test_summary2)
rownames(res_tax) <- res_tax[, 1] 
}
  {
    if(index_pvalue == "padj"){
res_tax_sig = subset(res_tax, q_value < alpha)
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
text_kruskal_wallis_test_summary<- paste("Total of ", nrow(res_tax_sig), " taxa were identified as significant.", sep="")
text_kruskal_wallis_test_summary<- paste("Total of ", nrow(res_tax_sig), " taxa were identified as significant.", sep="")
}
else if(dim(res_tax_sig)[1]== 0)
  {
  text_kruskal_wallis_test_summary<- paste("No. significant were identified. Adjust the FDR value or change the test correction method")
  #print("No. significant were identified")
}
}
}
   
    else if(index_pvalue == "pvalue"){
res_tax_sig = subset(res_tax, PValue < alpha)
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
text_kruskal_wallis_test_summary<- paste("Total of ", nrow(res_tax_sig), " taxa were identified as significant.", sep="")
text_kruskal_wallis_test_summary<- paste("Total of ", nrow(res_tax_sig), " taxa were identified as significant.", sep="")
}
else if(dim(res_tax_sig)[1]== 0)
  {
  text_kruskal_wallis_test_summary<- paste("No. significant were identified. Adjust the Pvalue value or change the test correction method")
  #print("No. significant were identified")
}
}
}
    else
    {
      print("error")
    }
  }
  if(plot_method == "1"){
    p1<-ggplot(data = df, aes(x=Taxa1, y=Value, color=Condition)) + 
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
   p1 <- ggplot(df,aes(Condition,Value,colour=Condition))+
   ylab("Relative frequency")+
   geom_boxplot()+geom_jitter()+theme_bw()+
   facet_wrap( ~ Taxa , scales="free", ncol=3)+
   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
      scale_color_manual(values = colors)
   #p2
      }
	else if(plot_method == "3")
   {
   #ha = HeatmapAnnotation(Group = group_index1$Condition, annotation_height = unit(4, "mm"))
   ha <- HeatmapAnnotation(
      Group = group_index1$Condition,
      col = list(Group = structure(brewer.pal(length(unique_conditions), kruskal_wallis_test_color_palette), names = unique_conditions))
    )
   p3<-Heatmap(df1, name = "Scale", show_row_names = TRUE, row_names_gp = gpar(fontsize = 7), show_column_names = TRUE, column_names_gp = gpar(fontsize = 7), show_row_dend = TRUE, show_column_dend = FALSE, cluster_columns = FALSE, top_annotation = ha, col = colors)
   p1 <- as.ggplot(p3)
   #q
    }
return(list(kruskal_wallis_test_result_summary=text_kruskal_wallis_test_summary, kruskal_wallis_test_result_significant=res_tax_sig, kruskal_wallis_test_result_table=res_tax, kruskal_wallis_test_relative_frequency=rel_freq, total_counts=parent_seq_count, plot=p1))

 }
 else
    {
    text_kruskal_wallis_test_summary <-  paste("If you have less than two condition or group use anyone of the two groups method for your analysis")
    return(list(kruskal_wallis_test_result_summary=text_kruskal_wallis_test_summary))

	}
	}
	
	


}	