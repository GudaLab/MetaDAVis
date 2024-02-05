library("ggplot2")
library("tibble")
library("qvalue")

anova_summary <- function(OTU_input, group_index, index_pvalue, ad_hoc,  plot_method, alpha){
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
    test_tmp <- aov(as.numeric(OTU_tmp)~as.factor(group_index))
    for(j in 1:length(names(table(group_index)))){
      test_summary[i, paste("Mean_relative_frequency", names(table(group_index))[j], sep = "_")] <- signif(mean(as.numeric(OTU_tmp[which(group_index == names(table(group_index))[j])])), 4)
      test_summary[i, paste("Present_in_no_of", names(table(group_index))[j], sep = "_")] <- length(which(as.numeric(OTU_tmp[which(group_index == names(table(group_index))[j])]) != 0))
    }
    test_summary[i, "OTU"] <- rownames(OTU_input)[i]
    test_summary[i, "PValue"] <- signif(summary(test_tmp)[[1]][["Pr(>F)"]][1], 4)
  }
}

test_summary2 <- test_summary[!is.na(test_summary[,"PValue"]), ]
test_summary2[ , "q_value"] <- signif(p.adjust(as.numeric(test_summary2[ ,"PValue"]), "BH"), 4)

#TukeyHSD
pairwise_p_matrix <- matrix(NA, nrow = nrow(test_summary2), ncol = cumsum(1: (length(names(table(group_index))) - 1))[(length(names(table(group_index))) - 1)])
for(i in 1: nrow(test_summary2)){
  OTU_tmp <- OTU_input[match(test_summary2[i, "OTU"], rownames(OTU_input)), ]
  #p_matrix_tmp2 <- dunn.test(as.numeric(OTU_tmp), group_index, method = "bh")
  p_matrix_tmp <- TukeyHSD((aov(as.numeric(OTU_tmp)~as.factor(group_index))), conf.level = 0.95)
  pairwise_p_tmp <- signif(p_matrix_tmp$`as.factor(group_index)`[,"p adj"], 4)
  names(pairwise_p_tmp) <- names(p_matrix_tmp$`as.factor(group_index)`[,1])
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
text_anova_summary<- paste("Total of ", nrow(res_tax_sig), " taxa were identified as significant.", sep="")
}
else if(dim(res_tax_sig)[1]== 0)
  {
  text_anova_summary<- paste("No. significant were identified. Adjust the FDR value or change the test correction method")
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
text_anova_summary<- paste("Total of ", nrow(res_tax_sig), " taxa were identified as significant.", sep="")
}
else if(dim(res_tax_sig)[1]== 0)
  {
  text_anova_summary<- paste("No. significant were identified. Adjust the Pvalue value or change the test correction method")
  #print("No. significant were identified")
}
}
}
    else
    {
      print("error")
    }
  }
return(list(anova_result_summary=text_anova_summary, anova_result_significant=res_tax_sig, anova_result_table=res_tax, anova_relative_frequency=rel_freq, total_counts=parent_seq_count))

 }
 else
    {
	print("If you have less than two condition or group use anyone of the two groups method for your analysis")
    text_anova_summary<- paste("If you have less than two condition or group use anyone of the two groups method for your analysis")
	return(list(anova_result_summary=text_anova_summary))

    }
	}
	
	


}	