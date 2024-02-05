library("ggplot2")
library("tibble")
library("edgeR")

edger_summary <- function(OTU_input, group_index, index_pvalue, group1, group2, plot_method, alpha){
  {
if(group1==group2)
{
  print("Please select different group  in condition1 or condition2")
  text_edger_summary<- paste("Please select different group in condition1 or condition2")
  return(list(edger_result_summary=text_edger_summary))
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
  
  
  setupinfo <- data.frame("sample_ID" = group_index$Samples, "Condition" = group_index$Condition)
  setupinfo$Condition <- factor(setupinfo$Condition, levels <- c(group1, group2))
  #alpha = 0.05
  fold = 1
  
  design <- model.matrix(~factor(setupinfo$Condition, levels = c(group1, group2)))
  
  dge <- DGEList(OTU_input, group = factor(setupinfo$Condition, levels = c(group1, group2)))
  remove <- is.na(dge$samples$group)
  dge <- dge[, !remove]
  dge <- edgeR::calcNormFactors(dge)
  dge <- estimateGLMRobustDisp(dge, design = design)
  
  fit <- glmFit(dge, design = design)
  lrt <- glmLRT(fit)
  
  results <- topTags(lrt, n = Inf)
  results <- results$table
  results$foldChange <- 2 ^ results$logFC
  res_tax = cbind(as.data.frame(results), OTU = rownames(results))
  res_tax <- res_tax[, c("OTU","foldChange", "logFC", "PValue", "FDR")]
  colnames(res_tax) <- c("OTU","foldChange", "log2FC", "PValue", "FDR")
  
  res_tax <- res_tax[order(res_tax$FDR, res_tax$PValue), ]
  #qv <- qvalue(res_tax$PValue)
  
  {
    if(index_pvalue == "padj"){
      #filter significant
      res_tax_sig = subset(res_tax, FDR < alpha & fold < abs(log2FC))
      {
	  if (dim(res_tax_sig)[1] >= 1) {
      #mark significant
      res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig) , "Yes", "No")
      OTU_input1 <- OTU_input[rownames(OTU_input) %in% rownames(res_tax_sig),]
      df1 <- scale(OTU_input1)
	  
      ##Apply normalisation (either use relative or log-relative transformation)
      abund_table<-t(OTU_input)
      data<-log((abund_table+1)/(rowSums(abund_table)+dim(abund_table)[2]))
      data<-as.data.frame(data)
      
      #Now we plot taxa significantly different between the categories
      df<-NULL
      groups<-group_index[,2]
      
      for(i in res_tax[rownames(res_tax_sig),"OTU"]){
        tmp<-data.frame(data[,i],groups,rep(paste(i)),rep(paste(i," FDR = ",round(res_tax[i,"FDR"],4),sep=""),dim(data)[1]))
        if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)} 
      }
      colnames(df)<-c("Value","Condition","Taxa1","Taxa")  
      text_edger_summary<- paste("Total of ", nrow(res_tax_sig), " taxa were identified as significant.", sep="")
}
else if(dim(res_tax_sig)[1]== 0)
  {
  text_edger_summary<- paste("No. significant were identified. Adjust the FDR value or change the test correction method")
  #print("No. significant were identified")
}
}
}
    
    else if(index_pvalue == "pvalue"){
      #filter significant
      res_tax_sig = subset(res_tax, PValue < alpha & fold < abs(log2FC))
      {
      if (dim(res_tax_sig)[1] >= 1) {
      #mark significant
      res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig) , "Yes", "No")
	  OTU_input1 <- OTU_input[rownames(OTU_input) %in% rownames(res_tax_sig),]
      df1 <- scale(OTU_input1)
      
      ##Apply normalisation (either use relative or log-relative transformation)
      abund_table<-t(OTU_input)
      data<-log((abund_table+1)/(rowSums(abund_table)+dim(abund_table)[2]))
      data<-as.data.frame(data)
      
      #Now we plot taxa significantly different between the categories
      df<-NULL
      groups<-group_index[,2]
      
      for(i in res_tax[rownames(res_tax_sig),"OTU"]){
        tmp<-data.frame(data[,i],groups,rep(paste(i)),rep(paste(i," PValue = ",round(res_tax[i,"PValue"],4),sep=""),dim(data)[1]))
        if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)} 
      }
      colnames(df)<-c("Value","Condition","Taxa1","Taxa") 
    text_edger_summary<- paste("Total of ", nrow(res_tax_sig), " taxa were identified as significant.", sep="")
}
else if(dim(res_tax_sig)[1]== 0)
  {
  text_edger_summary<- paste("No. significant were identified. Adjust the Pvalue value or change the test correction method")
  #print("No. significant were identified")
}
}
}
    else
    {
      print("error")
    }
  }
   #text_edger_summary<- paste("Total of ", nrow(res_tax_sig), " taxa were identified as significant.", sep="")
  }
  return(list(edger_result_summary=text_edger_summary, edger_result_significant=res_tax_sig, edger_result_table=res_tax, total_counts=parent_seq_count))
}
}
  
  