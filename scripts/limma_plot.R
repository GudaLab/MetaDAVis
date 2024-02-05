library("limma")
library("ggplot2")
library("tibble")
library("ComplexHeatmap")
library("scales")
library("ggplotify")
limma_plot <- function(OTU_input, group_index, index_pvalue, group1, group2, plot_method, alpha){
  {
if(group1==group2)
{
  print("Please select different group  in condition1 or condition2")
  text_limma_summary<- paste("Please select different group in condition1 or condition2")
  return(list(limma_result_summary=text_limma_summary))
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
  
  limma.design <- model.matrix(~factor(setupinfo$Condition, levels = c(group1, group2)))
  
  y <- DGEList(counts = OTU_input, group = factor(setupinfo$Condition, levels = c(group1, group2)))
  y <- edgeR::calcNormFactors(y)
  remove <- is.na(y$samples$group)
  y <- y[, !remove]
  v <- voom(y, limma.design, plot = FALSE)
  fit <- lmFit(v, limma.design)
  fit <- eBayes(fit) #this is the hierarchical test procedure
  
  results <- topTable(fit, coef = 2, number = nrow(fit), adjust.method = "BH", sort.by = "none")
  results$logFC <- results$logFC
  results$foldChange <- 2 ^ results$logFC
  res_tax = cbind(as.data.frame(results), OTU = rownames(results))
  res_tax <- res_tax[, c("OTU","foldChange", "logFC", "P.Value", "adj.P.Val")]
  
  colnames(res_tax) <- c("OTU", "foldChange", "log2FC", "PValue", "FDR")
  setupinfo <- setupinfo[complete.cases(setupinfo), ]
  group.exp <- round(2 ^ (as.data.frame.EList(v, row.names = rownames(1:nrow(OTU_input)))), digits = 0)
  
  Group1.Mean <- rowMeans(group.exp[, which(setupinfo$Condition == group1)])
  Group2.Mean <- rowMeans(group.exp[, which(setupinfo$Condition == group2)])
  All_Mean <- (Group1.Mean + Group2.Mean) / 2
  Means <- cbind(Group1.Mean, Group2.Mean, All_Mean)
  
  colnames(Means) <- c(paste0("Mean_", group1), paste0("Mean_", group2), "All_Mean")
  
  res_tax <- res_tax[order(res_tax$FDR, res_tax$PValue), ]
  
  #qv <- qvalue(res_tax$PValue)
  #res_tax$qValue <- qv$qvalues
  #res_tax$Local.FDR <- qv$lfdr
  
  res_tax <- cbind(res_tax,Means)
  
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
      text_limma_summary<- paste("Total of ", nrow(res_tax_sig), " taxa were identified as significant.", sep="")
}
else if(dim(res_tax_sig)[1]== 0)
  {
  text_limma_summary<- paste("No. significant were identified. Adjust the FDR value or change the test correction method")
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
    text_limma_summary<- paste("Total of ", nrow(res_tax_sig), " taxa were identified as significant.", sep="")
}
else if(dim(res_tax_sig)[1]== 0)
  {
  text_limma_summary<- paste("No. significant were identified. Adjust the Pvalue value or change the test correction method")
  #print("No. significant were identified")
}
}
}
    else
    {
      print("error")
    }
  }

    
  
  {
    if(plot_method == "1"){
    p1 <-ggplot(data = df, aes(x=Taxa1, y=Value, color=Condition)) + 
    geom_boxplot() +
    #scale_y_log10() +
    ggtitle("Significant Taxa's") + 
    xlab("Taxa's") +ylab("Log-relative normalised") +
    theme_bw() +
    theme(legend.position = "right")+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
	p1
    }
    else if(plot_method == "2")
    {
    p2 <- ggplot(df,aes(Condition,Value,colour=Condition))+
    ylab("Log-relative normalised")+
    geom_boxplot()+geom_jitter()+theme_bw()+
    facet_wrap( ~ Taxa , scales="free", ncol=3)+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
	p2
    }
	else if(plot_method == "3")
    {
	p3 <- ggplot(data = res_tax, aes(x = All_Mean, y = log2FC, color = Significant)) +
    geom_point(size = 2) +
    scale_x_log10() +
    scale_color_manual(values=c("black", "red")) +
    labs(x = "log10(Mean abundance)", y = "Log2FC")+theme_bw()+  
    geom_text(data = subset(res_tax, Significant == "Yes"), aes(label = OTU), size = 4, vjust = 1)
    p3
    }
	else if(plot_method == "4")
   {
   ha = HeatmapAnnotation(Group = group_index$Condition, annotation_height = unit(4, "mm"))
   p4<-Heatmap(df1, name = "Scale", show_row_names = TRUE, row_names_gp = gpar(fontsize = 7), show_column_names = TRUE, column_names_gp = gpar(fontsize = 7), show_row_dend = TRUE, show_column_dend = FALSE, cluster_columns = FALSE, top_annotation = ha)
   q <- as.ggplot(p4)
   q
    }
    }
}
}
}
