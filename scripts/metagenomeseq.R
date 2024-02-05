library("metagenomeSeq")
library("ggplot2")
library("tibble")

metagenomeseq_summary <- function(OTU_input, group_index, index_pvalue, group1, group2, plot_method, alpha){

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
  text_metagenomeseq_summary<- paste("Please select different group in condition1 or condition2")
  return(list(metagenomeseq_result_summary=text_metagenomeseq_summary))
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
  
  setupinfo <- data.frame("sample_ID" = group_index$Samples, "Condition" = group_index$Condition)
  setupinfo$Condition <- factor(setupinfo$Condition, levels <- c(group1, group2))
  setupinfo <- setupinfo %>% remove_rownames %>% column_to_rownames(var="sample_ID")
  

  OTU_input <- OTU_input[, colnames(OTU_input) %in% rownames(setupinfo)]
  OTU_input <- OTU_input[, order(colnames(OTU_input))]
  # rownames(OTU_input) <- x$OTU
  obj <- newMRexperiment(OTU_input, phenoData = new("AnnotatedDataFrame", data = setupinfo))
  
  # phenoData(obj) = AnnotatedDataFrame(setupinfo) #new("AnnotatedDataFrame", data = setupinfo))
  p <- cumNormStatFast(obj)
  obj <- cumNorm(obj, p = p)
  
  mod <- model.matrix(~factor(setupinfo$Condition, levels = c(group1, group2)), data = pData(obj))
  res <- fitFeatureModel(obj = obj, mod = mod, coef = 2)
  
  results <- MRtable(obj = res, number = Inf, by = 2, adjustMethod = "BH")
  # results$logFC <- results$logFC
  results$foldChange <- 2 ^ results$logFC
  All_Mean <- (results$`counts in group 0` + results$`counts in group 1`) /  nrow(setupinfo)
  results <- cbind(results, All_Mean)
  res_tax <- data.frame(OTU = rownames(results), results)
  colnames(res_tax) <- c("OTU", paste0("Present_in_no_of_", group1), paste0("Present_in_no_of_", group2), 
                         paste0("Counts_in_", group1), paste0("Counts_in_", group2), "log2FC", "standard_error", "PValue", "FDR", "foldChange", "All_Mean")

  res_tax <- res_tax[order(res_tax$FDR, res_tax$PValue), ]
  
  #name = paste0(saveprefix, ".", group1, ".vs.", group2, ".", saveprefix, ".res_tax.csv")
  
  #qv <- qvalue(res_tax$PValue)
  
  #res_tax$qValue <- qv$qvalues
  
fold=1
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
    text_metagenomeseq_summary<- paste("Total of ", nrow(res_tax_sig), " taxa were identified as significant.", sep="")
}
else if(dim(res_tax_sig)[1]== 0)
  {
  text_metagenomeseq_summary<- paste("No. significant were identified. Adjust the FDR value or change the test correction method")
  #print("No. significant were identified")
}
}
}

  else if(index_pvalue == "pvalue"){
   
      #filter significant
    res_tax_sig = subset(res_tax, PValue < alpha & ((log2FC >= 1) | (log2FC <= -1)))
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
  text_metagenomeseq_summary<- paste("Total of ", nrow(res_tax_sig), " taxa were identified as significant.", sep="")
}
else if(dim(res_tax_sig)[1]== 0)
  {
  text_metagenomeseq_summary<- paste("No. significant were identified. Adjust the Pvalue value or change the test correction method")
  #print("No. significant were identified")
}
}
}
      else
      {
        print("error")
      }
  }
}
return(list(metagenomeseq_result_summary=text_metagenomeseq_summary, metagenomeseq_result_significant=res_tax_sig, metagenomeseq_result_table=res_tax, total_counts=parent_seq_count))
}
}


