data_input_OTU <- function(Input, Index){
  data_tmp <- read.delim(Input, header = T, sep = "\t")
  #metadata file
  data_index <- read.delim(Index, header = T, sep = "\t")
  
  {
    #type <- readline(prompt="Enter type: ")
    type <- as.integer(type)
    
    if(type == "1"){
      #kingdom
      data_kingdom <- data_tmp[,c(1,8: ncol(data_tmp))]
      data_kingdom_sum = aggregate(. ~ Kingdom, data_kingdom, sum)
      data_OTU <- data_kingdom_sum[,2: ncol(data_kingdom_sum)]
      #rownames(data_OTU) <- data_kingdom_sum[, 1]
    }
    else if(type == "2"){
      #Phylum
      data_phylum <- data_tmp[,c(2,8: ncol(data_tmp))]
      data_phylum_sum = aggregate(. ~ Phylum, data_phylum, sum)
      data_OTU <- data_phylum_sum[,2: ncol(data_phylum_sum)]
      #rownames(data_OTU) <- data_phylum_sum[, 1]
    }
    else if(type == "3"){
      #Class
      data_class <- data_tmp[,c(3,8: ncol(data_tmp))]
      data_class_sum = aggregate(. ~ Class, data_class, sum)
      data_OTU <- data_class_sum[,2: ncol(data_class_sum)]
      #rownames(data_OTU) <- data_class_sum[, 1]
    }
    else if(type == "4"){
      #Order
      data_order <- data_tmp[,c(4,8: ncol(data_tmp))]
      data_order_sum = aggregate(. ~ Order, data_order, sum)
      data_OTU <- data_order_sum[,2: ncol(data_order_sum)]
      #rownames(data_OTU) <- data_order_sum[, 1]
    }
    else if(type == "5"){
      #Family
      data_family <- data_tmp[,c(5,8: ncol(data_tmp))]
      data_family_sum = aggregate(. ~ Family, data_family, sum)
      data_OTU <- data_family_sum[,2: ncol(data_family_sum)]
      #rownames(data_OTU) <- data_family_sum[, 1]
    }
    else if(type == "6"){
      #Genus
      data_genus <- data_tmp[,c(6,8: ncol(data_tmp))]
      data_genus_sum = aggregate(. ~ Genus, data_genus, sum)
      data_OTU <- data_genus_sum[,2: ncol(data_genus_sum)]
     #rownames(data_OTU) <- data_genus_sum[, 1]
    }
    else if(type == "7"){
      #Species
      data_species <- data_tmp[,c(7,8: ncol(data_tmp))]
      data_species_sum = aggregate(. ~ Species, data_species, sum)
      data_OTU <- data_species_sum[,2: ncol(data_species_sum)]
      #rownames(data_OTU) <- data_species_sum[, 1]
    }
    else {
      print("Incorrect file")
    }
  }
  
  
  
  #summary
  text_OTU_summary<- paste("There are ", nrow(data_OTU), " bacterial taxa at the", taxa_tmp[as.numeric(type)], " level.", sep="")
  #text_OTU_phylum_summary<- paste("There are ", nrow(data_phylum_sum), " bacterial taxa at the Phylum level.", sep="")
  #text_OTU_class_summary<- paste("There are ", nrow(data_class_sum), " bacterial taxa at the Class level.", sep="")
  #text_OTU_order_summary<- paste("There are ", nrow(data_order_sum), " bacterial taxa at the Order level.", sep="")
  #text_OTU_family_summary<- paste("There are ", nrow(data_family_sum), " bacterial taxa at the Family level.", sep="")
  #text_OTU_genus_summary<- paste("There are ", nrow(data_genus_sum), " bacterial taxa at the Genus level.", sep="")
  #text_OTU_species_summary<- paste("There are ", nrow(data_species_sum), " bacterial taxa at the Species level.", sep="")
  
  #metadata summary
  text_metadata_summary <- paste("Number of ", names(table(data_index[,2])), ": ", table(data_index[,2]), collapse = ", ", sep ="" )
  number_samples <- nrow(data_index)
  
  #taxa <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  #taxa_type <-c(1,2,3,4,5,6,7)
  #type=data.frame(taxa, taxa_type)
  
  #return(data_index)
  return(list(text_OTU_summary = text_OTU_summary, metadata_summary = text_metadata_summary, Number_of_samples = number_samples, Data_OTU = data_OTU, Data_Index = data_index[,2]))
  
}
