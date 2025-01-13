library("tidyr")
library("dplyr")
data_input_RA <- function(file_type, Input, Index, type, sep, sep1){

file_type <- as.character(file_type)
#type <- as.integer(type)
#type <- readline(prompt="Enter type: ")
taxa <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#metadata file
if (file_type == "qiime_format" | file_type == "Megan" | file_type == "check" ){
data_index <- read.delim(Index, header = T, sep = sep1)
colnames(data_index)[1:2] = c("Samples","Condition")
data_index1 <- read.delim(Index, header = T, row.names=1, sep = sep1)
colnames(data_index1) = "Condition"
text_data_index<- paste("There are ", nrow(data_index), " Samples.", sep="")
data_index2 <- unique(data_index[c("Condition")])
rownames(data_index2) <- NULL
}
else if (file_type == "example"){
data_index <- read.delim(Index, header = T, sep = "\t")
colnames(data_index)[1:2] = c("Samples","Condition")
data_index1 <- read.delim(Index, header = T, row.names=1, sep = "\t")
colnames(data_index1) = "Condition"
text_data_index<- paste("There are ", nrow(data_index), " Samples.", sep="")
data_index2 <- unique(data_index[c("Condition")])
rownames(data_index2) <- NULL
}
else {
    stop("Invalid file type. Please check your input.")
  }


{
if (file_type == "qiime_format")
{
data_tmp <- as.data.frame(t(read.delim(Input, header = F, sep = sep)))
#covert 1st row to header
names(data_tmp ) <- as.character(data_tmp[1,])
data_tmp  <- data_tmp [-1,]

colnames(data_tmp)[1] = "index"
#data_tmp$index<-gsub("[kpcofgs]__","", data_tmp$index) #remove leading characters from GG
data_tmp$index<-gsub("[;]__","", data_tmp$index)
data_tmp<-suppressWarnings(data_tmp %>% separate(index, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep="; |;"))
rownames(data_tmp)<-NULL
}
else if(file_type == "Megan"){
data_tmp <- read.delim(Input, header = T, sep = sep)
colnames(data_tmp)[1:7] = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
}
else if(file_type == "example"){
#Input <- read.delim("www/example_data/Megan_WGS_output.tsv", header = T, sep = "\t")
data_tmp <- read.delim(Input, header = T, sep = "\t")
colnames(data_tmp)[1:7] = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
}
else if(file_type == "check"){
data_tmp <- read.delim(Input, header = T, sep = sep)
colnames(data_tmp)[1:7] = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
}
else{
print("check your file format and is it tab  or comma or space separated. Or check your count file sample names and metadata sample names are same.")
}

data_tmp[is.na(data_tmp)] <- ''
data_tmp <- data_tmp %>%  mutate(Phylum = ifelse(Phylum == 'p__',Kingdom,Phylum))
data_tmp <- data_tmp %>%  mutate(Phylum = ifelse(Phylum == '',Kingdom,Phylum))
data_tmp <- data_tmp %>%  mutate(Class = ifelse(Class == 'c__',Phylum,Class))
data_tmp <- data_tmp %>%  mutate(Class = ifelse(Class == '',Phylum,Class))
data_tmp <- data_tmp %>%  mutate(Order = ifelse(Order == 'o__',Class,Order))
data_tmp <- data_tmp %>%  mutate(Order = ifelse(Order == '',Class,Order))
data_tmp <- data_tmp %>%  mutate(Family = ifelse(Family == 'f__',Order,Family))
data_tmp <- data_tmp %>%  mutate(Family = ifelse(Family == '',Order,Family))
data_tmp <- data_tmp %>%  mutate(Genus = ifelse(Genus == 'g__',Family,Genus))
data_tmp <- data_tmp %>%  mutate(Genus = ifelse(Genus == '',Family,Genus))
data_tmp <- data_tmp %>%  mutate(Species = ifelse(Species == 's__',Genus,Species))
data_tmp <- data_tmp %>%  mutate(Species = ifelse(Species == '',Genus,Species))
}

{
   
  if(type == "1"){
    #kingdom
    data_kingdom <- data_tmp[,c(1,8: ncol(data_tmp))]
	data_kingdom[, c(2: ncol(data_kingdom))] <- sapply(data_kingdom[,  c(2: ncol(data_kingdom))], as.numeric)
    data_kingdom_sum = aggregate(. ~ Kingdom, data_kingdom, sum)
	data_kingdom_sum <- data_kingdom_sum[rowSums(data_kingdom_sum[,  c(2: ncol(data_kingdom_sum))])>0,]
    data_OTU <- data_kingdom_sum[,2: ncol(data_kingdom_sum)]
    rownames(data_OTU) <- data_kingdom_sum[, 1]
	Tax_data <- as.data.frame(data_kingdom_sum[ ,1])
  }
  else if(type == "2"){
    #Phylum
    data_phylum <- data_tmp[,c(2,8: ncol(data_tmp))]
	data_phylum[, c(2: ncol(data_phylum))] <- sapply(data_phylum[,  c(2: ncol(data_phylum))], as.numeric)
    data_phylum_sum = aggregate(. ~ Phylum, data_phylum, sum)
	data_phylum_sum <- data_phylum_sum[rowSums(data_phylum_sum[,  c(2: ncol(data_phylum_sum))])>0,]
    data_OTU <- data_phylum_sum[,2: ncol(data_phylum_sum)]
    rownames(data_OTU) <- data_phylum_sum[, 1]
	Tax_data <- as.data.frame(data_phylum_sum[ ,1])
  }
  else if(type == "3"){
    #Class
    data_class <- data_tmp[,c(3,8: ncol(data_tmp))]
	data_class[, c(2: ncol(data_class))] <- sapply(data_class[,  c(2: ncol(data_class))], as.numeric)
    data_class_sum = aggregate(. ~ Class, data_class, sum)
	data_class_sum <- data_class_sum[rowSums(data_class_sum[,  c(2: ncol(data_class_sum))])>0,]
    data_OTU <- data_class_sum[,2: ncol(data_class_sum)]
    rownames(data_OTU) <- data_class_sum[, 1]
	Tax_data <- as.data.frame(data_class_sum[ ,1])
  }
  else if(type == "4"){
    #Order
    data_order <- data_tmp[,c(4,8: ncol(data_tmp))]
	data_order[, c(2: ncol(data_order))] <- sapply(data_order[,  c(2: ncol(data_order))], as.numeric)
    data_order_sum = aggregate(. ~ Order, data_order, sum)
	data_order_sum <- data_order_sum[rowSums(data_order_sum[,  c(2: ncol(data_order_sum))])>0,]
    data_OTU <- data_order_sum[,2: ncol(data_order_sum)]
    rownames(data_OTU) <- data_order_sum[, 1]
	Tax_data <- as.data.frame(data_order_sum[ ,1])
  }
  else if(type == "5"){
    #Family
    data_family <- data_tmp[,c(5,8: ncol(data_tmp))]
	data_family[, c(2: ncol(data_family))] <- sapply(data_family[,  c(2: ncol(data_family))], as.numeric)
    data_family_sum = aggregate(. ~ Family, data_family, sum)
	data_family_sum <- data_family_sum[rowSums(data_family_sum[,  c(2: ncol(data_family_sum))])>0,]
    data_OTU <- data_family_sum[,2: ncol(data_family_sum)]
    rownames(data_OTU) <- data_family_sum[, 1]
	Tax_data <- as.data.frame(data_family_sum[ ,1])
  }
  else if(type == "6"){
    #Genus
    data_genus <- data_tmp[,c(6,8: ncol(data_tmp))]
	data_genus[, c(2: ncol(data_genus))] <- sapply(data_genus[,  c(2: ncol(data_genus))], as.numeric)
    data_genus_sum = aggregate(. ~ Genus, data_genus, sum)
	data_genus_sum <- data_genus_sum[rowSums(data_genus_sum[,  c(2: ncol(data_genus_sum))])>0,]
    data_OTU <- data_genus_sum[,2: ncol(data_genus_sum)]
    rownames(data_OTU) <- data_genus_sum[, 1]
	Tax_data <- as.data.frame(data_genus_sum[ ,1])
  }
  else if(type == "7"){
    #Species
    data_species <- data_tmp[,c(7,8: ncol(data_tmp))]
	data_species[, c(2: ncol(data_species))] <- sapply(data_species[,  c(2: ncol(data_species))], as.numeric)
    data_species_sum = aggregate(. ~ Species, data_species, sum)
	data_species_sum <- data_species_sum[rowSums(data_species_sum[,  c(2: ncol(data_species_sum))])>0,]
    data_OTU <- data_species_sum[,2: ncol(data_species_sum)]
    rownames(data_OTU) <- data_species_sum[, 1]
	Tax_data <- as.data.frame(data_species_sum[ ,1])
	}
   else {
    print("Incorrect file")
  }
}

parent_seq_count <- as.data.frame(colSums(data_OTU))
parent_seq_count <- tibble::rownames_to_column(parent_seq_count, "Samples")
colnames(parent_seq_count)[2] <- "Total_counts"

  
  #summary
  text_OTU_summary<- paste("There are ", nrow(data_OTU), " bacterial taxa at the ", taxa[as.numeric(type)], " level.", sep="")
  #text_OTU_phylum_summary<- paste("There are ", nrow(data_phylum_sum), " bacterial taxa at the Phylum level.", sep="")
  #text_OTU_class_summary<- paste("There are ", nrow(data_class_sum), " bacterial taxa at the Class level.", sep="")
  #text_OTU_order_summary<- paste("There are ", nrow(data_order_sum), " bacterial taxa at the Order level.", sep="")
  #text_OTU_family_summary<- paste("There are ", nrow(data_family_sum), " bacterial taxa at the Family level.", sep="")
  #text_OTU_genus_summary<- paste("There are ", nrow(data_genus_sum), " bacterial taxa at the Genus level.", sep="")
  #text_OTU_species_summary<- paste("There are ", nrow(data_species_sum), " bacterial taxa at the Species level.", sep="")
  
  #metadata summary
  text_metadata_summary <- paste("Number of ", names(table(data_index[,2])), ": ", table(data_index[,2]), collapse = ", ", sep ="" )
  
  number_samples <- nrow(data_index)

 
  #return(data_index)
  return(list(text_OTU_summary = text_OTU_summary, metadata_summary = text_metadata_summary, Number_of_samples = number_samples, Data_OTU = data_OTU, Data_Index = data_index[,2], Type = taxa[as.numeric(type)], Data_Index1 = data_index1, Tax_data=Tax_data, Data_Index2 = data_index, Data_Index3 = data_index2, total_counts = parent_seq_count))
  
}
