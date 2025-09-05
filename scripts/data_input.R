library("tidyr")
library("dplyr")

data_input_RA <- function(file_type, Input, Index, type, sep, sep1, show_head = TRUE, head_n = 5){
  
  log_head <- function(x, label){
    if (show_head) {
      msg <- paste0(
        "\n--- ", label, " (head ", head_n, ") ---\n",
        paste(utils::capture.output(utils::head(x, head_n)), collapse = "\n"),
        "\n"
      )
      cat(msg, file = stderr())
      flush.console()
    }
  }

  
  file_type <- as.character(file_type)
  taxa <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  # metadata file
  if (file_type == "qiime_format" | file_type == "Megan" | file_type == "check"){
    data_index  <- read.delim(Index, header = TRUE, sep = sep1)
    colnames(data_index)[1:2] <- c("Samples","Condition")
    log_head(data_index, "data_index (metadata)")
    
    data_index1 <- read.delim(Index, header = TRUE, row.names = 1, sep = sep1)
    colnames(data_index1) <- "Condition"
    log_head(data_index1, "data_index1 (metadata, row.names=Samples)")
    
    text_data_index <- paste("There are ", nrow(data_index), " Samples.", sep = "")
    data_index2 <- unique(data_index[c("Condition")])
    rownames(data_index2) <- NULL
    log_head(data_index2, "data_index2 (unique conditions)")
  }
  else if (file_type == "example"){
    data_index  <- read.delim(Index, header = TRUE, sep = "\t")
    colnames(data_index)[1:2] <- c("Samples","Condition")
    log_head(data_index, "data_index (metadata)")
    
    data_index1 <- read.delim(Index, header = TRUE, row.names = 1, sep = "\t")
    colnames(data_index1) <- "Condition"
    log_head(data_index1, "data_index1 (metadata, row.names=Samples)")
    
    text_data_index <- paste("There are ", nrow(data_index), " Samples.", sep = "")
    data_index2 <- unique(data_index[c("Condition")])
    rownames(data_index2) <- NULL
    log_head(data_index2, "data_index2 (unique conditions)")
  }
  else {
    stop("Invalid file type. Please check your input.")
  }
  
  # counts / taxonomy input
  if (file_type == "qiime_format"){
    data_tmp <- as.data.frame(t(read.delim(Input, header = FALSE, sep = sep)))
    names(data_tmp) <- as.character(data_tmp[1,])
    data_tmp <- data_tmp[-1,]
    colnames(data_tmp)[1] <- "index"
    data_tmp$index <- gsub("[;]__","", data_tmp$index)
    data_tmp <- suppressWarnings(
      data_tmp %>% separate(index, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = "; |;")
    )
    rownames(data_tmp) <- NULL
    log_head(data_tmp, "data_tmp (QIIME format after split)")
  }
  else if (file_type == "Megan"){
    data_tmp <- read.delim(Input, header = TRUE, sep = sep)
    colnames(data_tmp)[1:7] <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    log_head(data_tmp, "data_tmp (MEGAN format)")
  }
  else if (file_type == "example"){
    data_tmp <- read.delim(Input, header = TRUE, sep = "\t")
    colnames(data_tmp)[1:7] <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    log_head(data_tmp, "data_tmp (example)")
  }
  else if (file_type == "check"){
    data_tmp <- read.delim(Input, header = TRUE, sep = sep)
    colnames(data_tmp)[1:7] <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    log_head(data_tmp, "data_tmp (check)")
  }
  else{
    print("Check your file format and separators. Also ensure sample names match between count and metadata.")
  }
  
  data_tmp[is.na(data_tmp)] <- ''
  data_tmp <- data_tmp %>%  mutate(Phylum = ifelse(Phylum == 'p__' | Phylum == '', Kingdom, Phylum))
  data_tmp <- data_tmp %>%  mutate(Class  = ifelse(Class  == 'c__' | Class  == '', Phylum,  Class))
  data_tmp <- data_tmp %>%  mutate(Order  = ifelse(Order  == 'o__' | Order  == '', Class,   Order))
  data_tmp <- data_tmp %>%  mutate(Family = ifelse(Family == 'f__' | Family == '', Order,   Family))
  data_tmp <- data_tmp %>%  mutate(Genus  = ifelse(Genus  == 'g__' | Genus  == '', Family,  Genus))
  data_tmp <- data_tmp %>%  mutate(Species= ifelse(Species== 's__' | Species== '', Genus,   Species))
  log_head(data_tmp, "data_tmp (filled lineage)")
  
  # collapse to requested rank
  {
    if(type == "1"){
      data_kingdom <- data_tmp[, c(1, 8:ncol(data_tmp))]
      data_kingdom[, 2:ncol(data_kingdom)] <- sapply(data_kingdom[, 2:ncol(data_kingdom)], as.numeric)
      log_head(data_kingdom, "data_kingdom (pre-aggregate)")
      data_kingdom_sum <- aggregate(. ~ Kingdom, data_kingdom, sum)
      data_kingdom_sum <- data_kingdom_sum[rowSums(data_kingdom_sum[, 2:ncol(data_kingdom_sum)]) > 0, ]
      log_head(data_kingdom_sum, "data_kingdom_sum (post-aggregate)")
      data_OTU <- data_kingdom_sum[, 2:ncol(data_kingdom_sum)]
      rownames(data_OTU) <- data_kingdom_sum[, 1]
      Tax_data <- as.data.frame(data_kingdom_sum[, 1])
    }
    else if(type == "2"){
      data_phylum <- data_tmp[, c(2, 8:ncol(data_tmp))]
      data_phylum[, 2:ncol(data_phylum)] <- sapply(data_phylum[, 2:ncol(data_phylum)], as.numeric)
      log_head(data_phylum, "data_phylum (pre-aggregate)")
      data_phylum_sum <- aggregate(. ~ Phylum, data_phylum, sum)
      data_phylum_sum <- data_phylum_sum[rowSums(data_phylum_sum[, 2:ncol(data_phylum_sum)]) > 0, ]
      log_head(data_phylum_sum, "data_phylum_sum (post-aggregate)")
      data_OTU <- data_phylum_sum[, 2:ncol(data_phylum_sum)]
      rownames(data_OTU) <- data_phylum_sum[, 1]
      Tax_data <- as.data.frame(data_phylum_sum[, 1])
    }
    else if(type == "3"){
      data_class <- data_tmp[, c(3, 8:ncol(data_tmp))]
      data_class[, 2:ncol(data_class)] <- sapply(data_class[, 2:ncol(data_class)], as.numeric)
      log_head(data_class, "data_class (pre-aggregate)")
      data_class_sum <- aggregate(. ~ Class, data_class, sum)
      data_class_sum <- data_class_sum[rowSums(data_class_sum[, 2:ncol(data_class_sum)]) > 0, ]
      log_head(data_class_sum, "data_class_sum (post-aggregate)")
      data_OTU <- data_class_sum[, 2:ncol(data_class_sum)]
      rownames(data_OTU) <- data_class_sum[, 1]
      Tax_data <- as.data.frame(data_class_sum[, 1])
    }
    else if(type == "4"){
      data_order <- data_tmp[, c(4, 8:ncol(data_tmp))]
      data_order[, 2:ncol(data_order)] <- sapply(data_order[, 2:ncol(data_order)], as.numeric)
      log_head(data_order, "data_order (pre-aggregate)")
      data_order_sum <- aggregate(. ~ Order, data_order, sum)
      data_order_sum <- data_order_sum[rowSums(data_order_sum[, 2:ncol(data_order_sum)]) > 0, ]
      log_head(data_order_sum, "data_order_sum (post-aggregate)")
      data_OTU <- data_order_sum[, 2:ncol(data_order_sum)]
      rownames(data_OTU) <- data_order_sum[, 1]
      Tax_data <- as.data.frame(data_order_sum[, 1])
    }
    else if(type == "5"){
      data_family <- data_tmp[, c(5, 8:ncol(data_tmp))]
      data_family[, 2:ncol(data_family)] <- sapply(data_family[, 2:ncol(data_family)], as.numeric)
      log_head(data_family, "data_family (pre-aggregate)")
      data_family_sum <- aggregate(. ~ Family, data_family, sum)
      data_family_sum <- data_family_sum[rowSums(data_family_sum[, 2:ncol(data_family_sum)]) > 0, ]
      log_head(data_family_sum, "data_family_sum (post-aggregate)")
      data_OTU <- data_family_sum[, 2:ncol(data_family_sum)]
      rownames(data_OTU) <- data_family_sum[, 1]
      Tax_data <- as.data.frame(data_family_sum[, 1])
    }
    else if(type == "6"){
      data_genus <- data_tmp[, c(6, 8:ncol(data_tmp))]
      data_genus[, 2:ncol(data_genus)] <- sapply(data_genus[, 2:ncol(data_genus)], as.numeric)
      log_head(data_genus, "data_genus (pre-aggregate)")
      data_genus_sum <- aggregate(. ~ Genus, data_genus, sum)
      data_genus_sum <- data_genus_sum[rowSums(data_genus_sum[, 2:ncol(data_genus_sum)]) > 0, ]
      log_head(data_genus_sum, "data_genus_sum (post-aggregate)")
      data_OTU <- data_genus_sum[, 2:ncol(data_genus_sum)]
      rownames(data_OTU) <- data_genus_sum[, 1]
      Tax_data <- as.data.frame(data_genus_sum[, 1])
    }
    else if(type == "7"){
      data_species <- data_tmp[, c(7, 8:ncol(data_tmp))]
      data_species[, 2:ncol(data_species)] <- sapply(data_species[, 2:ncol(data_species)], as.numeric)
      log_head(data_species, "data_species (pre-aggregate)")
      data_species_sum <- aggregate(. ~ Species, data_species, sum)
      data_species_sum <- data_species_sum[rowSums(data_species_sum[, 2:ncol(data_species_sum)]) > 0, ]
      log_head(data_species_sum, "data_species_sum (post-aggregate)")
      data_OTU <- data_species_sum[, 2:ncol(data_species_sum)]
      rownames(data_OTU) <- data_species_sum[, 1]
      Tax_data <- as.data.frame(data_species_sum[, 1])
    }
    else {
      print("Incorrect file")
    }
  }
  
  log_head(data_OTU, "data_OTU (final counts matrix)")
  
  parent_seq_count <- as.data.frame(colSums(data_OTU))
  parent_seq_count <- tibble::rownames_to_column(parent_seq_count, "Samples")
  colnames(parent_seq_count)[2] <- "Total_counts"
  log_head(parent_seq_count, "parent_seq_count (per-sample totals)")
  
  # summaries
  text_OTU_summary <- paste("There are ", nrow(data_OTU), " bacterial taxa at the ", taxa[as.numeric(type)], " level.", sep = "")
  text_metadata_summary <- paste("Number of ", names(table(data_index[,2])), ": ", table(data_index[,2]), collapse = ", ", sep = "")
  
  number_samples <- nrow(data_index)
  
  return(list(
    text_OTU_summary  = text_OTU_summary,
    metadata_summary  = text_metadata_summary,
    Number_of_samples = number_samples,
    Data_OTU          = data_OTU,
    Data_Index        = data_index[,2],
    Type              = taxa[as.numeric(type)],
    Data_Index1       = data_index1,
    Tax_data          = Tax_data,
    Data_Index2       = data_index,
    Data_Index3       = data_index2,
    total_counts      = parent_seq_count
  ))
}
