library("phyloseq")
library("vegan")
library("reshape2")
library("dplyr")
beta_diversity <- function(OTU_input, group_index, tax_index, index_method = "bray", plot_method = "PcoA"){
  #index_method <- readline(prompt="Enter index_method: ")
  #plot_method <- readline(prompt="Enter plot_method: ")
  
  colnames(tax_index)[1]  <- "Species"
  rownames(tax_index) <- tax_index[, 1]
  taxmat =  as.matrix(tax_index)
  OTU = otu_table(OTU_input, taxa_are_rows = TRUE)
  TAX = tax_table(taxmat)
  sampledata = sample_data(group_index)
  physeq1 = phyloseq(OTU, TAX, sampledata)
  sample_data(physeq1)$Group <- factor((sample_data(physeq1)$Condition), levels=unique(sampledata$Condition))
  # we need to transform the data to relative abundance and create a new phyloseq object:
  relab_genera = transform_sample_counts(physeq1, function(x) x / sum(x) * 100) 
  #  plot_method <- readline(prompt="Enter plot_method: ")
  
  
  #summarize_phyloseq(physeq1)
  #print_ps(physeq1)
  #summary(sample_sums(physeq1))
  {
    if(index_method == "bray"){
      diveristy_distance <- phyloseq::distance(relab_genera, method = "bray", digits=2)
    }
    else if(index_method == "jaccard"){
      diveristy_distance <- phyloseq::distance(relab_genera, method = "jaccard")
    }
    else if(index_method == "manhattan"){
      diveristy_distance <- phyloseq::distance(relab_genera, method = "manhattan")
    }
    else if(index_method == "euclidean"){
      diveristy_distance <- phyloseq::distance(relab_genera, method = "euclidean")
    }
    else if(index_method == "canberra"){
      diveristy_distance <- phyloseq::distance(relab_genera, method = "canberra")
    }
    else if(index_method == "kulczynski"){
      diveristy_distance <- phyloseq::distance(relab_genera, method = "kulczynski")
    }
    else if(index_method == "gower"){
      diveristy_distance <- phyloseq::distance(relab_genera, method = "gower")
    }
    else if(index_method == "altGower"){
      diveristy_distance <- phyloseq::distance(relab_genera, method = "altGower")
    }
    else if(index_method == "morisita"){
      diveristy_distance <- phyloseq::distance(relab_genera, method = "morisita")
    }
    else if(index_method == "horn"){
      diveristy_distance <- phyloseq::distance(relab_genera, method = "horn")
    }
    else if(index_method == "mountford"){
      diveristy_distance <- phyloseq::distance(relab_genera, method = "mountford")
    }
    else if(index_method == "raup"){
      diveristy_distance <- phyloseq::distance(relab_genera, method = "raup")
    }
    else if(index_method == "binomial"){
      diveristy_distance <- phyloseq::distance(relab_genera, method = "binomial")
    }
    else if(index_method == "chao"){
      diveristy_distance <- phyloseq::distance(relab_genera, method = "chao")
    }
    else if(index_method == "cao"){
      diveristy_distance <- phyloseq::distance(relab_genera, method = "cao")
    }
	else if(index_method == "w"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "w")
      }
	else if(index_method == "-1"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "-1")
      }
	else if(index_method == "c"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "c")
      }
	else if(index_method == "wb"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "wb")
      }
	else if(index_method == "r"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "r")
      }
	else if(index_method == "I"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "I")
      }
	else if(index_method == "e"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "e")
      }
	else if(index_method == "t"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "t")
      }
	else if(index_method == "me"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "me")
      }
	else if(index_method == "j"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "j")
      }
	else if(index_method == "sor"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "sor")
      }
	else if(index_method == "m"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "m")
      }
	else if(index_method == "-2"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "-2")
      }
	else if(index_method == "co"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "co")
      }
	else if(index_method == "cc"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "cc")
      }
	else if(index_method == "g"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "g")
      }
	else if(index_method == "-3"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "-3")
      }
	else if(index_method == "l"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "l")
      }
	else if(index_method == "19"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "19")
      }
	else if(index_method == "hk"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "hk")
      }
	else if(index_method == "rlb"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "rlb")
      }
	else if(index_method == "sim"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "sim")
      }
	else if(index_method == "gl"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "gl")
      }
	else if(index_method == "z"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "z")
      }
	else if(index_method == "maximum"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "maximum")
      }
	else if(index_method == "binary"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "binary")
      }
	else if(index_method == "minkowski"){
        diveristy_distance <- phyloseq::distance(relab_genera, method = "minkowski")
      }
    else
    {
      print("Select the distance")
    }
    
    diveristy_distance <- as.matrix(diveristy_distance)
    #head(diveristy_distance)[,1:6]
    
    #The “diveristy_distance” matrix presents the distance between all samples, but since we wanna generate a boxplot with distances considering the metagenomes in each group separately, we need to filter this matrix. With the code below we end with a dataframe “df.bray” storing Bray Curtis Distances between metagenomes from the same groups.
    sub_dist <- list()
    groups_all <- sample_data(relab_genera)$Group
    
    for (group in levels(groups_all)) { 
      row_group <- which(groups_all == group)
      sample_group <- sample_names(relab_genera)[row_group]
      sub_dist[[group]] <- diveristy_distance[ sample_group, sample_group]
      sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
    }
    braygroups<- melt(sub_dist)
    df.bray <- braygroups[complete.cases(braygroups), ]
    colnames(df.bray)[4]  <- "Conditions"
    
    df.bray$Conditions <- factor(df.bray$Conditions, levels=names(sub_dist))
    diveristy_distance <- as.data.frame(diveristy_distance)
	diveristy_distance <- diveristy_distance %>% mutate(across(where(is.numeric), ~ round(.,2)))
    #head(df.bray)
    
    if(plot_method == "PCoA"){
      ord = ordinate(relab_genera, method="PCoA", distance = index_method)
    }
    else if(plot_method == "NMDS"){
      ord = ordinate(relab_genera, method="NMDS", distance = index_method)
    }
    else if(plot_method == "DCA"){
      ord = ordinate(relab_genera, method="DCA", distance = index_method)
    }
    else if(plot_method == "CCA"){
      ord = ordinate(relab_genera, method="CCA", distance = index_method)
    }
    else if(plot_method == "RDA"){
      ord = ordinate(relab_genera, method="RDA", distance = index_method)
    }
    else if(plot_method == "MDS"){
      ord = ordinate(relab_genera, method="MDS", distance = index_method)
    }
    else
    {
      print("Select the method")
    }
    #Now we generate the boxplot:
  }
  

  return(list(beta_index=diveristy_distance, beta_index1=df.bray, beta_index2=ord))
}
