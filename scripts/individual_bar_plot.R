Bar_Plot_Individual_OTU <- function(OTU_input, group_index, color_palette_individual, n_top = 15, names="OTU", plot_type){
  OTU_input <- OTU_input[rowSums(OTU_input[])>0,]
  {
    if(plot_type == "1" || plot_type == "2"){
      OTU_input <- OTU_input
    }
    else if (plot_type == "3"){
      #calculate relative frequency
      rel_freq <- (OTU_input)/(colSums(OTU_input))*100
      OTU_input<-as.data.frame(rel_freq)
    }
  }
  {
  group_index1 = group_index[,2]
 if(n_top > nrow(OTU_input)){
    print("Number of selected top OTU number is greater than the total OTU numbers")
  }
  else{
    if(ncol(OTU_input)!=length(group_index1)){
      print("The length of the group index does not match the number of columns in OTU table")
    }
    else{
      Overall_mean <- apply(OTU_input, 1, mean, na.rm = T)
      Top_index <- order(Overall_mean, decreasing = T)[1: n_top]
      Top_OTU_names <- c(rownames(OTU_input)[Top_index], "OTHERS")
      group_tmp <- names(table(group_index1))
      Samples <- c()
      OTU <- c()
      Abundance <- c()
      for(i in 1: length(group_tmp)){
        matrix_tmp <- OTU_input[ , which(group_index1 == group_tmp[i])]
        for(j in 1: ncol(matrix_tmp)){
          top_OTU_abd_tmp <- matrix_tmp[Top_index, j] 
          #others_OTU_abd_tmp <- mean(matrix_tmp[setdiff(c(1: nrow(OTU_input)), Top_index), j], na.rm = T)
          others_OTU_abd_tmp <- 1 - sum(top_OTU_abd_tmp)
          if(others_OTU_abd_tmp < 0){others_OTU_abd_tmp <- 0}
          Samples <- c(Samples, rep(colnames(matrix_tmp)[j], n_top + 1))
          OTU <- c(OTU, Top_OTU_names)
          Abundance <- c(Abundance, c(top_OTU_abd_tmp, others_OTU_abd_tmp))
        }
        #top_means <- apply(matrix_tmp[Top_index, ], 1, mean, na.rm = T)
        #other_mean <- mean(unlist(matrix_tmp[setdiff(c(1: nrow(OTU_input)), Top_index), ]))
        #group <- c(group, rep(group_tmp[i], n_top + 1))
        #OTU <- c(OTU, Top_OTU_names)
        #Abundance <- c(Abundance, c(top_means, other_mean))
      }
      ## data structure for ggplot ##
      #colourCount = length(unique(mtcars$hp))
      #getPalette = colorRampPalette(brewer.pal(9, "Set1"))
      data_gg = data.frame(Samples, OTU, Abundance)
      level_order <- unique(as.character(data_gg$Samples))
      level_OTU <- Top_OTU_names[length(Top_OTU_names):1]
      result<- merge(data_gg, group_index, by = "Samples") 
    }
  }
  }
  colors <- colorRampPalette(brewer.pal(min(12, length(level_OTU)), color_palette_individual))(length(level_OTU))
  {
      
      if(plot_type == "1"){
        ggplot(result, aes(fill = factor(OTU, level = level_OTU), y = Abundance, x = Samples)) + 
          geom_bar( stat="identity", position="fill") + 
          facet_wrap(.~ Condition, scales = "free")+
          theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
          labs(y = "Abundance(%)", fill = names)+
		  scale_fill_manual(values = colors) +
        scale_y_continuous(labels = scales::percent_format())
      }
      else if (plot_type == "2"){
        ggplot(result, aes(fill = factor(OTU, level = level_OTU), y = Abundance, x = Samples)) + 
          geom_bar( stat="identity", position="stack") + 
          facet_wrap(.~ Condition, scales = "free")+
		  scale_fill_manual(values = colors) +
          theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
          labs(y = "Abundance value", fill = names)
      }
    else if (plot_type == "3"){
      ggplot(result, aes(fill = factor(OTU, level = level_OTU), y = Abundance, x = Samples)) +
        geom_bar( stat="identity", position="fill") + 
        facet_wrap(.~ Condition, scales = "free")+
		scale_fill_manual(values = colors) +
        theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
        labs(y = "Relative frequency", fill = names)+
        scale_y_continuous(labels = scales::percent_format())
    }
      else 
      {
        print ("select plot type")
      }
      
    }
  }
