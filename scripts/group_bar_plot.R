#library("RColorBrewer")

Bar_Plot_Group_OTU <- function(OTU_input, group_index, n_top = 15, names="OTU", plot_type){
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
  if(n_top > nrow(OTU_input)){
      print("Number of selected top OTU number is greater than the total OTU numbers")
    }
    else{
      if(ncol(OTU_input)!=length(group_index)){
        print("The length of the group index does not match the number of columns in OTU table")
      }
      else{
        Overall_mean <- apply(OTU_input, 1, mean, na.rm = T)
        Top_index <- order(Overall_mean, decreasing = T)[1: n_top]
        Top_OTU_names <- c(rownames(OTU_input)[Top_index], "OTHERS")
        group_tmp <- names(table(group_index))
        group <- c()
        OTU <- c()
        Abundance <- c()
        for(i in 1: length(group_tmp)){
          matrix_tmp <- OTU_input[ , which(group_index == group_tmp[i])]
          top_means <- apply(matrix_tmp[Top_index, ], 1, mean, na.rm = T)
          other_mean <- mean(unlist(matrix_tmp[setdiff(c(1: nrow(OTU_input)), Top_index), ]))
          #other_mean <- 1 -sum(top_means)
          group <- c(group, rep(group_tmp[i], n_top + 1))
          OTU <- c(OTU, Top_OTU_names)
          Abundance <- c(Abundance, c(top_means, other_mean))
        }
        ## data structure for ggplot ##
        #colourCount = length(unique(mtcars$hp))
        #getPalette = colorRampPalette(brewer.pal(9, "Set1"))
        data_gg = data.frame(group, OTU, Abundance)
        level_OTU <- Top_OTU_names[length(Top_OTU_names):1]
         
      }
      }
	  }

	  {
      if(plot_type == "1"){
        ggplot(data_gg, aes(fill = factor(OTU, level = level_OTU), y = Abundance, x = group)) + 
          geom_bar(stat="identity" , position="fill") + 
          labs(y = "Abundance(%)", fill = names) +
          scale_y_continuous(labels = scales::percent_format())+
          theme(axis.text=element_text(size=15), axis.title=element_text(size=20), legend.title = element_text(size=15), legend.text = element_text(size=10))+
          theme_bw()
              }
      else if (plot_type == "2"){
          ggplot(data_gg, aes(fill = factor(OTU, level = level_OTU), y = Abundance, x = group)) + 
            geom_bar(stat="identity" , position="stack") + 
            labs(y = "Abundance value", fill = names) +
            #scale_fill_manual(values = getPalette(colourCount)) +
            theme(axis.text=element_text(size=15), axis.title=element_text(size=20), legend.title = element_text(size=15), legend.text = element_text(size=10))+
            theme_bw()
      }
	 else if (plot_type == "3"){
          ggplot(data_gg, aes(fill = factor(OTU, level = level_OTU), y = Abundance, x = group)) + 
            geom_bar(stat="identity" , position="fill") + 
            labs(y = "Relative frequency", fill = names) +
			scale_y_continuous(labels = scales::percent_format())+
			theme(axis.text=element_text(size=15), axis.title=element_text(size=20), legend.title = element_text(size=15), legend.text = element_text(size=10))+
            theme_bw()
      }
      else 
      {
        print ("select plot type")
      }
      }
    }
