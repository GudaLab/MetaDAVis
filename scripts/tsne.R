library("phyloseq")
library("bluster")
library("patchwork")
library("scater")
library("mia")

tsne_plot_table <- function(OTU_input, group_index, tax_index, method, dimension){
    
  colnames(tax_index)[1]  <- "Species"
  rownames(tax_index) <- tax_index[, 1]
  taxmat =  as.matrix(tax_index)
  OTU = otu_table(OTU_input, taxa_are_rows = TRUE)
  TAX = tax_table(taxmat)
  sampledata = sample_data(group_index)
  physeq1 = phyloseq(OTU, TAX, sampledata)

tse <- makeTreeSummarizedExperimentFromPhyloseq(physeq1) 
tse
dimension=as.numeric(dimension)
{
 if(method == "counts"){
    tse <- runTSNE(tse, name = "TSNE", exprs_values = "counts", ncomponents = dimension)
	plottsne <-plotReducedDim(tse, "TSNE", colour_by = "Condition", ncomponents = c(1:dimension))
    plottsne
	}
 else if(method == "rclr"){
    tse <- transformCounts(tse, method = "rclr")
	tse <- runTSNE(tse, name = "TSNE", exprs_values = "rclr", ncomponents = dimension)
	plottsne <-plotReducedDim(tse, "TSNE", colour_by = "Condition", ncomponents = c(1:dimension))
    plottsne
	}
 else if(method == "hellinger"){
    tse <- transformCounts(tse, method = "hellinger")
	tse <- runTSNE(tse, name = "TSNE", exprs_values = "hellinger", ncomponents = dimension)
	plottsne <-plotReducedDim(tse, "TSNE", colour_by = "Condition", ncomponents = c(1:dimension))
    plottsne
	}
 else if(method == "pa"){
    tse <- transformCounts(tse, method = "pa")
	tse <- runTSNE(tse, name = "TSNE", exprs_values = "pa", ncomponents = dimension)
	plottsne <-plotReducedDim(tse, "TSNE", colour_by = "Condition", ncomponents = c(1:dimension))
    plottsne
	}
 else if(method == "rank"){
    tse <- transformCounts(tse, method = "rank")
	tse <- runTSNE(tse, name = "TSNE", exprs_values = "rank", ncomponents = dimension)
	plottsne <-plotReducedDim(tse, "TSNE", colour_by = "Condition", ncomponents = c(1:dimension))
    plottsne
	}
 else if(method == "relabundance"){
    tse <- transformCounts(tse, method = "relabundance")
	tse <- runTSNE(tse, name = "TSNE", exprs_values = "relabundance", ncomponents = dimension)
	plottsne <-plotReducedDim(tse, "TSNE", colour_by = "Condition", ncomponents = c(1:dimension))
    plottsne
	}
 else
    {
    print("Select the method")
    }
	}
	plottsne$data <- plottsne$data %>% mutate(across(where(is.numeric), ~ round(.,2)))
return(list(plot=plottsne, data=plottsne$data))
}