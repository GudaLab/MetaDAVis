library("phyloseq")
library("bluster")
library("patchwork")
library("scater")
library("mia")

umap_plot_table <- function(OTU_input, group_index, tax_index, method, kvalue){
    
  colnames(tax_index)[1]  <- "Species"
  rownames(tax_index) <- tax_index[, 1]
  taxmat =  as.matrix(tax_index)
  OTU = otu_table(OTU_input, taxa_are_rows = TRUE)
  TAX = tax_table(taxmat)
  sampledata = sample_data(group_index)
  physeq1 = phyloseq(OTU, TAX, sampledata)

tse <- makeTreeSummarizedExperimentFromPhyloseq(physeq1) 
tse
kvalue=as.numeric(kvalue)
{
 if(method == "counts"){
    tse <- runUMAP(tse, name="UMAP", exprs_values="counts")
	graph_clusters <- clusterRows(t(assays(tse)$counts), NNGraphParam(k=kvalue))
	plots<-plotUMAP(tse, colour_by = "Condition") + labs(title = paste0("counts"))
	plots1<-plotUMAP(tse, colour_by = I(graph_clusters)) + labs(title = paste0("k = ", kvalue))
	umapplot <- (plots + plots1)
	umapplot
	}
 else if(method == "rclr"){
    tse <- transformCounts(tse, method = "rclr")
	tse <- runUMAP(tse, name="UMAP", exprs_values="rclr")
	graph_clusters <- clusterRows(t(assays(tse)$rclr), NNGraphParam(k=kvalue))
	plots<-plotUMAP(tse, colour_by = "Condition") + labs(title = paste0("rclr"))
	plots1<-plotUMAP(tse, colour_by = I(graph_clusters)) + labs(title = paste0("k = ", kvalue))
	umapplot <- (plots + plots1)
	umapplot
	}
 else if(method == "hellinger"){
    tse <- transformCounts(tse, method = "hellinger")
	tse <- runUMAP(tse, name="UMAP", exprs_values="hellinger")
	graph_clusters <- clusterRows(t(assays(tse)$hellinger), NNGraphParam(k=kvalue))
	plots<-plotUMAP(tse, colour_by = "Condition") + labs(title = paste0("hellinger"))
	plots1<-plotUMAP(tse, colour_by = I(graph_clusters)) + labs(title = paste0("k = ", kvalue))
	umapplot <- (plots + plots1)
	umapplot
	}
 else if(method == "pa"){
    tse <- transformCounts(tse, method = "pa")
	tse <- runUMAP(tse, name="UMAP", exprs_values="pa")
	graph_clusters <- clusterRows(t(assays(tse)$pa), NNGraphParam(k=kvalue))
	plots<-plotUMAP(tse, colour_by = "Condition") + labs(title = paste0("pa"))
	plots1<-plotUMAP(tse, colour_by = I(graph_clusters)) + labs(title = paste0("k = ", kvalue))
	umapplot <- (plots + plots1)
	umapplot
	}
 else if(method == "rank"){
    tse <- transformCounts(tse, method = "rank")
	tse <- runUMAP(tse, name="UMAP", exprs_values="rank")
	graph_clusters <- clusterRows(t(assays(tse)$rank), NNGraphParam(k=kvalue))
	plots<-plotUMAP(tse, colour_by = "Condition") + labs(title = paste0("rank"))
	plots1<-plotUMAP(tse, colour_by = I(graph_clusters)) + labs(title = paste0("k = ", kvalue))
	umapplot <- (plots + plots1)
	umapplot
	}
 else if(method == "relabundance"){
    tse <- transformCounts(tse, method = "relabundance")
	tse <- runUMAP(tse, name="UMAP", exprs_values="relabundance")
	graph_clusters <- clusterRows(t(assays(tse)$relabundance), NNGraphParam(k=kvalue))
	plots<-plotUMAP(tse, colour_by = "Condition") + labs(title = paste0("relabundance"))
	plots1<-plotUMAP(tse, colour_by = I(graph_clusters)) + labs(title = paste0("k = ", kvalue))
	umapplot <- (plots + plots1)
	umapplot
	}
 else
    {
    print("Select the method")
    }
	}
	plots$data <- plots$data %>% mutate(across(where(is.numeric), ~ round(.,2)))
	plots1$data <- plots1$data %>% mutate(across(where(is.numeric), ~ round(.,2)))
return(list(plot=umapplot, data=plots$data, data1=plots1$data))
}