library("phyloseq")
library("microbiomeutilities")
alpha_diversity_boxplot <- function(OTU_input, group_index, tax_index, index_method = "1", plot_method = "1", pvalue){
  colnames(tax_index)[1]  <- "Species"
  rownames(tax_index) <- tax_index[, 1]
  taxmat =  as.matrix(tax_index)
  OTU = otu_table(OTU_input, taxa_are_rows = TRUE)
  TAX = tax_table(taxmat)
  sampledata = sample_data(group_index)
  physeq1 = phyloseq(OTU, TAX, sampledata)
  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
  #summarize_phyloseq(physeq1)
  #print_ps(physeq1)
  #summary(sample_sums(physeq1))
  {
    if(index_method == "1"){
      #alpha_index <- estimate_richness(physeq1, measures="Observed")
      p <- plot_richness(physeq1, x="Condition", measures="Observed", color="Condition") + geom_boxplot() + theme_bw()
      q <- plot_richness(physeq1, x="Condition", measures="Observed", color="Condition") + geom_violin() + theme_bw()
    }
    else if(index_method == "2"){
      #alpha_index <- estimate_richness(physeq1, measures="Chao1")
      p <- plot_richness(physeq1, x="Condition", measures="Chao1", color="Condition") + geom_boxplot() + theme_bw()
      q <- plot_richness(physeq1, x="Condition", measures="Chao1", color="Condition") + geom_violin() + theme_bw()
    }
    else if(index_method == "3"){
      #alpha_index <- estimate_richness(physeq1, measures="ACE")
      p <- plot_richness(physeq1, x="Condition", measures="ACE", color="Condition") + geom_boxplot() + theme_bw()
      q <- plot_richness(physeq1, x="Condition", measures="ACE", color="Condition") + geom_violin() + theme_bw()
    }
    else if(index_method == "4"){
      #alpha_index <- estimate_richness(physeq1, measures="Shannon")
      p <- plot_richness(physeq1, x="Condition", measures="Shannon", color="Condition") + geom_boxplot() + theme_bw()
      q <- plot_richness(physeq1, x="Condition", measures="Shannon", color="Condition") + geom_violin() + theme_bw()
    }
    else if(index_method == "5"){
      #alpha_index <- estimate_richness(physeq1, measures="Simpson")
      p <- plot_richness(physeq1, x="Condition", measures="Simpson", color="Condition") + geom_boxplot() + theme_bw()
      q <- plot_richness(physeq1, x="Condition", measures="Simpson", color="Condition") + geom_violin() + theme_bw()
    }
    else if(index_method == "6"){
      #alpha_index <- estimate_richness(physeq1, measures="InvSimpson")
      p <- plot_richness(physeq1, x="Condition", measures="InvSimpson", color="Condition") + geom_boxplot() + theme_bw()
      q <- plot_richness(physeq1, x="Condition", measures="InvSimpson", color="Condition") + geom_violin() + theme_bw()
    }
    else if(index_method == "7"){
      #alpha_index <- estimate_richness(physeq1, measures="Fisher")
      p <- plot_richness(physeq1, x="Condition", measures="Fisher", color="Condition") + geom_boxplot() + theme_bw()
      q <- plot_richness(physeq1, x="Condition", measures="Fisher", color="Condition") + geom_violin() + theme_bw()
    }
    else if(index_method == "8"){
      #alpha_index <- estimate_richness(physeq1, measures=c("Observed","Chao1","ACE","Shannon", "simpson", "InvSimpson","Fisher"))
      p <- plot_richness(physeq1, x="Condition", measures=c("Observed","Chao1","ACE","Shannon", "simpson", "InvSimpson","Fisher"), color="Condition") + geom_boxplot() + theme_bw()
      q <- plot_richness(physeq1, x="Condition", measures=c("Observed","Chao1","ACE","Shannon", "simpson", "InvSimpson","Fisher"), color="Condition") + geom_violin() + theme_bw()
    }
    else
    {
      print("error")
    }
  }
  if(plot_method == "1" & pvalue== "No"){
    p
  }
  else if(plot_method == "1" & pvalue== "Yes"){
      p
	  comps <- make_pairs(sample_data(physeq1)$Condition)
      p+stat_compare_means(comparisons = comps,  label = "p.format",tip.length = 0.01, method = "wilcox.test")
  }
  else if(plot_method == "2" & pvalue== "Yes"){
      q
	  comps <- make_pairs(sample_data(physeq1)$Condition)
      q+stat_compare_means(comparisons = comps,  label = "p.format",tip.length = 0.01, method = "wilcox.test")
  }
    else if(plot_method == "1" & pvalue== "star"){
      p
	  comps <- make_pairs(sample_data(physeq1)$Condition)
      p+stat_compare_means(comparisons = comps,  label = "p.format",tip.length = 0.01, method = "wilcox.test", symnum.args = symnum.args )
  }
  else if(plot_method == "2" & pvalue== "star"){
      q
	  comps <- make_pairs(sample_data(physeq1)$Condition)
      q+stat_compare_means(comparisons = comps,  label = "p.format",tip.length = 0.01, method = "wilcox.test", symnum.args = symnum.args)
  }
    else
  {
    q
  }
}