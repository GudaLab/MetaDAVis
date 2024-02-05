library("phyloseq")
library("dplyr")
alpha_diversity <- function(OTU_input, group_index, tax_index, index_method = "1", plot_method = "1", pvalue){
  colnames(tax_index)[1]  <- "Species"
  rownames(tax_index) <- tax_index[, 1]
  taxmat =  as.matrix(tax_index)
  OTU = otu_table(OTU_input, taxa_are_rows = TRUE)
  TAX = tax_table(taxmat)
  sampledata = sample_data(group_index)
  physeq1 = phyloseq(OTU, TAX, sampledata)
  #summarize_phyloseq(physeq1)
  #print_ps(physeq1)
  #summary(sample_sums(physeq1))
  {
  if(index_method == "1"){
      alpha_index <- estimate_richness(physeq1, measures="Observed")
	  alpha_index <- alpha_index %>% mutate(across(where(is.numeric), ~ round(.,2)))
	#plot_richness(physeq1, x="Condition", measures="Observed", color="Condition") + geom_boxplot() + theme_bw()
  }
  else if(index_method == "2"){
    alpha_index <- estimate_richness(physeq1, measures="Chao1")
	alpha_index <- alpha_index %>% mutate(across(where(is.numeric), ~ round(.,2)))
    #plot_richness(physeq1, x="Condition", measures="Chao1", color="Condition") + geom_boxplot() + theme_bw()
  }
  else if(index_method == "3"){
    alpha_index <- estimate_richness(physeq1, measures="ACE")
	alpha_index <- alpha_index %>% mutate(across(where(is.numeric), ~ round(.,2)))
	#alpha_index <-formattable(alpha_index,format="f",digits=2)
    #plot_richness(physeq1, x="Condition", measures="ACE", color="Condition") + geom_boxplot() + theme_bw()
  }
  else if(index_method == "4"){
    alpha_index <- estimate_richness(physeq1, measures="Shannon")
	alpha_index <- alpha_index %>% mutate(across(where(is.numeric), ~ round(.,2)))
   # plot_richness(physeq1, x="Condition", measures="Shannon", color="Condition") + geom_boxplot() + theme_bw()
  }
  else if(index_method == "5"){
    alpha_index <- estimate_richness(physeq1, measures="Simpson")
	alpha_index <- alpha_index %>% mutate(across(where(is.numeric), ~ round(.,2)))
    #plot_richness(physeq1, x="Condition", measures="Simpson", color="Condition") + geom_boxplot() + theme_bw()
  }
  else if(index_method == "6"){
    alpha_index <- estimate_richness(physeq1, measures="InvSimpson")
	alpha_index <- alpha_index %>% mutate(across(where(is.numeric), ~ round(.,2)))
    #plot_richness(physeq1, x="Condition", measures="InvSimpson", color="Condition") + geom_boxplot() + theme_bw()
  }
  else if(index_method == "7"){
    alpha_index <- estimate_richness(physeq1, measures="Fisher")
	alpha_index <- alpha_index %>% mutate(across(where(is.numeric), ~ round(.,2)))
    #plot_richness(physeq1, x="Condition", measures="Fisher", color="Condition") + geom_boxplot() + theme_bw()
  }
    else if(index_method == "8"){
      alpha_index <- estimate_richness(physeq1, measures=c("Observed","Chao1","ACE","Shannon", "simpson", "InvSimpson","Fisher"))
	  alpha_index <- alpha_index %>% mutate(across(where(is.numeric), ~ round(.,2)))
      #plot_richness(physeq1, x="Condition", measures=c("Observed","Chao1","ACE","Shannon", "simpson", "InvSimpson","Fisher"), color="Condition") + geom_boxplot() + theme_bw()
    }
 else
 {
   print("error")
 }
  }
  return(list(alpha_index=alpha_index, taxmat=taxmat ))
}