library("ggpubr")
library("GGally")
taxa_based_correlation_plot_table <- function(OTU_input, group_index, method, labe_size){
OTU_input <- as.data.frame(t(OTU_input))
p<-ggcorr(OTU_input, label = FALSE, label_alpha = TRUE, size = labe_size, geom = "circle", method = c("pairwise", method))
p
p$data <- p$data %>% mutate(across(where(is.numeric), ~ round(.,2)))
return(list(plot=p, data=p$data))
}