library("ggpubr")
library("GGally")
samples_based_correlation_plot_table <- function(OTU_input, group_index, method, labe_size){
p<-ggcorr(OTU_input, label = FALSE, label_alpha = TRUE, size = labe_size, geom = "circle", method = c("pairwise", method))
p
p$data <- p$data %>% mutate(across(where(is.numeric), ~ round(.,2)))
return(list(plot=p, data=p$data))
}