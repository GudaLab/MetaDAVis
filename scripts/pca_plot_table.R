library(ggfortify)
pca_plot_table <- function(OTU_input, group_index, pca_label1, pca_label_size, pca_frame){
OTU_input <- OTU_input[rowSums(OTU_input[])>0,]
  pcDat <- prcomp((t(OTU_input)))
  theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
  p <- autoplot(pcDat, data = group_index, colour="Condition", Size =20)+theme # plot PCA
  p
  pca_data <- as.data.frame(p$data[,c(1:3)])
  pca_data <- pca_data %>% mutate(across(where(is.numeric), ~ round(.,2)))
  return(list(PCA_data = pca_data, PCA_data = pcDat))
}
