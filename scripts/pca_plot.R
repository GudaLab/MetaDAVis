library("ggfortify")
pca_plot <- function(OTU_input, group_index, pca_label1, pca_label_size, pca_frame){
OTU_input <- OTU_input[rowSums(OTU_input[])>0,]
  pcDat <- prcomp((t(OTU_input)))
  theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
if(pca_label1==TRUE)
{
autoplot(pcDat, data = group_index, colour="Condition", Size = 20, label = TRUE, label.size = pca_label_size, frame = pca_frame, frame.type = 'norm')+theme # plot PCA
}
else {
autoplot(pcDat, data = group_index, colour="Condition", Size = 20, label = FALSE, label.size = pca_label_size, frame = pca_frame, frame.type = 'norm')+theme # plot PCA
}
  #pca_data <- P$data[,c(1:3)]
  #return(list(plot=p, pca_data=pca_data))
}
