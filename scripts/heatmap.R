library("ComplexHeatmap")
library("scales")
library("ggplotify")
heatmap_abundance <- function(OTU_input, group_index, heatmap_row_names, heatmap_row_names_size, heatmap_column_names, heatmap_column_names_size, heatmap_row_dend, heatmap_column_dend){
heatmap_row_dend <- as.logical(heatmap_row_dend)
heatmap_column_dend <- as.logical(heatmap_column_dend)
OTU_input <- OTU_input[rowSums(OTU_input[])>0,] 
df1 <- scale(OTU_input)
ha = HeatmapAnnotation(Group = group_index$Condition, annotation_height = unit(4, "mm"))
p<-Heatmap(df1, name = "Scale", show_row_names = heatmap_row_names, row_names_gp = gpar(fontsize = heatmap_row_names_size), show_column_names = heatmap_column_names, column_names_gp = gpar(fontsize = heatmap_column_names_size), show_row_dend = heatmap_row_dend, show_column_dend = heatmap_column_dend, top_annotation = ha)
p
q <- as.ggplot(p)
return (list(plot=p, plot1=q))
}