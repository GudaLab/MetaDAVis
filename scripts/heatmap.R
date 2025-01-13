library("ComplexHeatmap")
library("scales")
library("ggplotify")
library("circlize")
set.seed(123)
heatmap_abundance <- function(OTU_input, group_index, heatmap_clustering_method_rows, heatmap_clustering_method_columns, heatmap_color_palette, heatmap_normalization, heatmap_row_names, heatmap_row_names_size, heatmap_column_names, heatmap_column_names_size, heatmap_row_dend, heatmap_column_dend){
heatmap_row_dend <- as.logical(heatmap_row_dend)
heatmap_column_dend <- as.logical(heatmap_column_dend)
OTU_input <- OTU_input[rowSums(OTU_input[])>0,] 
#df1 <- scale(OTU_input)
#normalization_method = "none"
if (heatmap_normalization == "scale") {
  df1 <- scale(OTU_input)  # Z-score normalization
} 
else if (heatmap_normalization == "minmax") {
  df1 <- apply(OTU_input, 2, function(x) (x - min(x)) / (max(x) - min(x)))  # Min-max normalization
} 
else if (heatmap_normalization == "log") {
  df1 <- log(OTU_input + 1)  # Log transformation
} 
else if (heatmap_normalization == "row") {
  df1 <- t(apply(OTU_input, 1, function(x) x / sum(x)))  # Row normalization
} 
else if (heatmap_normalization == "column") {
  df1 <- apply(OTU_input, 2, function(x) x / sum(x))  # Column normalization
} 
else {
  df1 <- OTU_input  # No normalization
}
#color_palette = "RdYlBu"
# Create color scheme using the specified palette
  color_scheme <- colorRamp2(
    c(-4, 0, 4), # Breakpoints for scaling
    brewer.pal(3, heatmap_color_palette) # Colors from the chosen palette
  )
ha = HeatmapAnnotation(Group = group_index$Condition, annotation_height = unit(4, "mm"),col = list(Group = structure(brewer.pal(length(unique(group_index$Condition)), heatmap_color_palette), 
                                                     names = unique(group_index$Condition))))

p<-Heatmap(df1, name = "Scale", 
show_row_names = heatmap_row_names, 
row_names_gp = gpar(fontsize = heatmap_row_names_size), 
show_column_names = heatmap_column_names, 
column_names_gp = gpar(fontsize = heatmap_column_names_size), 
show_row_dend = heatmap_row_dend, 
show_column_dend = heatmap_column_dend, 
top_annotation = ha,
col = color_scheme,
clustering_method_columns = heatmap_clustering_method_columns,
clustering_method_rows = heatmap_clustering_method_rows)
p
q <- as.ggplot(p)
return (list(plot=p, plot1=q))
}