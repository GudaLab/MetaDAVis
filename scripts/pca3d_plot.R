library("plotly")
pca3d_plot <- function(OTU_input, group_index){
OTU_input <- as.data.frame(t(OTU_input))
prin_comp <- prcomp(OTU_input, rank. = 3)
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3

components = merge(components, group_index, by = 'row.names', all = TRUE )
colnames(components)[1]<- "Samples"
components <- components %>% mutate(across(where(is.numeric), ~ round(.,2)))

fig <- plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3, text= ~Samples, color = ~Condition) %>% add_markers(size = 12)

fig <- fig %>%  layout(
    legend=list(title=list(text='Condition')),
    scene = list(bgcolor = "white")
  )
fig

return(list(plot=fig, data=components, data1=components))
}
