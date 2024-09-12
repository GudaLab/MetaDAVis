source("global.R")
#library("R.utils")
#options(warn=-1)
server <- function(input, output, session) {
  
  ###########################
  ##  Download sample data ##
  ###########################  
  ## *** Download sample data in csv format ***
  #output$Qiime2_output_level7.csv <- downloadfile("example_data/Qiime2_output_level7.csv")
  #output$Qiime2_metadata.csv <- downloadfile("example_data/Qiime2_metadata.csv")
  #output$Megan_WGS_output.tsv <- downloadfile("Megan_WGS_output.tsv")
  #output$Megan_WGS_metadata.tsv <- R.utils::downloadfile("Megan_WGS_metadata.tsv")
  #output$Taxa_count_file.tsv <- downloadfile("Taxa_count_file.tsv")
  #output$Qiime2_metadata_for_Silva.csv <- download.file("example_data/Qiime2_metadata_for_Silva.csv")

  
  
###########################
##         Input         ##
###########################  
  dataInput_RA_level<- eventReactive(input$action_level,{
    inFile1 <- input$file1
    inFile2 <- input$metadata
    sep=input$file1_split
    sep1=input$metadata_split
    type=input$select_file_type
    ext <- tools::file_ext(input$file1)
    req(input$file1)
    validate(need(((type=="Megan" | type=="check") & (ext == "tsv") & (sep=="\t")) | ((type=="qiime_format") & (ext == "csv") & (sep==",")), "Please select the proper input format and 'Fields separated by' "))
    ext1 <- tools::file_ext(input$metadata)
    req(input$metadata)    
    validate(need(((type=="Megan" | type=="check") & (ext1 == "tsv") & (sep1=="\t")) | ((type=="qiime_format") & (ext1 == "csv") & (sep1==",")), "Please select the proper input format and 'Fields separated by' "))
    
    #validate(need((ext1 == "csv") & (sep1==",") | (ext1 == "tsv") & (sep1=="\t"), "Please upload a .csv or .tsv file and select the proper menu 'Fields separated by' "))
    
    
    if(is.null(input$file1)|is.null(input$metadata)){
      return(list("Input is missing", "Input is missing", 30, "Input is missing"))
    }
    else{
      data_input_RA(Input = inFile1$datapath, Index = inFile2$datapath, type=input$select_RA_type, sep=input$file1_split, sep1=input$metadata_split, file_type=input$select_file_type)
    }  
  })
  
  output$text_level<- renderText({
    paste(dataInput_RA_level()[[1]])
  })
  output$text_metadata<- renderText({
    paste(dataInput_RA_level()[[2]])
  })
  output$taxonomy_table <- renderDataTable(DT::datatable(dataInput_RA_level()[[4]],options = list(pageLength = 15,scrollX = TRUE)))
  output$metadata_table <- renderDataTable(DT::datatable(dataInput_RA_level()[[7]],options = list(pageLength = 15,scrollX = TRUE)))
  output$conditions_table <- renderDataTable(DT::datatable(dataInput_RA_level()[[10]],options = list(pageLength = 15,scrollX = TRUE)))
  output$count_table <- renderDataTable(DT::datatable(dataInput_RA_level()[[11]],options = list(pageLength = 15,scrollX = TRUE)))
  
  output$download_taxonomy_table <- downloadHandler(
    filename = function() { 
      paste("taxonomy_table", '.csv', sep='') },
    content = function(file){
      write.csv(dataInput_RA_level()[[4]], file, row.names = TRUE)
    }
  )
  output$download_metadata_table <- downloadHandler(
    filename = function() { 
      paste("metadata_table", '.csv', sep='') },
    content = function(file){
      write.csv(dataInput_RA_level()[[7]], file, row.names = TRUE)
    }
  )
  output$download_conditions_table <- downloadHandler(
    filename = function() { 
      paste("conditions_table", '.csv', sep='') },
    content = function(file){
      write.csv(dataInput_RA_level()[[10]], file, row.names = TRUE)
    }
  )
  output$download_count_table <- downloadHandler(
    filename = function() { 
      paste("count_table", '.csv', sep='') },
    content = function(file){
      write.csv(dataInput_RA_level()[[11]], file, row.names = FALSE)
    }
  )

  
  
  ###########################
  ##    Bar-Plot: group    ##
  ###########################   

  observe({
    labels_data_type <- input$file1$name
    updateSelectInput(session, "input_RA_bar_plot_group", choices = labels_data_type)
    updateSelectInput(session, "input_RA_bar_plot_individual", choices = labels_data_type)
  })
  
  data_bar_plot_group <- eventReactive(input$action_m1_bar_plot_group,{
    source("scripts/group_bar_plot.R")
    Bar_Plot_Group_OTU(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[5]], n_top = input$top_n_bar_plot_group, names = dataInput_RA_level()[[6]], plot_type=input$select_plot_type_group)
  })
  
  output$bar_plot_group <- renderPlot({
    data_bar_plot_group()
  })
  
  output$download_bar_plot_group<- downloadHandler(
    filename = function(){
      paste("Bar_Plot_Group_Top_", input$top_n_bar_plot_group,"_", dataInput_RA_level()[[6]],input$select_image_type_group, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_bar_plot_group(), width = input$bar_plot_group_output_width, height = input$bar_plot_group_output_height, dpi = input$bar_plot_group_output_dpi, units = "in")
    }
  )
  
  
 
  ###########################
  ##  Bar-Plot: individual ##
  ###########################   

  
  data_bar_plot_individual <- eventReactive(input$action_m1_bar_plot_individual,{
    source("scripts/individual_bar_plot.R")
    Bar_Plot_Individual_OTU(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], n_top = input$top_n_bar_plot_individual, names = dataInput_RA_level()[[6]], plot_type=input$select_plot_type_individual)
  })
  
  output$bar_plot_individual <- renderPlot({
    data_bar_plot_individual()
  })
  
  
  output$download_bar_plot_individual<- downloadHandler(
    filename = function(){
      paste("Bar_Plot_Individual_Top_", input$top_n_bar_plot_individual,"_", dataInput_RA_level()[[6]],input$select_image_type_individual, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_bar_plot_individual(), width = input$bar_plot_individual_output_width, height = input$bar_plot_individual_output_height, dpi = input$bar_plot_individual_output_dpi, units = "in")
    }
  )
  
  
  ###########################
  ##        Heatmap        ##
  ########################### 
  
  #Bar-Plot: group
  observe({
    labels_data_type <- input$file1$name
    updateSelectInput(session, "input_heatmap", choices = labels_data_type)
    })
  
  data_heatmap <- eventReactive(input$action_heatmap,{
    source("scripts/heatmap.R")
    heatmap_abundance(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], heatmap_row_names = input$heatmap_row_names, heatmap_row_names_size = input$heatmap_row_names_size, heatmap_column_names = input$heatmap_column_names, heatmap_column_names_size = input$heatmap_column_names_size, heatmap_row_dend = input$heatmap_row_dend, heatmap_column_dend = input$heatmap_column_dend)
  })
  
  output$plot_heatmap <- renderPlot({
    data_heatmap()[[1]]
  })
  
  output$download_heatmap<- downloadHandler(
    filename = function(){
      paste("Heatmap",input$select_image_type_heatmap, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_heatmap()[[2]], width = input$heatmap_output_width, height = input$heatmap_output_height, dpi = input$heatmap_output_dpi, units = "in")
    }
  )
  
  
  ###########################
  ##    Alpha Diversity    ##
  ###########################   

  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_RA_Alpha", choices = labels_data_type)
  })
  
  data_Alpha_Div <- eventReactive(input$action_alpha_diversity,{
    source("scripts/alpha_diversity.R")
    labels_data_type<- input$file1$name
    alpha_diversity(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[7]], tax_index = dataInput_RA_level()[[8]], index_method = input$select_alpha, plot_method = input$select_plot, pvalue=input$select_alpha_pvalue)
  })
  
  output$alpha_table <- renderDataTable(DT::datatable(data_Alpha_Div()[[1]],options = list(pageLength = 15,scrollX = TRUE)))
  

  output$download_result_alpha <- downloadHandler(
    filename = function() { 
      paste("Alpha_diversity_result", '.csv', sep='') },
    content = function(file){
      write.csv(data_Alpha_Div()[[1]], file)
    }
  )

    #boxplot
  data_Alpha_Div_plot <- eventReactive(input$action_alpha_diversity,{
    source("scripts/alpha_diversity_boxplot.R")
    labels_data_type<- input$file1$name
    alpha_diversity_boxplot(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[7]], tax_index = dataInput_RA_level()[[8]], index_method = input$select_alpha, plot_method = input$select_plot, pvalue=input$select_alpha_pvalue)
  })
  
  output$boxplot_Alpha_Div <- renderPlot({
    data_Alpha_Div_plot()
  })
  
  
  output$download_Boxplot_Alpha_Div<- downloadHandler(
    filename = function(){
      method_tmp <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher", "All_Combined")
      paste("Alpha_Diversity_", method_tmp[as.numeric(input$select_alpha)],"_index", input$select_image_type_alpha, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_Alpha_Div_plot(), width = input$Boxplot_alpha_div_output_width, height = input$Boxplot_alpha_div_output_height, dpi = input$Boxplot_alpha_div_output_dpi, units = "in")
    }
  )
  
  

  
  ###########################
  ##     Beta Diversity    ##
  ########################### 
  

  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_RA_Beta", choices = labels_data_type)
  })
  
  data_beta_Div <- eventReactive(input$action_beta_diversity,{
    source("scripts/beta_diversity.R")
    labels_data_type<- input$file1$name
    beta_diversity(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[7]], tax_index = dataInput_RA_level()[[8]], index_method = input$select_beta, plot_method = input$select_method)
  })

  output$beta_table <- renderDataTable(DT::datatable(data_beta_Div()[[1]],options = list(pageLength = 15,scrollX = TRUE)))
  
  output$download_result_beta <- downloadHandler(
    filename = function() { 
      paste("Beta_diversity_result", '.csv', sep='') },
    content = function(file){
      write.csv(data_beta_Div()[[1]], file)
    }
  )
  
  
  #boxplot
  data_Beta_Div_plot <- eventReactive(input$action_beta_diversity,{
    source("scripts/beta_diversity_boxplot.R")
    labels_data_type<- input$file1$name
    beta_diversity_boxplot(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[7]], tax_index = dataInput_RA_level()[[8]], index_method = input$select_beta, plot_method = input$select_method)
  })
  
  output$boxplot_Beta_Div <- renderPlot({
    data_Beta_Div_plot()
  })
  
  output$download_Boxplot_beta_Div<- downloadHandler(
    filename = function(){
      paste("Beta_Diversity_", input$select_beta, "_",input$select_method, input$select_image_type_beta, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_Beta_Div_plot(), width = input$Boxplot_beta_div_output_width, height = input$Boxplot_beta_div_output_height, dpi = input$Boxplot_beta_div_output_dpi, units = "in")
    }
  )


  
  
  
  ###########################
  ##         PCA-2D        ##
  ###########################  
  

  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_RA_pca", choices = labels_data_type)
  })
  
  data_pca_table <- eventReactive(input$action_pca,{
    source("scripts/pca_plot_table.R")
    labels_data_type<- input$file1$name
    pca_plot_table(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[7]], pca_label1 = input$select_pca_label, pca_label_size = input$select_pca_label_size, pca_frame = input$select_pca_frame)
  })
  
 output$pca_table <- renderDataTable(DT::datatable(data_pca_table()[[1]],options = list(pageLength = 15,scrollX = TRUE)))
 
 output$download_result_pca <- downloadHandler(
   filename = function() { 
     paste("pca_result", '.csv', sep='') },
   content = function(file){
     write.csv(data_pca_table()[[1]], file)
   }
 )
 
 
  #boxplot
  data_pca_plot <- eventReactive(input$action_pca,{
    source("scripts/pca_plot.R")
    labels_data_type<- input$file1$name
    pca_plot(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[7]], pca_label1 = input$select_pca_label, pca_label_size = input$select_pca_label_size, pca_frame = input$select_pca_frame)
  })
  
  output$plot_pca <- renderPlot({
    data_pca_plot()
  })
  
  output$download_plot_pca <- downloadHandler(
    filename = function(){
     paste("PCA_plot", input$select_image_type_pca, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_pca_plot(), width = input$pca_plot_output_width, height = input$pca_plot_output_height, dpi = input$pca_plot_output_dpi, units = "in")
    }
  )  
  
  
    
  ###########################
  ##         PCA-3D        ##
  ###########################  
  
  
  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_RA_pca3d", choices = labels_data_type)
  })
  
  data_pca3d_table <- eventReactive(input$action_pca3d,{
    source("scripts/pca3d_plot.R")
    labels_data_type<- input$file1$name
    pca3d_plot(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[7]])
  })
  
  output$pca3d_table <- renderDataTable(DT::datatable(data_pca3d_table()[[2]],options = list(pageLength = 15,scrollX = TRUE)))
  
  output$download_result_pca3d <- downloadHandler(
    filename = function() { 
      paste("pca3d_result", '.csv', sep='') },
    content = function(file){
      write.csv(data_pca3d_table()[[2]], file)
    }
  )
  
    output$plot_pca3d <- renderPlotly({
      data_pca3d_table()[[1]]
  })
  
  
  
  ###########################
  ##          t-SNE        ##
  ###########################  
  
  
  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_tsne", choices = labels_data_type)
  })
  
  data_tsne_table <- eventReactive(input$action_tsne,{
    source("scripts/tsne.R")
    labels_data_type<- input$file1$name
    tsne_plot_table(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[7]], tax_index = dataInput_RA_level()[[8]], method = input$select_tsne_method, dimension = input$select_tsne_dimension)
  })
  
  output$tsne_table <- renderDataTable(DT::datatable(data_tsne_table()[[2]],options = list(pageLength = 15,scrollX = TRUE)))
  
  output$download_result_tsne <- downloadHandler(
    filename = function() { 
      paste("tsne_result_", input$select_tsne_method, '.csv', sep='') },
    content = function(file){
      write.csv(data_tsne_table()[[2]], file)
    }
  )
  
  output$plot_tsne <- renderPlot({
    data_tsne_table()[[1]]
  })
  
  output$download_plot_tsne <- downloadHandler(
    filename = function(){
      paste("tsne_plot_", input$select_tsne_method, input$select_image_type_tsne, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_tsne_table()[[1]], width = input$tsne_plot_output_width, height = input$tsne_plot_output_height, dpi = input$tsne_plot_output_dpi, units = "in")
    }
  )  
  
  
  
  
  
  ###########################
  ##          UMAP         ##
  ###########################  
  
  
  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_umap", choices = labels_data_type)
  })
  
  data_umap_table <- eventReactive(input$action_umap,{
    source("scripts/umap.R")
    labels_data_type<- input$file1$name
    umap_plot_table(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[7]], tax_index = dataInput_RA_level()[[8]], method = input$select_umap_method, kvalue = input$select_umap_kvalue)
  })
  
  output$umap_table <- renderDataTable(DT::datatable(data_umap_table()[[2]],options = list(pageLength = 15,scrollX = TRUE)))
  
  output$download_result_umap <- downloadHandler(
    filename = function() { 
      paste("umap_result_based_on_condition_", input$select_umap_method, '.csv', sep='') },
    content = function(file){
      write.csv(data_umap_table()[[2]], file)
    }
  )
  
  output$umap_table1 <- renderDataTable(DT::datatable(data_umap_table()[[3]],options = list(pageLength = 15,scrollX = TRUE)))
  
  output$download_result_umap1 <- downloadHandler(
    filename = function() { 
      paste("umap_resul_based_on_cluster_", input$select_umap_method, '.csv', sep='') },
    content = function(file){
      write.csv(data_umap_table()[[3]], file)
    }
  )
  
  output$plot_umap <- renderPlot({
    data_umap_table()[[1]]
  })
  
  output$download_plot_umap <- downloadHandler(
    filename = function(){
      paste("umap_plot_", input$select_umap_method, input$select_image_type_umap, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_umap_table()[[1]], width = input$umap_plot_output_width, height = input$umap_plot_output_height, dpi = input$umap_plot_output_dpi, units = "in")
    }
  )  
  
  
  #############################
  ## Correlation: Taxa-based ##
  #############################  
  
  
  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_taxa_based_correlation", choices = labels_data_type)
  })
  
  data_taxa_based_correlation_table <- eventReactive(input$action_taxa_based_correlation,{
    source("scripts/taxa_based_correlation.R")
    labels_data_type<- input$file1$name
    taxa_based_correlation_plot_table(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[7]], method = input$select_taxa_based_correlation_method, labe_size =input$select_taxa_based_correlation_label_size)
  })
  
  output$taxa_based_correlation_table <- renderDataTable(DT::datatable(data_taxa_based_correlation_table()[[2]],options = list(pageLength = 15,scrollX = TRUE)))
  
  output$download_result_taxa_based_correlation <- downloadHandler(
    filename = function() { 
      paste("taxa_based_correlation_result_", input$select_taxa_based_correlation_method, '.csv', sep='') },
    content = function(file){
      write.csv(data_taxa_based_correlation_table()[[2]], file)
    }
  )
  
  output$plot_taxa_based_correlation <- renderPlot({
    data_taxa_based_correlation_table()[[1]]
  })
  
  output$download_plot_taxa_based_correlation <- downloadHandler(
    filename = function(){
      paste("taxa_based_correlation_plot_", input$select_taxa_based_correlation_method, input$select_image_type_taxa_based_correlation, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_taxa_based_correlation_table()[[1]], width = input$taxa_based_correlation_output_width, height = input$taxa_based_correlation_output_height, dpi = input$taxa_based_correlation_output_dpi, units = "in")
    }
  )  
  
  ################################
  ## Correlation: Samples-based ##
  ################################ 
  
  
  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_samples_based_correlation", choices = labels_data_type)
  })
  
  data_samples_based_correlation_table <- eventReactive(input$action_samples_based_correlation,{
    source("scripts/samples_based_correlation.R")
    labels_data_type<- input$file1$name
    samples_based_correlation_plot_table(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[7]], method = input$select_samples_based_correlation_method, labe_size =input$select_samples_based_correlation_label_size)
  })
  
  output$samples_based_correlation_table <- renderDataTable(DT::datatable(data_samples_based_correlation_table()[[2]],options = list(pageLength = 15,scrollX = TRUE)))
  
  output$download_result_samples_based_correlation <- downloadHandler(
    filename = function() { 
      paste("samples_based_correlation_result_", input$select_samples_based_correlation_method, '.csv', sep='') },
    content = function(file){
      write.csv(data_samples_based_correlation_table()[[2]], file)
    }
  )
  
  output$plot_samples_based_correlation <- renderPlot({
    data_samples_based_correlation_table()[[1]]
  })
  
  output$download_plot_samples_based_correlation <- downloadHandler(
    filename = function(){
      paste("samples_based_correlation_plot_", input$select_samples_based_correlation_method, input$select_image_type_samples_based_correlation, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_samples_based_correlation_table()[[1]], width = input$samples_based_correlation_output_width, height = input$samples_based_correlation_output_height, dpi = input$samples_based_correlation_output_dpi, units = "in")
    }
  )  
  
  
  ###########################
  ##      wilcox-test      ##
  ###########################  
  
  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_wilcoxtest", choices = labels_data_type)
  })
  
  #Update conditions
  observeEvent(input$action_level,{
    label_condition1 <- "Please upload metadata in upload page"
    if(is.null(dataInput_RA_level()[[10]])){
      label_type1 <- "Please upload metadata in upload page"
    }
    else{
      label_condition1<- dataInput_RA_level()[[10]]
    }
    updateSelectInput(session, "group1_wilcoxtest", choices = label_condition1)
  })
  
  observeEvent(input$action_level,{
    label_condition2 <- "Please upload metadata in upload page"
    if(is.null(dataInput_RA_level()[[10]])){
      label_type2 <- "Please upload metadata in upload page"
    }
    else{
      label_condition2<- dataInput_RA_level()[[10]]
    }
    updateSelectInput(session, "group2_wilcoxtest", choices = label_condition2)
  })
  
  data_wilcoxtest <- eventReactive(input$action_wilcoxtest,{
    source("scripts/wilcoxtest.R")
    labels_data_type<- input$file1$name
    wilcoxtest_summary(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_wilcoxtest_pvalue, group1 = input$group1_wilcoxtest, group2 = input$group2_wilcoxtest,  plot_method = input$select_wilcoxtest_plot, alpha = input$wilcoxtest_pvalue)
    })
  
  output$text_wilcoxtest_level<- renderText({
    paste(data_wilcoxtest()[[1]])
  })
  
  output$wilcoxtest_table <- renderDataTable(DT::datatable((data_wilcoxtest()[[2]]),options = list(pageLength = 15,scrollX = TRUE)))
  
  output$download_result_wilcoxtest_1 <- downloadHandler(
    filename = function() { 
      paste("wilcoxtest_result_significant", '.csv', sep='') },
    content = function(file){
      write.csv(data_wilcoxtest()[[2]], file, row.names = FALSE)
    }
  )
  output$download_result_wilcoxtest_2 <- downloadHandler(
    filename = function() { 
      paste("wilcoxtest_result_all", '.csv', sep='') },
    content = function(file){
      write.csv(data_wilcoxtest()[[3]], file, row.names = FALSE)
    }
  )
  output$download_result_wilcoxtest_3 <- downloadHandler(
    filename = function() { 
      paste("wilcoxtest_relative_frequency", '.csv', sep='') },
    content = function(file){
      write.csv(data_wilcoxtest()[[4]], file)
    }
  )
  output$download_result_wilcoxtest_4 <- downloadHandler(
    filename = function() { 
      paste("total_counts_in_each_samples", '.csv', sep='') },
    content = function(file){
      write.csv(data_wilcoxtest()[[5]], file, row.names = FALSE)
    }
  )
  
  #boxplot
  wilcoxtest_data_plot <- eventReactive(input$action_wilcoxtest,{
    source("scripts/wilcoxtest_plot.R")
    labels_data_type<- input$file1$name
    wilcoxtest_plot(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_wilcoxtest_pvalue, group1 = input$group1_wilcoxtest, group2 = input$group2_wilcoxtest,  plot_method = input$select_wilcoxtest_plot, alpha = input$wilcoxtest_pvalue)
  })
  
  output$boxplot_wilcoxtest <- renderPlot({
    wilcoxtest_data_plot()
  })
  
  
  output$download_Boxplot_wilcoxtest<- downloadHandler(
    filename = function(){
      paste("wilcoxtest_plot", input$select_image_type_wilcoxtest, sep="")
    },
    content = function(file){
      ggsave(file,plot = wilcoxtest_data_plot(), width = input$Boxplot_wilcoxtest_output_width, height = input$Boxplot_wilcoxtest_output_height, dpi = input$Boxplot_wilcoxtest_output_dpi, units = "in")
    }
  )    
  
  
  
  
  
  
  ###########################
  ##         T-test        ##
  ###########################  
  
  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_ttest", choices = labels_data_type)
  })
  
  #Update conditions
  observeEvent(input$action_level,{
    label_condition1 <- "Please upload metadata in upload page"
    if(is.null(dataInput_RA_level()[[10]])){
      label_type1 <- "Please upload metadata in upload page"
    }
    else{
      label_condition1<- dataInput_RA_level()[[10]]
    }
    updateSelectInput(session, "group1_ttest", choices = label_condition1)
  })
  
  observeEvent(input$action_level,{
    label_condition2 <- "Please upload metadata in upload page"
    if(is.null(dataInput_RA_level()[[10]])){
      label_type2 <- "Please upload metadata in upload page"
    }
    else{
      label_condition2<- dataInput_RA_level()[[10]]
    }
    updateSelectInput(session, "group2_ttest", choices = label_condition2)
  })
  
  data_ttest <- eventReactive(input$action_ttest,{
    source("scripts/ttest.R")
    labels_data_type<- input$file1$name
    ttest_summary(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_ttest_pvalue, group1 = input$group1_ttest, group2 = input$group2_ttest,  plot_method = input$select_ttest_plot, alpha = input$ttest_pvalue)
  })
  
  output$text_ttest_level<- renderText({
    paste(data_ttest()[[1]])
  })
  
  output$ttest_table <- renderDataTable(DT::datatable((data_ttest()[[2]]), options = list( pageLength = 15, scrollX = TRUE)))
  
  output$download_result_ttest_1 <- downloadHandler(
    filename = function() { 
      paste("ttest_result_significant", '.csv', sep='') },
    content = function(file){
      write.csv(data_ttest()[[2]], file, row.names = FALSE)
    }
  )
  output$download_result_ttest_2 <- downloadHandler(
    filename = function() { 
      paste("ttest_result_all", '.csv', sep='') },
    content = function(file){
      write.csv(data_ttest()[[3]], file, row.names = FALSE)
    }
  )
  output$download_result_ttest_3 <- downloadHandler(
    filename = function() { 
      paste("ttest_relative_frequency", '.csv', sep='') },
    content = function(file){
      write.csv(data_ttest()[[4]], file)
    }
  )
  output$download_result_ttest_4 <- downloadHandler(
    filename = function() { 
      paste("total_counts_in_each_samples", '.csv', sep='') },
    content = function(file){
      write.csv(data_ttest()[[5]], file, row.names = FALSE)
    }
  )
  
  #boxplot
  ttest_data_plot <- eventReactive(input$action_ttest,{
    source("scripts/ttest_plot.R")
    labels_data_type<- input$file1$name
    ttest_plot(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_ttest_pvalue, group1 = input$group1_ttest, group2 = input$group2_ttest,  plot_method = input$select_ttest_plot, alpha = input$ttest_pvalue)
  })
  
  output$boxplot_ttest <- renderPlot({
    ttest_data_plot()
  })
  
  
  output$download_Boxplot_ttest<- downloadHandler(
    filename = function(){
      paste("ttest_plot", input$select_image_type_ttest, sep="")
    },
    content = function(file){
      ggsave(file,plot = ttest_data_plot(), width = input$Boxplot_ttest_output_width, height = input$Boxplot_ttest_output_height, dpi = input$Boxplot_ttest_output_dpi, units = "in")
    }
  )    
  
  
  
  ###########################
  ##     metagenomeSeq     ##
  ###########################  
  
  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_RA_metagenomeseq", choices = labels_data_type)
  })
  
  #Update conditions
  observeEvent(input$action_level,{
    label_condition1 <- "Please upload metadata in upload page"
    if(is.null(dataInput_RA_level()[[10]])){
      label_type1 <- "Please upload metadata in upload page"
    }
    else{
      label_condition1<- dataInput_RA_level()[[10]]
    }
    updateSelectInput(session, "group1_metagenomeseq", choices = label_condition1)
  })
  
  observeEvent(input$action_level,{
    label_condition2 <- "Please upload metadata in upload page"
    if(is.null(dataInput_RA_level()[[10]])){
      label_type2 <- "Please upload metadata in upload page"
    }
    else{
      label_condition2<- dataInput_RA_level()[[10]]
    }
    updateSelectInput(session, "group2_metagenomeseq", choices = label_condition2)
  })
  
  
  data_metagenomeseq <- eventReactive(input$action_metagenomeseq,{
    source("scripts/metagenomeseq.R")
    labels_data_type<- input$file1$name
    metagenomeseq_summary(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_metagenomeseq_pvalue, group1 = input$group1_metagenomeseq, group2 = input$group2_metagenomeseq,  plot_method = input$select_metagenomeseq_plot, alpha = input$metagenomeseq_pvalue)
  })
  
  output$text_metagenomeseq_level<- renderText({
    paste(data_metagenomeseq()[[1]])
  })
  
  output$metagenomeseq_table <- renderDataTable(DT::datatable((data_metagenomeseq()[[2]]), options = list( pageLength = 15, scrollX = TRUE)))
  
  output$download_result_metagenomeseq_1 <- downloadHandler(
    filename = function() { 
      paste("metagenomeseq_result_significant", '.csv', sep='') },
    content = function(file){
      write.csv(data_metagenomeseq()[[2]], file, row.names = FALSE)
    }
  )
  output$download_result_metagenomeseq_2 <- downloadHandler(
    filename = function() { 
      paste("metagenomeseq_result_all", '.csv', sep='') },
    content = function(file){
      write.csv(data_metagenomeseq()[[3]], file, row.names = FALSE)
    }
  )
  output$download_result_metagenomeseq_3 <- downloadHandler(
    filename = function() { 
      paste("total_counts_in_each_samples", '.csv', sep='') },
    content = function(file){
      write.csv(data_metagenomeseq()[[4]], file, row.names = FALSE)
    }
  )
  
  #boxplot
  metagenomeseq_data_plot <- eventReactive(input$action_metagenomeseq,{
    source("scripts/metagenomeseq_plot.R")
    labels_data_type<- input$file1$name
    metagenomeseq_plot(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_metagenomeseq_pvalue, group1 = input$group1_metagenomeseq, group2 = input$group2_metagenomeseq,  plot_method = input$select_metagenomeseq_plot, alpha = input$metagenomeseq_pvalue)
  })
  
  output$boxplot_metagenomeseq <- renderPlot({
    metagenomeseq_data_plot()
  })
  
  
  output$download_Boxplot_metagenomeseq<- downloadHandler(
    filename = function(){
      paste("metagenomeseq_plot", input$select_image_type_metagenomeseq, sep="")
    },
    content = function(file){
      ggsave(file,plot = metagenomeseq_data_plot(), width = input$Boxplot_metagenomeseq_output_width, height = input$Boxplot_metagenomeseq_output_height, dpi = input$Boxplot_metagenomeseq_output_dpi, units = "in")
    }
  )    
  
  
  
  ###########################
  ##         DESeq2        ##
  ###########################  
  
  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_RA_deseq2", choices = labels_data_type)
  })
  
  #Update conditions
  observeEvent(input$action_level,{
    label_condition1 <- "Please upload metadata in upload page"
    if(is.null(dataInput_RA_level()[[10]])){
      label_type1 <- "Please upload metadata in upload page"
    }
    else{
      label_condition1<- dataInput_RA_level()[[10]]
    }
    updateSelectInput(session, "group1_deseq2", choices = label_condition1)
  })
  
  observeEvent(input$action_level,{
    label_condition2 <- "Please upload metadata in upload page"
    if(is.null(dataInput_RA_level()[[10]])){
      label_type2 <- "Please upload metadata in upload page"
    }
    else{
      label_condition2<- dataInput_RA_level()[[10]]
    }
    updateSelectInput(session, "group2_deseq2", choices = label_condition2)
  })
  
  
  
  data_deseq2 <- eventReactive(input$action_deseq2,{
    source("scripts/deseq2.R")
    labels_data_type<- input$file1$name
    deseq2_summary(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_deseq2_pvalue, group1 = input$group1_deseq2, group2 = input$group2_deseq2,  plot_method = input$select_deseq2_plot, alpha = input$deseq2_pvalue)
  })
  
  output$text_deseq2_level<- renderText({
    paste(data_deseq2()[[1]])
  })
 
  output$deseq2_table <- renderDataTable(DT::datatable((data_deseq2()[[2]]), options = list( pageLength = 15, scrollX = TRUE)))
  
  output$download_result_deseq2_1 <- downloadHandler(
    filename = function() { 
      paste("deseq2_result_significant", '.csv', sep='') },
    content = function(file){
      write.csv(data_deseq2()[[2]], file, row.names = FALSE)
    }
  )
  output$download_result_deseq2_2 <- downloadHandler(
    filename = function() { 
      paste("deseq2_result_all", '.csv', sep='') },
    content = function(file){
      write.csv(data_deseq2()[[3]], file, row.names = FALSE)
    }
  )
  output$download_result_deseq2_3 <- downloadHandler(
    filename = function() { 
      paste("deseq2_normalized_count", '.csv', sep='') },
    content = function(file){
      write.csv(data_deseq2()[[4]], file)
    }
  )
  output$download_result_deseq2_4 <- downloadHandler(
    filename = function() { 
      paste("total_counts_in_each_samples", '.csv', sep='') },
    content = function(file){
      write.csv(data_deseq2()[[5]], file, row.names = FALSE)
    }
  )
  #boxplot
  deseq2_data_plot <- eventReactive(input$action_deseq2,{
    source("scripts/deseq2_plot.R")
    labels_data_type<- input$file1$name
    deseq2_plot(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_deseq2_pvalue, group1 = input$group1_deseq2, group2 = input$group2_deseq2,  plot_method = input$select_deseq2_plot, alpha = input$deseq2_pvalue)
  })
  
  output$boxplot_deseq2 <- renderPlot({
    deseq2_data_plot()
  })
  
  
  output$download_Boxplot_deseq2<- downloadHandler(
    filename = function(){
      paste("deseq2_plot", input$select_image_type_deseq2, sep="")
    },
    content = function(file){
      ggsave(file,plot = deseq2_data_plot(), width = input$Boxplot_deseq2_output_width, height = input$Boxplot_deseq2_output_height, dpi = input$Boxplot_deseq2_output_dpi, units = "in")
    }
  )    
  
  
  
  ###########################
  ##      limma-voom       ##
  ###########################  
  
  

  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_RA_limma", choices = labels_data_type)
  })
  
  #Update conditions
  observeEvent(input$action_level,{
    label_condition1 <- "Please upload metadata in upload page"
    if(is.null(dataInput_RA_level()[[10]])){
      label_type1 <- "Please upload metadata in upload page"
    }
    else{
      label_condition1<- dataInput_RA_level()[[10]]
    }
    updateSelectInput(session, "group1_limma", choices = label_condition1)
  })
  
  observeEvent(input$action_level,{
    label_condition2 <- "Please upload metadata in upload page"
    if(is.null(dataInput_RA_level()[[10]])){
      label_type2 <- "Please upload metadata in upload page"
    }
    else{
      label_condition2<- dataInput_RA_level()[[10]]
    }
    updateSelectInput(session, "group2_limma", choices = label_condition2)
  })
  
  
  
  data_limma <- eventReactive(input$action_limma,{
    source("scripts/limma.R")
    labels_data_type<- input$file1$name
    limma_summary(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_limma_pvalue, group1 = input$group1_limma, group2 = input$group2_limma,  plot_method = input$select_limma_plot, alpha = input$limma_pvalue)
  })
  
  output$text_limma_level<- renderText({
    paste(data_limma()[[1]])
  })
  
  output$limma_table <- renderDataTable(DT::datatable((data_limma()[[2]]), options = list( pageLength = 15, scrollX = TRUE)))
  
  output$download_result_limma_1 <- downloadHandler(
    filename = function() { 
      paste("limma_result_significant", '.csv', sep='') },
    content = function(file){
      write.csv(data_limma()[[2]], file, row.names = FALSE)
    }
  )
  output$download_result_limma_2 <- downloadHandler(
    filename = function() { 
      paste("limma_result_all", '.csv', sep='') },
    content = function(file){
      write.csv(data_limma()[[3]], file, row.names = FALSE)
    }
  )
  output$download_result_limma_3 <- downloadHandler(
    filename = function() { 
      paste("total_counts_in_each_samples", '.csv', sep='') },
    content = function(file){
      write.csv(data_limma()[[4]], file, row.names = FALSE)
    }
  )
 
  #boxplot
  limma_data_plot <- eventReactive(input$action_limma,{
    source("scripts/limma_plot.R")
    labels_data_type<- input$file1$name
    limma_plot(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_limma_pvalue, group1 = input$group1_limma, group2 = input$group2_limma,  plot_method = input$select_limma_plot, alpha = input$limma_pvalue)
  })
  
  output$boxplot_limma <- renderPlot({
    limma_data_plot()
  })
  
  
  output$download_Boxplot_limma<- downloadHandler(
    filename = function(){
      paste("limma_plot", input$select_image_type_limma, sep="")
    },
    content = function(file){
      ggsave(file,plot = limma_data_plot(), width = input$Boxplot_limma_output_width, height = input$Boxplot_limma_output_height, dpi = input$Boxplot_limma_output_dpi, units = "in")
    }
  )    
  
  
  
  
  ###########################
  ##         EdgeR         ##
  ###########################  
  
  
  
  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_RA_edger", choices = labels_data_type)
  })
  
  #Update conditions
  observeEvent(input$action_level,{
    label_condition1 <- "Please upload metadata in upload page"
    if(is.null(dataInput_RA_level()[[10]])){
      label_type1 <- "Please upload metadata in upload page"
    }
    else{
      label_condition1<- dataInput_RA_level()[[10]]
    }
    updateSelectInput(session, "group1_edger", choices = label_condition1)
  })
  
  observeEvent(input$action_level,{
    label_condition2 <- "Please upload metadata in upload page"
    if(is.null(dataInput_RA_level()[[10]])){
      label_type2 <- "Please upload metadata in upload page"
    }
    else{
      label_condition2<- dataInput_RA_level()[[10]]
    }
    updateSelectInput(session, "group2_edger", choices = label_condition2)
  })
  
  
  
  data_edger <- eventReactive(input$action_edger,{
    source("scripts/edger.R")
    labels_data_type<- input$file1$name
    edger_summary(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_edger_pvalue, group1 = input$group1_edger, group2 = input$group2_edger, plot_method = input$select_edger_plot, alpha = input$edger_pvalue)
  })
  
  output$text_edger_level<- renderText({
    paste(data_edger()[[1]])
  })
  
  output$edger_table <- renderDataTable(DT::datatable((data_edger()[[2]]), options = list( pageLength = 15, scrollX = TRUE)))
  
  output$download_result_edger_1 <- downloadHandler(
    filename = function() { 
      paste("edger_result_significant", '.csv', sep='') },
    content = function(file){
      write.csv(data_edger()[[2]], file, row.names = FALSE)
    }
  )
  output$download_result_edger_2 <- downloadHandler(
    filename = function() { 
      paste("edger_result_all", '.csv', sep='') },
    content = function(file){
      write.csv(data_edger()[[3]], file, row.names = FALSE)
    }
  )
  output$download_result_edger_3 <- downloadHandler(
    filename = function() { 
      paste("total_counts_in_each_samples", '.csv', sep='') },
    content = function(file){
      write.csv(data_edger()[[4]], file, row.names = FALSE)
    }
  )
  
  #boxplot
  edger_data_plot <- eventReactive(input$action_edger,{
    source("scripts/edger_plot.R")
    labels_data_type<- input$file1$name
    edger_plot(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_edger_pvalue, group1 = input$group1_edger, group2 = input$group2_edger, plot_method = input$select_edger_plot, alpha = input$edger_pvalue)
  })
  
  output$boxplot_edger <- renderPlot({
    edger_data_plot()
  })
  
  
  output$download_Boxplot_edger<- downloadHandler(
    filename = function(){
      paste("edger_plot", input$select_image_type_edger, sep="")
    },
    content = function(file){
      ggsave(file,plot = edger_data_plot(), width = input$Boxplot_edger_output_width, height = input$Boxplot_edger_output_height, dpi = input$Boxplot_edger_output_dpi, units = "in")
    }
  )   
  
  ###########################
  ##  Kruskal-Wallis Test  ##
  ###########################  
  
  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_kruskal_wallis_test", choices = labels_data_type)
  })
  
  data_kruskal_wallis_test <- eventReactive(input$action_kruskal_wallis_test,{
    source("scripts/kruskal_wallis_test.R")
    labels_data_type<- input$file1$name
    kruskal_wallis_test_summary(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_kruskal_wallis_test_pvalue, ad_hoc = input$kruskal_wallis_test_ad_hoc, plot_method = input$kruskal_wallis_test_plot, alpha = input$kruskal_wallis_test_pvalue)
  })
  
  output$text_kruskal_wallis_test_level<- renderText({
    paste(data_kruskal_wallis_test()[[1]])
  })
  
  output$kruskal_wallis_test_table <- renderDataTable(DT::datatable((data_kruskal_wallis_test()[[2]]), options = list( pageLength = 15, scrollX = TRUE)))
  
  output$download_result_kruskal_wallis_test_1 <- downloadHandler(
    filename = function() { 
      paste("kruskal_wallis_test_result_significant", '.csv', sep='') },
    content = function(file){
      write.csv(data_kruskal_wallis_test()[[2]], file, row.names = FALSE)
    }
  )
  output$download_result_kruskal_wallis_test_2 <- downloadHandler(
    filename = function() { 
      paste("kruskal_wallis_test_result_all", '.csv', sep='') },
    content = function(file){
      write.csv(data_kruskal_wallis_test()[[3]], file, row.names = FALSE)
    }
  )
  output$download_result_kruskal_wallis_test_3 <- downloadHandler(
    filename = function() { 
      paste("kruskal_wallis_test_relative_frequency", '.csv', sep='') },
    content = function(file){
      write.csv(data_kruskal_wallis_test()[[4]], file)
    }
  )
  output$download_result_kruskal_wallis_test_4 <- downloadHandler(
    filename = function() { 
      paste("total_counts_in_each_samples", '.csv', sep='') },
    content = function(file){
      write.csv(data_kruskal_wallis_test()[[5]], file, row.names = FALSE)
    }
  )
  
  #boxplot
  kruskal_wallis_test_data_plot <- eventReactive(input$action_kruskal_wallis_test,{
    source("scripts/kruskal_wallis_test_plot.R")
    labels_data_type<- input$file1$name
    kruskal_wallis_test_plot(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_kruskal_wallis_test_pvalue, ad_hoc = input$kruskal_wallis_test_ad_hoc, plot_method = input$kruskal_wallis_test_plot, alpha = input$kruskal_wallis_test_pvalue)
  })
  
  output$boxplot_kruskal_wallis_test <- renderPlot({
    kruskal_wallis_test_data_plot()
  })
  
  
  output$download_Boxplot_kruskal_wallis_test<- downloadHandler(
    filename = function(){
      paste("kruskal_wallis_test_plot", input$select_image_type_kruskal_wallis_test, sep="")
    },
    content = function(file){
      ggsave(file,plot = kruskal_wallis_test_data_plot(), width = input$Boxplot_kruskal_wallis_test_output_width, height = input$Boxplot_kruskal_wallis_test_output_height, dpi = input$Boxplot_kruskal_wallis_test_output_dpi, units = "in")
    }
  ) 
  
  ###########################
  ##         ANOVA         ##
  ###########################  
  
  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_anova", choices = labels_data_type)
  })
  
  data_anova <- eventReactive(input$action_anova,{
    source("scripts/anova.R")
    labels_data_type<- input$file1$name
    anova_summary(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_anova_pvalue, ad_hoc = input$anova_ad_hoc, plot_method = input$anova_plot, alpha = input$anova_pvalue)
  })
  
  output$text_anova_level<- renderText({
    paste(data_anova()[[1]])
  })
  
  output$anova_table <- renderDataTable(DT::datatable((data_anova()[[2]]), options = list( pageLength = 15, scrollX = TRUE)))
  
  output$download_result_anova_1 <- downloadHandler(
    filename = function() { 
      paste("anova_result_significant", '.csv', sep='') },
    content = function(file){
      write.csv(data_anova()[[2]], file, row.names = FALSE)
    }
  )
  output$download_result_anova_2 <- downloadHandler(
    filename = function() { 
      paste("anova_result_all", '.csv', sep='') },
    content = function(file){
      write.csv(data_anova()[[3]], file, row.names = FALSE)
    }
  )
  output$download_result_anova_3 <- downloadHandler(
    filename = function() { 
      paste("anova_relative_frequency", '.csv', sep='') },
    content = function(file){
      write.csv(data_anova()[[4]], file)
    }
  )
  output$download_result_anova_4 <- downloadHandler(
    filename = function() { 
      paste("total_counts_in_each_samples", '.csv', sep='') },
    content = function(file){
      write.csv(data_anova()[[5]], file, row.names = FALSE)
    }
  )
  
  #boxplot
  anova_data_plot <- eventReactive(input$action_anova,{
    source("scripts/anova_plot.R")
    labels_data_type<- input$file1$name
    anova_plot(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_anova_pvalue, ad_hoc = input$anova_ad_hoc, plot_method = input$anova_plot, alpha = input$anova_pvalue)
  })
  
  output$boxplot_anova <- renderPlot({
    anova_data_plot()
  })
  
  
  output$download_Boxplot_anova<- downloadHandler(
    filename = function(){
      paste("anova_plot", input$select_image_type_anova, sep="")
    },
    content = function(file){
      ggsave(file,plot = anova_data_plot(), width = input$Boxplot_anova_output_width, height = input$Boxplot_anova_output_height, dpi = input$Boxplot_anova_output_dpi, units = "in")
    }
  )    
  
  
  
 
  
}