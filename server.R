source("global.R")
#library("R.utils")
#options(warn=-1)
server <- function(input, output, session) {
  #Timeout
  observeEvent(input$timeOut, { 
    print(paste0("Session (", session$token, ") timed out at: ", Sys.time()))
    showModal(modalDialog(
      title = "Timeout",
      paste("Session timeout due to", input$timeOut, "inactivity -", Sys.time()),
      footer = NULL
    ))
    session$close()
  })
  
  points <- eventReactive(input$recalc, {
    cbind(rnorm(40) * 2 + 13, rnorm(40) + 48)
  }, ignoreNULL = FALSE)
  

    # Increment exactly once when the session first renders
  session$onFlushed(function() {
    current <- increment_count()
    output$view_count <- renderText(format(current, big.mark = ","))
  }, once = TRUE)

  default_if_missing <- function(x, y) {
    if (is.null(x) || length(x) == 0) {
      return(y)
    }

    if (is.character(x) && length(x) == 1 && !nzchar(x)) {
      return(y)
    }

    x
  }

  format_timestamp <- function(x) {
    format(x, "%Y-%m-%d %I:%M:%S %p")
  }

  format_scope <- function(menu, tab) {
    menu <- default_if_missing(menu, "")
    tab <- default_if_missing(tab, "")

    if (!nzchar(menu)) {
      return(tab)
    }

    if (!nzchar(tab) || identical(menu, tab)) {
      return(menu)
    }

    paste(menu, tab, sep = " > ")
  }

  status_class <- function(status) {
    tolower(gsub("[^a-z0-9]+", "-", status))
  }

  extract_ui_input_labels <- function(path = "ui.R") {
    if (!file.exists(path)) {
      return(character(0))
    }

    lines <- readLines(path, warn = FALSE)
    labels <- list()
    input_pattern <- "(?:selectInput|numericInput|radioButtons|sliderInput|checkboxInput|checkboxGroupInput|dateInput|fileInput|textInput)"
    patterns <- c(
      paste0(".*", input_pattern, "\\(\"([^\"]+)\"\\s*,\\s*label\\s*=\\s*h[1-6]\\(\"([^\"]+)\".*"),
      paste0(".*", input_pattern, "\\(\"([^\"]+)\"\\s*,\\s*label\\s*=\\s*\"([^\"]+)\".*"),
      paste0(".*", input_pattern, "\\(\"([^\"]+)\"\\s*,\\s*h[1-6]\\(\"([^\"]+)\".*"),
      paste0(".*", input_pattern, "\\(\"([^\"]+)\"\\s*,\\s*\"([^\"]+)\".*")
    )

    for (line in lines) {
      for (pattern in patterns) {
        match <- regmatches(line, regexec(pattern, line, perl = TRUE))[[1]]

        if (length(match) == 3) {
          labels[[match[2]]] <- trimws(match[3])
          break
        }
      }
    }

    unlist(labels, use.names = TRUE)
  }

  ui_input_labels <- extract_ui_input_labels("ui.R")

  format_input_value <- function(value) {
    if (is.null(value) || length(value) == 0) {
      return("Not selected")
    }

    if (is.data.frame(value) && "name" %in% names(value)) {
      value <- value$name
    } else if (is.list(value) && !is.data.frame(value)) {
      value <- unlist(value, use.names = FALSE)
    }

    if (is.null(value) || length(value) == 0) {
      return("Not selected")
    }

    value <- value[!is.na(value)]

    if (!length(value)) {
      return("Not selected")
    }

    paste(as.character(value), collapse = ", ")
  }

  format_input_label <- function(input_id) {
    default_if_missing(
      ui_input_labels[[input_id]],
      tools::toTitleCase(gsub("_", " ", input_id))
    )
  }

  collect_parameters <- function(param_ids = character(0)) {
    if (!length(param_ids)) {
      return(character(0))
    }

    vapply(
      param_ids,
      function(input_id) {
        paste0(
          format_input_label(input_id),
          ": ",
          format_input_value(isolate(input[[input_id]]))
        )
      },
      character(1)
    )
  }

  build_run_card <- function(entry) {
    params_text <- if (length(entry$parameters)) {
      paste(entry$parameters, collapse = "\n")
    } else {
      "No parameters recorded."
    }

    div(
      class = paste("run-log-card", status_class(entry$status)),
      div(class = "run-log-title", format_scope(entry$menu, entry$tab)),
      div(
        class = "run-log-meta",
        paste(
          c(
            format_timestamp(entry$timestamp),
            if (!is.na(entry$run_id)) paste0("Run #", sprintf("%03d", entry$run_id)),
            paste("Status:", entry$status)
          ),
          collapse = " | "
        )
      ),
      div(
        class = "run-log-message",
        span(class = paste("run-status-badge", status_class(entry$status)), entry$status),
        span(entry$message)
      ),
      tags$pre(class = "run-log-params", params_text)
    )
  }

  run_history <- reactiveVal(list())
  run_state <- reactiveValues(
    counter = 0L,
    current_status = "Idle",
    current_menu = "",
    current_tab = "No active analysis",
    current_message = "Load data or submit an analysis to see progress here.",
    current_parameters = character(0),
    current_timestamp = Sys.time()
  )
  action_run_cache <- reactiveValues()

  set_current_run_state <- function(status, menu = "", tab = "", message = "", parameters = character(0)) {
    run_state$current_status <- status
    run_state$current_menu <- default_if_missing(menu, "")
    run_state$current_tab <- default_if_missing(tab, "No active analysis")
    run_state$current_message <- default_if_missing(message, "")
    run_state$current_parameters <- default_if_missing(parameters, character(0))
    run_state$current_timestamp <- Sys.time()
  }

  append_run_log <- function(status, menu, tab, message, parameters = character(0), run_id = NA_integer_) {
    entry <- list(
      run_id = run_id,
      timestamp = Sys.time(),
      status = status,
      menu = menu,
      tab = tab,
      message = message,
      parameters = parameters
    )

    run_history(c(list(entry), run_history()))
    invisible(entry)
  }

  tracked_actions <- list(
    action_level = list(
      menu = "Upload files",
      tab = "Upload files",
      task = "Data upload and preprocessing",
      params = c("select_file_type", "file1", "file1_split", "metadata", "metadata_split", "select_RA_type"),
      loading_detail = "Validating the uploaded inputs",
      running_detail = "Loading and preparing the selected data"
    ),
    action_m1_bar_plot_group = list(
      menu = "Distribution",
      tab = "Group",
      task = "Group abundance distribution",
      params = c("input_RA_bar_plot_group", "select_plot_type_group", "color_palette_group", "top_n_bar_plot_group", "select_image_type_group")
    ),
    action_m1_bar_plot_individual = list(
      menu = "Distribution",
      tab = "Individual",
      task = "Sample abundance distribution",
      params = c("input_RA_bar_plot_individual", "select_plot_type_individual", "color_palette_individual", "top_n_bar_plot_individual", "select_image_type_individual")
    ),
    action_alpha_diversity = list(
      menu = "Diversity",
      tab = "Alpha",
      task = "Alpha diversity analysis",
      params = c("input_RA_Alpha", "select_alpha", "select_alpha_pvalue", "select_plot", "select_alpha_color_palette", "select_image_type_alpha")
    ),
    action_beta_diversity = list(
      menu = "Diversity",
      tab = "Beta",
      task = "Beta diversity analysis",
      params = c("input_RA_Beta", "select_beta", "adonis_permutations", "adonis_dissimilarities", "select_method", "select_beta_color_palette", "select_image_type_beta")
    ),
    action_pca = list(
      menu = "Dimension reduction",
      tab = "PCA-2D",
      task = "PCA-2D analysis",
      params = c("input_RA_pca", "select_pca_label", "select_pca_label_size", "select_pca_frame", "select_pca_color_palette", "select_image_type_pca")
    ),
    action_pca3d = list(
      menu = "Dimension reduction",
      tab = "PCA-3D",
      task = "PCA-3D analysis",
      params = c("input_RA_pca3d", "select_pca3d_color_palette")
    ),
    action_tsne = list(
      menu = "Dimension reduction",
      tab = "t-SNE",
      task = "t-SNE analysis",
      params = c("input_tsne", "select_tsne_method", "select_tsne_dimension", "select_tsne_color_palette", "select_image_type_tsne")
    ),
    action_umap = list(
      menu = "Dimension reduction",
      tab = "UMAP",
      task = "UMAP analysis",
      params = c("input_umap", "select_umap_method", "select_umap_kvalue", "select_umap_color_palette", "select_image_type_umap")
    ),
    action_taxa_condition_based_correlation = list(
      menu = "Correlation",
      tab = "Taxa-based",
      task = "Taxa-condition correlation analysis",
      params = c("input_taxa_condition_based_correlation", "input_taxa_condition_based_correlation2", "select_taxa_condition_based_correlation_method", "select_taxa_condition_based_correlation_label_size", "select_taxa_condition_geom_shape", "select_image_taxa_condition_based_correlation")
    ),
    action_samples_based_correlation = list(
      menu = "Correlation",
      tab = "Sample-based",
      task = "Sample correlation analysis",
      params = c("input_samples_based_correlation", "input_samples_based_correlation2", "select_samples_based_correlation_method", "select_samples_based_correlation_label_size", "select_sample_geom_shape", "select_image_type_samples_based_correlation")
    ),
    action_heatmap = list(
      menu = "Heatmap",
      tab = "Heatmap",
      task = "Heatmap analysis",
      params = c("input_heatmap", "heatmap_clustering_method_rows", "heatmap_clustering_method_columns", "heatmap_normalization", "heatmap_color_palette", "heatmap_row_names", "heatmap_row_names_size", "heatmap_column_names", "heatmap_column_names_size", "heatmap_row_dend", "heatmap_column_dend", "select_image_type_heatmap")
    ),
    action_wilcoxtest = list(
      menu = "Differential abundance",
      tab = "Wilcoxon Rank Sum test",
      task = "Wilcoxon Rank Sum test",
      params = c("input_wilcoxtest", "group1_wilcoxtest", "group2_wilcoxtest", "select_wilcoxtest_pvalue", "wilcoxtest_pvalue", "wilcox_color_palette", "select_wilcoxtest_plot", "select_image_type_wilcoxtest")
    ),
    action_ttest = list(
      menu = "Differential abundance",
      tab = "t-test",
      task = "t-test",
      params = c("input_ttest", "group1_ttest", "group2_ttest", "select_ttest_pvalue", "ttest_pvalue", "ttest_color_palette", "select_ttest_plot", "select_image_type_ttest")
    ),
    action_metagenomeseq = list(
      menu = "Differential abundance",
      tab = "metagenomeSeq",
      task = "metagenomeSeq analysis",
      params = c("input_RA_metagenomeseq", "group1_metagenomeseq", "group2_metagenomeseq", "select_metagenomeseq_pvalue", "metagenomeseq_pvalue", "metagenomeseq_color_palette", "select_metagenomeseq_plot", "select_image_type_metagenomeseq")
    ),
    action_deseq2 = list(
      menu = "Differential abundance",
      tab = "DESeq2",
      task = "DESeq2 analysis",
      params = c("input_RA_deseq2", "group1_deseq2", "group2_deseq2", "select_deseq2_pvalue", "deseq2_pvalue", "deseq2_color_palette", "select_deseq2_plot", "select_image_type_deseq2")
    ),
    action_LEfSe = list(
      menu = "Differential abundance",
      tab = "LEfSe",
      task = "LEfSe analysis",
      params = c("input_RA_LEfSe", "group1_LEfSe", "group2_LEfSe", "select_LEfSe_method", "select_LEfSe_pvalue", "select_LEfSe_threshold", "select_LEfSe_color_palette", "select_image_type_LEfSe")
    ),
    action_MaAsLin3 = list(
      menu = "Differential abundance",
      tab = "MaAsLin3",
      task = "MaAsLin3 analysis",
      params = c("input_RA_MaAsLin3", "group1_MaAsLin3", "group2_MaAsLin3", "select_MaAsLin3_normalization", "select_MaAsLin3_transformation", "select_MaAsLin3_correction", "select_MaAsLin3_pvalue")
    ),
    action_limma = list(
      menu = "Differential abundance",
      tab = "Limma-Voom",
      task = "Limma-Voom analysis",
      params = c("input_RA_limma", "group1_limma", "group2_limma", "select_limma_pvalue", "limma_pvalue", "limma_color_palette", "select_limma_plot", "select_image_type_limma")
    ),
    action_edger = list(
      menu = "Differential abundance",
      tab = "edgeR",
      task = "edgeR analysis",
      params = c("input_RA_edger", "group1_edger", "group2_edger", "select_edger_pvalue", "edger_pvalue", "edger_color_palette", "select_edger_plot", "select_image_type_edger")
    ),
    action_kruskal_wallis_test = list(
      menu = "Differential abundance",
      tab = "Kruskal-Wallis test",
      task = "Kruskal-Wallis test",
      params = c("input_kruskal_wallis_test", "select_kruskal_wallis_test_pvalue", "kruskal_wallis_test_pvalue", "kruskal_wallis_test_ad_hoc", "kruskal_wallis_test_color_palette", "kruskal_wallis_test_plot", "select_image_type_kruskal_wallis_test")
    ),
    action_anova = list(
      menu = "Differential abundance",
      tab = "ANOVA",
      task = "ANOVA",
      params = c("input_anova", "select_anova_pvalue", "anova_pvalue", "anova_ad_hoc", "anova_color_palette", "anova_plot", "select_image_type_anova")
    )
  )

  run_with_tracking <- function(action_id, runner) {
    cfg <- tracked_actions[[action_id]]

    if (is.null(cfg)) {
      stop("No tracking configuration found for action: ", action_id)
    }

    action_click <- isolate(input[[action_id]])
    cache_key <- paste(action_id, action_click, sep = "::")
    failure_notification_id <- paste0("run-failed-", action_id)
    cached_run <- isolate(action_run_cache[[cache_key]])

    if (!is.null(cached_run)) {
      if (identical(cached_run$status, "success")) {
        return(cached_run$value)
      }

      validate(need(FALSE, cached_run$message))
    }

    run_id <- isolate(run_state$counter) + 1L
    params <- collect_parameters(cfg$params)
    loading_message <- default_if_missing(cfg$loading_detail, paste("Loading", tolower(cfg$task), "..."))
    running_message <- default_if_missing(cfg$running_detail, paste("Running", tolower(cfg$task), "..."))
    completed_message <- paste(cfg$task, "completed successfully.")
    progress_title <- paste(format_scope(cfg$menu, cfg$tab), "Running")

    run_state$counter <- run_id
    removeNotification(failure_notification_id)
    set_current_run_state("Loading", cfg$menu, cfg$tab, loading_message, params)
    append_run_log("Loading", cfg$menu, cfg$tab, loading_message, params, run_id)

    tryCatch({
      progress <- shiny::Progress$new(session, min = 0, max = 1, style = "notification")
      on.exit(progress$close(), add = TRUE)
      progress$set(message = progress_title, detail = loading_message, value = 0)
      progress$inc(0.2, detail = loading_message)
      
      result <- {
        set_current_run_state("Running", cfg$menu, cfg$tab, running_message, params)
        append_run_log("Running", cfg$menu, cfg$tab, running_message, params, run_id)
        progress$inc(0.55, detail = running_message)
        result_value <- runner()
        progress$inc(0.25, detail = "Completed successfully.")
        result_value
      }

      action_run_cache[[cache_key]] <- list(status = "success", value = result)
      set_current_run_state("Completed successfully", cfg$menu, cfg$tab, completed_message, params)
      append_run_log("Completed successfully", cfg$menu, cfg$tab, completed_message, params, run_id)
      removeNotification(failure_notification_id)
      showNotification(completed_message, type = "message", duration = 5)
      result
    }, error = function(err) {
      failed_message <- if (inherits(err, "shiny.silent.error")) {
        paste(cfg$task, "needs valid input before it can run.")
      } else {
        paste(cfg$task, "failed:", conditionMessage(err))
      }

      action_run_cache[[cache_key]] <- list(status = "failed", message = failed_message)
      set_current_run_state("Failed", cfg$menu, cfg$tab, failed_message, params)
      append_run_log("Failed", cfg$menu, cfg$tab, failed_message, params, run_id)
      showNotification(
        failed_message,
        id = failure_notification_id,
        type = if (inherits(err, "shiny.silent.error")) "warning" else "error",
        duration = NULL,
        closeButton = TRUE
      )
      
      validate(need(FALSE, failed_message))
    })
  }

  output$run_status_panel <- renderUI({
    build_run_card(list(
      run_id = if (run_state$counter > 0) run_state$counter else NA_integer_,
      timestamp = run_state$current_timestamp,
      status = run_state$current_status,
      menu = run_state$current_menu,
      tab = run_state$current_tab,
      message = run_state$current_message,
      parameters = run_state$current_parameters
    ))
  })

  output$run_log_ui <- renderUI({
    entries <- run_history()

    if (!length(entries)) {
      return(div(class = "run-log-empty", "Run history will appear here after data loading or analysis submission."))
    }

    do.call(tagList, lapply(entries, build_run_card))
  })

  observeEvent(input$main_navbar, ignoreInit = TRUE, {
    selected_tab <- input$main_navbar

    if (is.null(selected_tab) || !nzchar(selected_tab) || identical(selected_tab, "Run Log")) {
      return()
    }

    append_run_log(
      status = "Viewed",
      menu = "Navigation",
      tab = as.character(selected_tab),
      message = paste("Opened", as.character(selected_tab), "tab."),
      parameters = character(0)
    )
  })
  
######session info
  
  # --- server ---
  # Show session info in the tab and enable download as .txt
  sess_txt <- reactive({
    paste(capture.output(utils::sessionInfo()), collapse = "\n")
  })
  
  output$sess <- renderPrint({
    cat(sess_txt())
  })
  
  output$download_sess <- downloadHandler(
    filename = function() {
      paste0("MetaDAVis_session-info_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".txt")
    },
    content = function(file) {
      writeLines(sess_txt(), con = file, useBytes = TRUE)
    }
  )
  
###########################
##         Input         ##
###########################  
  observe({
    if (input$select_file_type == "example") {
      shinyjs::hide("box1")
      shinyjs::hide("box2")
      }
    else {
      shinyjs::show("box1")
      shinyjs::show("box2") 
    }
  })
  
  dataInput_RA_level <- eventReactive(input$action_level, {
    run_with_tracking("action_level", function() {
      # Check the selected file type
      if (input$select_file_type %in% c("Megan", "check", "qiime_format")) {
        inFile1 <- input$file1
        inFile2 <- input$metadata
        sep <- input$file1_split
        sep1 <- input$metadata_split
        type <- input$select_file_type
        
        req(inFile1, inFile2)  # Ensure files are uploaded
        
        # Check file extension and validate
        ext <- tools::file_ext(inFile1$name)
        ext1 <- tools::file_ext(inFile2$name)

        validate(
          need(
            ((type %in% c("Megan", "check") && ext == "tsv" && sep == "\t") ||
               (type == "qiime_format" && ext == "csv" && sep == ",")),
            "Please select the proper input format and 'Fields separated by'."
          ),
          need(
            ((type %in% c("Megan", "check") && ext1 == "tsv" && sep1 == "\t") ||
               (type == "qiime_format" && ext1 == "csv" && sep1 == ",")),
            "Please select the proper input format and 'Fields separated by'."
          )
        )
      } else if (input$select_file_type == "example") {
        # Load example data
        inFile1 <- "www/example_data/Megan_WGS_output.tsv"
        inFile2 <- "www/example_data/Megan_WGS_metadata.tsv"
      } else {
        return(list("Input is missing", "Input is missing", 30, "Input is missing"))
      }
      
      # Check if files are uploaded
      if (is.null(inFile1) || is.null(inFile2)) {
        return(list("Input is missing", "Input is missing", 30, "Input is missing"))
      }
      
      # Call the data_input_RA function with provided inputs
      if (input$select_file_type == "example") {
        data_input_RA(
          Input = inFile1,
          Index = inFile2,
          type = input$select_RA_type,
          # sep = sep,
          # sep1 = sep1,
          file_type = input$select_file_type
        )
      } else {
        data_input_RA(
          Input = inFile1$datapath,
          Index = inFile2$datapath,
          type = input$select_RA_type,
          sep = input$file1_split,
          sep1 = input$metadata_split,
          file_type = input$select_file_type
        )
      }
    })
  })
  
  
  output$text_level<- renderText({
    paste(dataInput_RA_level()[[1]])
  })
  output$text_metadata<- renderText({
    paste(dataInput_RA_level()[[2]])
  })
  output$taxonomy_table <- renderDataTable(DT::datatable(dataInput_RA_level()[[8]],options = list(pageLength = 15,scrollX = TRUE)))
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
    run_with_tracking("action_m1_bar_plot_group", function() {
      source("scripts/group_bar_plot.R")
      Bar_Plot_Group_OTU(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[5]], color_palette_group = input$color_palette_group, n_top = input$top_n_bar_plot_group, names = dataInput_RA_level()[[6]], plot_type=input$select_plot_type_group)
    })
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
    run_with_tracking("action_m1_bar_plot_individual", function() {
      source("scripts/individual_bar_plot.R")
      Bar_Plot_Individual_OTU(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], color_palette_individual = input$color_palette_individual, n_top = input$top_n_bar_plot_individual, names = dataInput_RA_level()[[6]], plot_type=input$select_plot_type_individual)
    })
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
    run_with_tracking("action_heatmap", function() {
      source("scripts/heatmap.R")
      heatmap_abundance(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], heatmap_clustering_method_rows = input$heatmap_clustering_method_rows, heatmap_clustering_method_columns = input$heatmap_clustering_method_columns, heatmap_color_palette = input$heatmap_color_palette, heatmap_normalization = input$heatmap_normalization, heatmap_row_names = input$heatmap_row_names, heatmap_row_names_size = input$heatmap_row_names_size, heatmap_column_names = input$heatmap_column_names, heatmap_column_names_size = input$heatmap_column_names_size, heatmap_row_dend = input$heatmap_row_dend, heatmap_column_dend = input$heatmap_column_dend)
    })
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
  
  data_Alpha_Div_plot <- eventReactive(input$action_alpha_diversity,{
    run_with_tracking("action_alpha_diversity", function() {
      source("scripts/alpha_diversity.R")
      labels_data_type<- input$file1$name
      alpha_diversity_boxplot(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[7]], tax_index = dataInput_RA_level()[[8]], index_method = input$select_alpha, select_alpha_color_palette = input$select_alpha_color_palette, plot_method = input$select_plot, pvalue=input$select_alpha_pvalue)
    })
  })
  
  output$boxplot_Alpha_Div <- renderPlot({
    data_Alpha_Div_plot()[1]
  })
  
  
  output$download_Boxplot_Alpha_Div<- downloadHandler(
    filename = function(){
      method_tmp <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher", "All_Combined")
      paste("Alpha_Diversity_", method_tmp[as.numeric(input$select_alpha)],"_index", input$select_image_type_alpha, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_Alpha_Div_plot()[[1]], width = input$Boxplot_alpha_div_output_width, height = input$Boxplot_alpha_div_output_height, dpi = input$Boxplot_alpha_div_output_dpi, units = "in")
    }
  )
  
  output$alpha_table <- renderDataTable(DT::datatable(data_Alpha_Div_plot()[[2]],options = list(pageLength = 15,scrollX = TRUE)))
  
  
  output$download_result_alpha <- downloadHandler(
    filename = function() { 
      paste("Alpha_diversity_result", '.csv', sep='') },
    content = function(file){
      write.csv(data_Alpha_Div_plot()[[2]], file)
    }
  )

  
  ###########################
  ##     Beta Diversity    ##
  ########################### 
  

  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_RA_Beta", choices = labels_data_type)
  })
  
  data_Beta_Div_plot <- eventReactive(input$action_beta_diversity,{
    run_with_tracking("action_beta_diversity", function() {
      source("scripts/beta_diversity.R")
      labels_data_type<- input$file1$name
      beta_diversity_boxplot(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[7]], tax_index = dataInput_RA_level()[[8]], select_beta_color_palette=input$select_beta_color_palette, index_method = input$select_beta, plot_method = input$select_method, adonis_dissimilarities = input$adonis_dissimilarities, adonis_permutations = input$adonis_permutations)
    })
  })
  
  output$boxplot_Beta_Div <- renderPlot({
    data_Beta_Div_plot()[1]
  })
  
  output$download_Boxplot_beta_Div<- downloadHandler(
    filename = function(){
      paste("Beta_Diversity_", input$select_beta, "_",input$select_method, input$select_image_type_beta, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_Beta_Div_plot()[[1]], width = input$Boxplot_beta_div_output_width, height = input$Boxplot_beta_div_output_height, dpi = input$Boxplot_beta_div_output_dpi, units = "in")
    }
  )


  output$beta_table <- renderDataTable(DT::datatable(data_Beta_Div_plot()[[2]],options = list(pageLength = 15,scrollX = TRUE)))
  
  output$download_result_beta <- downloadHandler(
    filename = function() { 
      paste("Beta_diversity_result", '.csv', sep='') },
    content = function(file){
      write.csv(data_Beta_Div_plot()[[2]], file)
    }
  )
  
  output$beta_table2 <- renderDataTable(DT::datatable(data_Beta_Div_plot()[[3]],options = list(pageLength = 15,scrollX = TRUE)))
  
  output$download_result_beta2 <- downloadHandler(
    filename = function() { 
      paste("Beta_diversity_adonis_result", '.csv', sep='') },
    content = function(file){
      write.csv(data_Beta_Div_plot()[[3]], file)
    }
  )
  
  
  ###########################
  ##         PCA-2D        ##
  ###########################  
  
  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_RA_pca", choices = labels_data_type)
  })
  
  
  data_pca_plot <- eventReactive(input$action_pca,{
    run_with_tracking("action_pca", function() {
      source("scripts/pca_plot.R")
      labels_data_type<- input$file1$name
      pca_plot(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[7]], pca_label1 = input$select_pca_label, pca_label_size = input$select_pca_label_size, pca_frame = input$select_pca_frame, select_pca_color_palette=input$select_pca_color_palette)
    })
  })
  
  output$plot_pca <- renderPlot({
    data_pca_plot()[1]
  })
  
  output$download_plot_pca <- downloadHandler(
    filename = function(){
     paste("PCA_plot", input$select_image_type_pca, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_pca_plot()[[1]], width = input$pca_plot_output_width, height = input$pca_plot_output_height, dpi = input$pca_plot_output_dpi, units = "in")
    }
  )  
  
  output$pca_table <- renderDataTable(DT::datatable(data_pca_plot()[[2]],options = list(pageLength = 15,scrollX = TRUE)))
  
  output$download_result_pca <- downloadHandler(
    filename = function() { 
      paste("pca_result", '.csv', sep='') },
    content = function(file){
      write.csv(data_pca_plot()[[2]], file)
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
    run_with_tracking("action_pca3d", function() {
      source("scripts/pca3d_plot.R")
      labels_data_type<- input$file1$name
      pca3d_plot(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[7]], select_pca3d_color_palette=input$select_pca3d_color_palette)
    })
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
    run_with_tracking("action_tsne", function() {
      source("scripts/tsne.R")
      labels_data_type<- input$file1$name
      tsne_plot_table(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[7]], tax_index = dataInput_RA_level()[[8]], method = input$select_tsne_method, dimension = input$select_tsne_dimension, select_tsne_color_palette = input$select_tsne_color_palette)
    })
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
    run_with_tracking("action_umap", function() {
      source("scripts/umap.R")
      labels_data_type<- input$file1$name
      umap_plot_table(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[7]], tax_index = dataInput_RA_level()[[8]], method = input$select_umap_method, kvalue = input$select_umap_kvalue, select_umap_color_palette=input$select_umap_color_palette)
    })
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
  
  
  # observe({
  #   labels_data_type<- input$file1$name
  #   updateSelectInput(session, "input_taxa_based_correlation", choices = labels_data_type)
  # })
  # 
  # data_taxa_based_correlation_table <- eventReactive(input$action_taxa_based_correlation,{
  #   source("scripts/taxa_based_correlation.R")
  #   labels_data_type<- input$file1$name
  #   taxa_based_correlation_plot_table(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[7]], method = input$select_taxa_based_correlation_method, labe_size =input$select_taxa_based_correlation_label_size,  select_taxa_geom_shape=input$select_taxa_geom_shape)
  # })
  # 
  # output$taxa_based_correlation_table <- renderDataTable(DT::datatable(data_taxa_based_correlation_table()[[2]],options = list(pageLength = 15,scrollX = TRUE)))
  # 
  # output$download_result_taxa_based_correlation <- downloadHandler(
  #   filename = function() { 
  #     paste("taxa_based_correlation_result_", input$select_taxa_based_correlation_method, '.csv', sep='') },
  #   content = function(file){
  #     write.csv(data_taxa_based_correlation_table()[[2]], file)
  #   }
  # )
  # 
  # output$plot_taxa_based_correlation <- renderPlot({
  #   data_taxa_based_correlation_table()[[1]]
  # })
  # 
  # output$download_plot_taxa_based_correlation <- downloadHandler(
  #   filename = function(){
  #     paste("taxa_based_correlation_plot_", input$select_taxa_based_correlation_method, input$select_image_type_taxa_based_correlation, sep="")
  #   },
  #   content = function(file){
  #     ggsave(file,plot = data_taxa_based_correlation_table()[[1]], width = input$taxa_based_correlation_output_width, height = input$taxa_based_correlation_output_height, dpi = input$taxa_based_correlation_output_dpi, units = "in")
  #   }
  # )  
  
  ################################
  ## Correlation: Samples-based ##
  ################################ 
  
  
  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_samples_based_correlation", choices = labels_data_type)
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
    updateSelectInput(session, "input_samples_based_correlation2", choices = label_condition1, selected = label_condition1[1])
  })
  
  data_samples_based_correlation_table <- eventReactive(input$action_samples_based_correlation,{
    run_with_tracking("action_samples_based_correlation", function() {
      source("scripts/samples_based_correlation.R")
      labels_data_type<- input$file1$name
      samples_based_correlation_plot_table(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], group1 = input$input_samples_based_correlation2, method = input$select_samples_based_correlation_method, labe_size =input$select_samples_based_correlation_label_size, select_sample_geom_shape=input$select_sample_geom_shape)
    })
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
  
  ##########################################
  ##   Correlation: Taxa-condition based  ##
  ##########################################  
  
  
  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_taxa_condition_based_correlation", choices = labels_data_type)
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
    updateSelectInput(session, "input_taxa_condition_based_correlation2", choices = label_condition1, selected = label_condition1[1])
  })
  
  data_taxa_condition_based_correlation_table <- eventReactive(input$action_taxa_condition_based_correlation,{
    run_with_tracking("action_taxa_condition_based_correlation", function() {
      source("scripts/taxa_condition_based_correlation.R")
      labels_data_type<- input$file1$name
      taxa_condition_based_correlation_plot_table(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], group1 = input$input_taxa_condition_based_correlation2, method = input$select_taxa_condition_based_correlation_method, labe_size =input$select_taxa_condition_based_correlation_label_size,  select_taxa_condition_geom_shape=input$select_taxa_condition_geom_shape)
    })
  })
  
  output$taxa_condition_based_correlation_table <- renderDataTable(DT::datatable(data_taxa_condition_based_correlation_table()[[2]],options = list(pageLength = 15,scrollX = TRUE)))
  
  output$download_result_taxa_condition_based_correlation <- downloadHandler(
    filename = function() { 
      paste("taxa_condition_based_correlation_result_", input$select_taxa_condition_based_correlation_method, '.csv', sep='') },
    content = function(file){
      write.csv(data_taxa_condition_based_correlation_table()[[2]], file)
    }
  )

  taxa_condition_based_correlation_plot_dims <- reactive({
    plot_data <- data_taxa_condition_based_correlation_table()[[2]]

    if (is.null(plot_data) || nrow(plot_data) == 0) {
      return(list(width = 1100, height = 1000))
    }

    axis_labels <- unique(c(as.character(plot_data$x), as.character(plot_data$y)))
    axis_labels <- axis_labels[!is.na(axis_labels) & nzchar(axis_labels)]

    label_count <- length(axis_labels)
    max_label_length <- if (label_count > 0) max(nchar(axis_labels), na.rm = TRUE) else 0

    list(
      width = max(1100, min(3200, 220 + label_count * 22 + max_label_length * 7)),
      height = max(1000, min(3200, 220 + label_count * 20 + max_label_length * 4))
    )
  })

  output$taxa_condition_based_correlation_plot_ui <- renderUI({
    req(input$action_taxa_condition_based_correlation > 0)

    dims <- taxa_condition_based_correlation_plot_dims()

    withSpinner(
      div(
        style = "width: 100%; height: 20%; overflow-x: auto;",
        plotOutput(
          "plot_taxa_condition_based_correlation",
          width = paste0(dims$width, "px"),
          height = paste0(dims$height, "px")
        )
      )
    )
  })
  
  output$plot_taxa_condition_based_correlation <- renderPlot({
    data_taxa_condition_based_correlation_table()[[1]]
  })
  
  output$download_plot_taxa_condition_based_correlation <- downloadHandler(
    filename = function(){
      paste("taxa_condition_based_correlation_plot_", input$select_taxa_condition_based_correlation_method, input$select_image_taxa_condition_based_correlation, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_taxa_condition_based_correlation_table()[[1]], width = input$taxa_condition_based_correlation_output_width, height = input$taxa_condition_based_correlation_output_height, dpi = input$taxa_condition_based_correlation_output_dpi, units = "in")
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
    run_with_tracking("action_wilcoxtest", function() {
      source("scripts/wilcoxtest.R")
      labels_data_type<- input$file1$name
      wilcoxtest_summary(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_wilcoxtest_pvalue, group1 = input$group1_wilcoxtest, group2 = input$group2_wilcoxtest,  plot_method = input$select_wilcoxtest_plot, alpha = input$wilcoxtest_pvalue, wilcox_color_palette=input$wilcox_color_palette)
    })
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
  

  output$boxplot_wilcoxtest <- renderPlot({
    data_wilcoxtest()[6]
  })
  
  
  output$download_Boxplot_wilcoxtest<- downloadHandler(
    filename = function(){
      paste("wilcoxtest_plot", input$select_image_type_wilcoxtest, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_wilcoxtest()[[6]], width = input$Boxplot_wilcoxtest_output_width, height = input$Boxplot_wilcoxtest_output_height, dpi = input$Boxplot_wilcoxtest_output_dpi, units = "in")
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
    run_with_tracking("action_ttest", function() {
      source("scripts/ttest.R")
      labels_data_type<- input$file1$name
      ttest_summary(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_ttest_pvalue, group1 = input$group1_ttest, group2 = input$group2_ttest,  plot_method = input$select_ttest_plot, alpha = input$ttest_pvalue, ttest_color_palette = input$ttest_color_palette)
    })
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
  

  output$boxplot_ttest <- renderPlot({
    data_ttest()[6]
  })
  
  
  output$download_Boxplot_ttest<- downloadHandler(
    filename = function(){
      paste("ttest_plot", input$select_image_type_ttest, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_ttest()[[6]], width = input$Boxplot_ttest_output_width, height = input$Boxplot_ttest_output_height, dpi = input$Boxplot_ttest_output_dpi, units = "in")
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
    run_with_tracking("action_metagenomeseq", function() {
      source("scripts/metagenomeseq.R")
      labels_data_type<- input$file1$name
      metagenomeseq_summary(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_metagenomeseq_pvalue, group1 = input$group1_metagenomeseq, group2 = input$group2_metagenomeseq,  plot_method = input$select_metagenomeseq_plot, alpha = input$metagenomeseq_pvalue, metagenomeseq_color_palette = input$metagenomeseq_color_palette)
    })
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
  
  output$boxplot_metagenomeseq <- renderPlot({
    data_metagenomeseq()[5]
  })
  
  
  output$download_Boxplot_metagenomeseq<- downloadHandler(
    filename = function(){
      paste("metagenomeseq_plot", input$select_image_type_metagenomeseq, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_metagenomeseq()[[5]], width = input$Boxplot_metagenomeseq_output_width, height = input$Boxplot_metagenomeseq_output_height, dpi = input$Boxplot_metagenomeseq_output_dpi, units = "in")
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
    run_with_tracking("action_deseq2", function() {
      source("scripts/deseq2.R")
      labels_data_type<- input$file1$name
      deseq2_summary(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_deseq2_pvalue, group1 = input$group1_deseq2, group2 = input$group2_deseq2,  plot_method = input$select_deseq2_plot, alpha = input$deseq2_pvalue, deseq2_color_palette = input$deseq2_color_palette)
    })
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
  output$boxplot_deseq2 <- renderPlot({
    data_deseq2()[6]
  })
  
  
  output$download_Boxplot_deseq2<- downloadHandler(
    filename = function(){
      paste("deseq2_plot", input$select_image_type_deseq2, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_deseq2()[[6]], width = input$Boxplot_deseq2_output_width, height = input$Boxplot_deseq2_output_height, dpi = input$Boxplot_deseq2_output_dpi, units = "in")
    }
  )    
  
  
  
  ###########################
  ##         LEfSe        ##
  ###########################  
  
  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_RA_LEfSe", choices = labels_data_type)
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
    updateSelectInput(session, "group1_LEfSe", choices = label_condition1)
  })
  
  observeEvent(input$action_level,{
    label_condition2 <- "Please upload metadata in upload page"
    if(is.null(dataInput_RA_level()[[10]])){
      label_type2 <- "Please upload metadata in upload page"
    }
    else{
      label_condition2<- dataInput_RA_level()[[10]]
    }
    updateSelectInput(session, "group2_LEfSe", choices = label_condition2)
  })
  
  
  
  data_LEfSe <- eventReactive(input$action_LEfSe,{
    run_with_tracking("action_LEfSe", function() {
      source("scripts/LEfSe.R")
      labels_data_type<- input$file1$name
      LEfSe_summary(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], group1 = input$group1_LEfSe, group2 = input$group2_LEfSe,  select_LEfSe_method = input$select_LEfSe_method, select_LEfSe_pvalue = input$select_LEfSe_pvalue, select_LEfSe_threshold = input$select_LEfSe_threshold, select_LEfSe_color_palette = input$select_LEfSe_color_palette)
    })
  })
  
  output$text_LEfSe_level<- renderText({
    paste(data_LEfSe()[[1]])
  })
  
  output$LEfSe_table <- renderDataTable(DT::datatable((data_LEfSe()[[2]]), options = list( pageLength = 15, scrollX = TRUE)))
  
  output$download_result_LEfSe_1 <- downloadHandler(
    filename = function() { 
      paste("LEfSe_result_significant", '.csv', sep='') },
    content = function(file){
      write.csv(data_LEfSe()[[2]], file, row.names = FALSE)
    }
  )
  
  output$download_result_LEfSe_4 <- downloadHandler(
    filename = function() { 
      paste("total_counts_in_each_samples", '.csv', sep='') },
    content = function(file){
      write.csv(data_LEfSe()[[3]], file, row.names = FALSE)
    }
  )
  output$boxplot_LEfSe <- renderPlot({
    data_LEfSe()[4]
  })
  
  
  output$download_Boxplot_LEfSe<- downloadHandler(
    filename = function(){
      paste("LEfSe_plot", input$select_image_type_LEfSe, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_LEfSe()[[4]], width = input$Boxplot_LEfSe_output_width, height = input$Boxplot_LEfSe_output_height, dpi = input$Boxplot_LEfSe_output_dpi, units = "in")
    }
  )  
  
  
  ###########################
  ##        MaAsLin3       ##
  ###########################  
  
  observe({
    labels_data_type<- input$file1$name
    updateSelectInput(session, "input_RA_MaAsLin3", choices = labels_data_type)
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
    updateSelectInput(session, "group1_MaAsLin3", choices = label_condition1)
  })
  
  observeEvent(input$action_level,{
    label_condition2 <- "Please upload metadata in upload page"
    if(is.null(dataInput_RA_level()[[10]])){
      label_type2 <- "Please upload metadata in upload page"
    }
    else{
      label_condition2<- dataInput_RA_level()[[10]]
    }
    updateSelectInput(session, "group2_MaAsLin3", choices = label_condition2)
  })
  
  data_MaAsLin3 <- eventReactive(input$action_MaAsLin3,{
    run_with_tracking("action_MaAsLin3", function() {
      files_to_delete <- c(
        "www/hmp2_output",
        "www/hmp2_output.zip"
      )

      for (file in files_to_delete) {
        full_path <- file.path(getwd(), file)

        if (file.exists(full_path)) {
          unlink(full_path, recursive = TRUE, force = TRUE)
          cat("Deleted:", full_path, "\n")
        } else {
          cat("File or folder not found, skipping:", full_path, "\n")
        }
      }

      source("scripts/MaAsLin3.R")
      labels_data_type<- input$file1$name
      MaAsLin3_summary(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], group1 = input$group1_MaAsLin3, group2 = input$group2_MaAsLin3,  select_MaAsLin3_normalization = input$select_MaAsLin3_normalization, select_MaAsLin3_transformation = input$select_MaAsLin3_transformation, select_MaAsLin3_correction = input$select_MaAsLin3_correction, select_MaAsLin3_pvalue = input$select_MaAsLin3_pvalue)
    })
  })
  
  output$text_MaAsLin3_level<- renderText({
    paste(data_MaAsLin3()[[1]])
  })

  
  #temp_dir <- getwd()
  zip_path <- file.path(getwd(),"/www/hmp2_output.zip")
  # Provide the ZIP file for download
  output$download_zip_MaAsLin3 <- downloadHandler(
    filename = function() {
      "Maaslin3_output.zip"
    },
    content = function(file) {
      # Copy the zip file to the location Shiny specifies
      file.copy(zip_path, file)
    },
    contentType = "application/zip"
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
    run_with_tracking("action_limma", function() {
      source("scripts/limma.R")
      labels_data_type<- input$file1$name
      limma_summary(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_limma_pvalue, group1 = input$group1_limma, group2 = input$group2_limma,  plot_method = input$select_limma_plot, alpha = input$limma_pvalue, limma_color_palette = input$limma_color_palette)
    })
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
 
  
  output$boxplot_limma <- renderPlot({
    data_limma()[5]
  })
  
  
  output$download_Boxplot_limma<- downloadHandler(
    filename = function(){
      paste("limma_plot", input$select_image_type_limma, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_limma()[[5]], width = input$Boxplot_limma_output_width, height = input$Boxplot_limma_output_height, dpi = input$Boxplot_limma_output_dpi, units = "in")
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
    run_with_tracking("action_edger", function() {
      source("scripts/edger.R")
      labels_data_type<- input$file1$name
      edger_summary(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_edger_pvalue, group1 = input$group1_edger, group2 = input$group2_edger, plot_method = input$select_edger_plot, alpha = input$edger_pvalue, edger_color_palette = input$edger_color_palette)
    })
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
  
 
  output$boxplot_edger <- renderPlot({
    data_edger()[5]
  })
  
  
  output$download_Boxplot_edger<- downloadHandler(
    filename = function(){
      paste("edger_plot", input$select_image_type_edger, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_edger()[[5]], width = input$Boxplot_edger_output_width, height = input$Boxplot_edger_output_height, dpi = input$Boxplot_edger_output_dpi, units = "in")
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
    run_with_tracking("action_kruskal_wallis_test", function() {
      source("scripts/kruskal_wallis_test.R")
      labels_data_type<- input$file1$name
      kruskal_wallis_test_summary(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_kruskal_wallis_test_pvalue, ad_hoc = input$kruskal_wallis_test_ad_hoc, plot_method = input$kruskal_wallis_test_plot, alpha = input$kruskal_wallis_test_pvalue, kruskal_wallis_test_color_palette = input$kruskal_wallis_test_color_palette)
    })
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
  
   output$boxplot_kruskal_wallis_test <- renderPlot({
    data_kruskal_wallis_test()[6]
  })
  
  
  output$download_Boxplot_kruskal_wallis_test<- downloadHandler(
    filename = function(){
      paste("kruskal_wallis_test_plot", input$select_image_type_kruskal_wallis_test, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_kruskal_wallis_test()[[6]], width = input$Boxplot_kruskal_wallis_test_output_width, height = input$Boxplot_kruskal_wallis_test_output_height, dpi = input$Boxplot_kruskal_wallis_test_output_dpi, units = "in")
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
    run_with_tracking("action_anova", function() {
      source("scripts/anova.R")
      labels_data_type<- input$file1$name
      anova_summary(OTU_input = dataInput_RA_level()[[4]], group_index = dataInput_RA_level()[[9]], index_pvalue = input$select_anova_pvalue, ad_hoc = input$anova_ad_hoc, plot_method = input$anova_plot, alpha = input$anova_pvalue, anova_color_palette = input$anova_color_palette)
    })
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
  
  
  output$boxplot_anova <- renderPlot({
    data_anova()[6]
  })
  
  
  output$download_Boxplot_anova<- downloadHandler(
    filename = function(){
      paste("anova_plot", input$select_image_type_anova, sep="")
    },
    content = function(file){
      ggsave(file,plot = data_anova()[[6]], width = input$Boxplot_anova_output_width, height = input$Boxplot_anova_output_height, dpi = input$Boxplot_anova_output_dpi, units = "in")
    }
  )    
  
  
  
 
  
}
