library(shiny)
library(shinythemes)
library(shinyFiles)
library(DT)
#library(bslib)
#options(warn=-1)
shinythemes::themeSelector()
shinyUI(
    navbarPage(
    theme = shinytheme("cerulean"),
    "",
    tabPanel(
      "MetaDAVis",
      mainPanel(
        h1("MetaDAVis",align = "center"),
        hr(),
        h3("Introduction"),
        hr(),
        p(strong("MetaDAVis"), " (interactive Metagenome Data Analysis and Visualization) is a browser-based and user-friendly R Shiny application for researchers without programming proficiency to analyze and visualize metagenomics results from kingdom to species level. It comprises six functional analyses."),
        
        HTML("<B>The package includes the following:</B><br>
        <ul>
          <li>Data summary and abundance distribution</li>
          The data can be visualized in the stacked bar plot for abundance percentage, value, and relative frequency from 2 to 100 Taxa.<br>
            <ul>
            <li>Group: Samples are grouped by given conditions based on the metadata. </li>
            <li>Individual: Sample-based plots.</li>
            </ul>
          <li>Diversity analysis</li>
            <ul>
            <li>Alpha: Seven different methods was used from the phyloseq package. The results were displayed in a box and violin plot with a summary table.</li>
            <li>Beta: A total of 42 different diversity metrics were integrated from the phyloseq (unlist(distanceMethodList)) package with six selection methods. Results are visualized by bar and ordination with a summary table.</li>
            </ul>
          <li>Dimension reduction</li>
            <ul>
            <li>PCA: The ggfortify and plotly package was used to display plots in 2D(with and without labels and frame) and 3D with their summary table.</li>
            <li>t-SNE: Six different methods was used from the scater package. The samples were displayed in a t-SNE plot (2 and 3 dimensions) with a summary table.</li>
            <li>UMAP: Six different methods was used from the scater package. Displays the UMAP plot in sample and cluster-based with a summary table.</li>
            </ul>
          <li>Correlation analysis</li>
            <ul>
            <li>Taxa-based: The ggfortify package was used to display plot (with and without labels and frame) and their summary table.</li>
            <li>Sample-based: Six different methods was used from the scater package. The result was displayed in 2 and 3 dimensions t-SNE plots with a summary table.</li>
            </ul>
          <li>Heatmap: It was integrated with the ComplexHeatmap package. Display heatmap with and without row and column dendrograms and names.</li>
          <li>Differential abundance</li>
            <ul>
            <li>Two groups: Six different analyses were provided using the Wilcoxon Rank Sum test, t-test, metagenomeSeq, DESeq2, Limma-Voom, and edgeR. These will perform statistical analysis and generate plots and summary tables based on the significant taxa.</li>
            <li>Multiple groups comparison: Two different analyses, such as Kruskal-Wallis test and ANOVA was used for more than multiple group comparisons. These will perform statistical analysis and generate plots and summary tables based on the significant taxa.</li>
            </ul>
          </ul>
        <p>*It provides publication quality plots in seven formats: JPG, TIFF, PDF, SVG, BMP, EPS, and PS and summary tables (.csv format) to visualize and download.</p>
        <hr>
        <h3> use MetaDAVis online</h3>
        <p>MetaDAVis is deployed at: <a href='https://www.gudalab-rtools.net/MetaDAVis'>https://www.gudalab-rtools.net/MetaDAVis</a></p>
        <hr>
        <h3> Launch MetaDAVis using R and GitHub</h3>
        <p> MetaDAVis were deposited under the GitHub repository: <br>
        Before running the app, the user must have R (>= 4.4.1), RStudio (>= 2024.09.0), Bioconductor (>= 3.19) and Shiny (>= 1.9.1) (Tested with this version).<br>
         If users use an older R version, they may encounter errors in installing packages, So the users are recommended to update their R version first.<br>
         Once the user opens the R in the command line or Rstudio, need to run the following command in R to install the shiny package.<br><br></p>
          
<pre>install.packages('shiny')<br>
library(shiny)</pre>
          <hr>
          <h3>Start the app</h3>
          Start the R session using RStudio and run these lines:<br><br>
<pre>shiny::runGitHub('MetaDAVis','GudaLab')</pre>
or
Alternatively, download the source code from GitHub and run the following command in the R session using RStudio:
<pre>
library(shiny)
runApp('/path/to/the/MetaDAVis-master', launch.browser=TRUE)</pre>
<hr>
<h3>Help manual for the usage of MetaDAVis <a href='manual/MetaDAVis_manual.pdf', target='_blank'>[Download]</a></h3>
<hr>
<h3> Developed and maintained by</h3>
<p>This application was created by Sankarasubramanian Jagdesan and Babu Guda.  We share the passion about developing an user-friendly tool for all biologists, especially those who do not have access to bioinformaticians or programming efficency.
</p>
<hr>
"),
        
    ),
   ),

       tabPanel(
      "Upload files",
      sidebarLayout(
        sidebarPanel(
          h3("Upload files"),
          selectInput("select_file_type", label = "Select Input format", choices = list("Qiime2" = "qiime_format", "Megan" ="Megan", "Taxa count file (prepare your own file based on examble)" = "check"), selected = "qiime_format"),
          h5("The file accepts .txt or .tsv (Megan and users own file) or .csv formats (Qiime2)"),
           fluidRow(
          column(width = 8, fileInput("file1", "Upload count file ", accept = c(".tsv", ".txt", ".csv"),multiple = F)),
          column(width = 4, selectInput("file1_split", "Fields separated by", c("tab" = "\t", "Comma" =","), selected = "tab")),
          ),
          fluidRow(
          column(width = 8, fileInput("metadata", "Upload meta-data ", accept = c(".tsv", ".txt", ".csv"),multiple = F)),
          column(width = 4, selectInput("metadata_split", "Fields separated by", c("tab" = "\t", "Comma" =","), selected = "tab")),
          ),
          radioButtons("select_RA_type", "Choose the level to display", choices = c("Kingdom" = 1, "Phylum" = 2, "Class" = 3, "Order" = 4, "Family" = 5,  "Genus" = 6, "Species" = 7), selected = 4),
          actionButton("action_level", "Submit")
          
        ),
        mainPanel(
          tabsetPanel(
            type = "tabs",
            tabPanel(
              "Summary",    
         h4("Download example data"),
         tags$b("Qiime2 format"),
         br(),
         a(href="example_data/Qiime2_Greengenes_output_level7.csv", "Qiime2 Greengenes Output format",download=NA, target="_blank"),
         br(),
         a(href="example_data/Qiime2_metadata_for_Greengenes.csv", "Qiime2 metadata for greengenes",download=NA, target="_blank"),
         br(),
         a(href="example_data/Qiime2_Silva_output_level7.csv", "Qiime2 Silva Outpus format",download=NA, target="_blank"),
         br(),
         a(href="example_data/Qiime2_metadata_for_Silva.csv", "Qiime2 metadata for Silva",download=NA, target="_blank"),
         br(),
         tags$b("MEGAN output format"),
         br(),
         a(href="example_data/Megan_WGS_output.tsv", "Megan WGS output format",download=NA, target="_blank"),
         br(),
         a(href="example_data/Megan_WGS_metadata.tsv", "Megan metadata",download=NA, target="_blank"),
         br(),
         tags$b("Taxa count file (Prepare your input accourding to our count and metadata format)",download=NA, target="_blank"),
         br(),
         a(href="example_data/Taxa_count_file.tsv", "Count files format",download=NA, target="_blank"),
         br(),
         a(href="example_data/Tax_metadata.tsv", "Metadata",download=NA, target="_blank"),
         #downloadButton("Qiime2_metadata_for_Silva.csv", "Metadata"),
         hr(),
         h4("After the data is uploaded and checked, it will be displayed in the table summary below."),
          hr(),
          h4("Number of OTUs"),
         withSpinner(verbatimTextOutput("text_level")),
          h4("Metadata"),
         withSpinner(verbatimTextOutput("text_metadata"))
            ),
          tabPanel(
            "Taxonomy table",
            h4("Display the taxonomy counts for each samples"),
            withSpinner(dataTableOutput("taxonomy_table")),
           downloadButton(outputId = "download_taxonomy_table", label = "Download as csv"),
          ),
            tabPanel(
              "Metadata table",
              h4("Display the metadata file"),
              fluidRow(
                column(
              withSpinner(dataTableOutput("metadata_table")),
              downloadButton(outputId = "download_metadata_table", label = "Download as csv"), width = 4),
          ),
        ),
          tabPanel(
            "No. of conditions",
            h4("Display the number of condition based on your metadata"),
              fluidRow(
              column(
                withSpinner(dataTableOutput("conditions_table")),
            downloadButton(outputId = "download_conditions_table", label = "Download as csv"), width = 4),
              ),
          ),
          tabPanel(
            "Counts in samples",
            h4("Display the total number of counts in each samples"),
            fluidRow(
              column(
                withSpinner(dataTableOutput("count_table")),
            downloadButton(outputId = "download_count_table", label = "Download as csv"), width = 4),
            ),
          ),
        ),
      ),
    ),
    ),
    #distribution
    navbarMenu(
      "Distribution",
      tabPanel(
        "Group",
        sidebarLayout(
          sidebarPanel(
            h3("Distribution of top bacterial taxa (groups)"), selectInput("input_RA_bar_plot_group", label = "Selected input", choices = "No data selected! please load the data first"),
            radioButtons("select_plot_type_group", "Types of plot", choices = c("Abundance (%) - stacked bar" = 1,"Abundance value - stacked bar" = 2, "Relative frequency - stacked bar" = 3), selected = 1),
            numericInput("top_n_bar_plot_group", label = "Number of top bacterial taxa (Max = 100)", value = 15),
            selectInput("select_image_type_group", label = "Output image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
            actionButton("action_m1_bar_plot_group", "Submit"),
          ),
          mainPanel(
            withSpinner(plotOutput("bar_plot_group", width = "50%", height = "500px")),
            h4("Relative abundance."),
            fluidRow(
            column(3, numericInput("bar_plot_group_output_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px")),
            column(3, numericInput("bar_plot_group_output_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px")),
            column(3, numericInput("bar_plot_group_output_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px")),
            ),
            downloadButton(outputId = "download_bar_plot_group", label = "Download plot"),
          ),
        ),
      ),
      tabPanel(
        "Individual",
        sidebarLayout(
          sidebarPanel(
            h3("Distribution of top bacterial taxa (samples)"),
            selectInput("input_RA_bar_plot_individual", label = "Selected input", choices = "No data selected! please load the data first"),
            radioButtons("select_plot_type_individual", "Types of plot", choices = c("Abundance (%) - stacked bar" = 1,"Abundance value - stacked bar" = 2, "Relative frequency - stacked bar" = 3), selected = 1),
            numericInput("top_n_bar_plot_individual", label = "Number of top bacterial taxa (Max = 100)", value = 15),
            selectInput("select_image_type_individual", label = "Output image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
            actionButton("action_m1_bar_plot_individual", "Submit"),
          ),
          mainPanel(
            withSpinner(plotOutput("bar_plot_individual")),
            fluidRow(
            column(3, numericInput("bar_plot_individual_output_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px")),
            column(3, numericInput("bar_plot_individual_output_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px")),
            column(3, numericInput("bar_plot_individual_output_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px")),
            ),
            downloadButton(outputId = "download_bar_plot_individual", label = "Download plot")
          ),
        ),
      ),
    ),   
    
    #Diversity
    navbarMenu("Diversity",
               tabPanel(
                 "Alpha",
                 sidebarLayout(
                   sidebarPanel(
                     h3("Alpha diversity"),
                     selectInput("input_RA_Alpha", label = "Selected input",choices = "No data selected! please load the data first"),
                     selectInput("select_alpha", label = "Select Method", choices = list("Observed" = 1, "Chao1" = 2, "ACE" = 3, "Shannon" = 4, "Simpson" = 5, "InvSimpson" = 6, "Fisher"=7, "All_Combined"=8), selected = 1),
                     radioButtons("select_alpha_pvalue", label = "Wilcoxon test", choices = c("Yes (show's Pvalue)" = "Yes", "No" = "No", "Show *" = "star"), selected = "No"),
                     radioButtons("select_plot",label = "Types of plot", choices = c("Box plot" = 1, "Violin plot" = 2), selected = 1),
                     selectInput("select_image_type_alpha", label = "Output image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
                     actionButton("action_alpha_diversity", "Submit"),
                   ),
                   mainPanel(
                     tabsetPanel(
                       type = "tabs",
                       tabPanel(
                         "Alpha diversity plot",
                         h3("Boxplot"),
                         withSpinner(plotOutput("boxplot_Alpha_Div")),
                         br(),
                         br(),
                         br(),
                         br(),
                         fluidRow(
                         column(3, numericInput("Boxplot_alpha_div_output_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px")),
                         column(3, numericInput("Boxplot_alpha_div_output_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px")),
                         column(3, numericInput("Boxplot_alpha_div_output_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px")),
                         ),
                         downloadButton(outputId = "download_Boxplot_Alpha_Div", label = "Download plot"),
                       ),
                       tabPanel(
                         "Summary Table",
                         h3("Result - alpha diversity estimates for each metagenome"),
                         hr(),
                         withSpinner(dataTableOutput("alpha_table")),
                         downloadButton(outputId = "download_result_alpha", label = "Download as csv"),
                       ),
                     ),
                   ),
                 ),
               ),
               
               tabPanel(
                 "Beta",
                 sidebarLayout(
                   sidebarPanel(
                     h3("Beta diversity"),
                     selectInput("input_RA_Beta", label = "Selected input", choices = "No data selected! please load the data first"),
                     selectInput("select_beta", label = "Select diversity metrics", choices = c("bray-curtis" = "bray", "jaccard" = "jaccard", "manhattan" = "manhattan", "euclidean" = "euclidean", "canberra" = "canberra", "kulczynski" = "kulczynski", "gower" = "gower", "altGower" = "altGower", "morisita" = "morisita", "horn" = "horn", "mountford" = "mountford", "raup" = "raup", "binomial" = "binomial", "chao" = "chao", "cao" ="cao", "w" = "w","-1" = "-1","c" = "c","wb" = "wb","r" = "r","I" = "I","e" = "e","t" = "t","me" = "me","j" = "j","sor" = "sor","m" = "m","-2" = "-2","co" = "co","cc" = "cc","g" = "g","-3" = "-3","l" = "l","19" = "19","hk" = "hk","rlb" = "rlb","sim" = "sim","gl" = "gl","z" = "z","maximum" = "maximum","binary" = "binary","minkowski" = "minkowski"), selected = "bray"),
                     selectInput("select_method", "Select method", choices = c("PCoA" = "PCoA", "NMDS" = "NMDS", "DCA" = "DCA", "CCA" = "CCA", "RDA" = "RDA", "MDS" = "MDS"), selected = "PCoA"),
                     selectInput("select_image_type_beta", label = "Output image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
                     actionButton("action_beta_diversity", "Submit"),
                   ),
                   mainPanel(
                     tabsetPanel(
                       type = "tabs",
                       tabPanel(
                         "Beta diversity Plot",
                         h3("Beta diversity plot"),
                         withSpinner(plotOutput("boxplot_Beta_Div")),
                         fluidRow(
                         column(3, numericInput("Boxplot_beta_div_output_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px")),
                         column(3, numericInput("Boxplot_beta_div_output_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px")),
                         column(3, numericInput("Boxplot_beta_div_output_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px")),
                         ),
                         #numericInput("Boxplot_beta_div_output_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 200, width = "300px"),
                         downloadButton(outputId = "download_Boxplot_beta_Div", label = "Download plot"),
                       ),
                       tabPanel(
                         "Summary Table",
                         h3("Result - distance between all the samples"),
                         hr(),
                         withSpinner(dataTableOutput("beta_table")),
                         downloadButton(outputId = "download_result_beta", label = "Download as csv"),
                       ),
                   ),
                 ),
               ),


),
),
navbarMenu("Dimension reduction",
           tabPanel(
             "PCA-2D",
             sidebarLayout(
               sidebarPanel(
                 h3("PCA-2D"),
                 selectInput("input_RA_pca", label = "Selected input", choices = "No data selected! please load the data first"),
                 selectInput("select_pca_label", "Label", choices = c("TRUE"="TRUE", "FALSE"="FALSE"), selected = FALSE),
                 numericInput("select_pca_label_size", label = "Label size", value = 3),
                 selectInput("select_pca_frame", "Frame", choices = c("TRUE"="TRUE", "FALSE"="FALSE"), selected = FALSE),
                 selectInput("select_image_type_pca", label = "Output image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
                 actionButton("action_pca", "Submit"),
               ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel(
                     "PCA 2D Plot",
                     h3("Principal Component Analysis (PCA)"),
                     withSpinner(plotOutput("plot_pca", width = "70%", height = "500px")),
                     fluidRow(
                       column(3, numericInput("pca_plot_output_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("pca_plot_output_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("pca_plot_output_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px")),
                     ),
                     downloadButton(outputId = "download_plot_pca", label = "Download plot"),
                   ),
                   tabPanel(
                     "Summary Table",
                     fluidRow(
                       column(
                         withSpinner(dataTableOutput("pca_table")),
                     downloadButton(outputId = "download_result_pca", label = "Download as csv"), width = 6),
                     ),
                   ),
                 ),
               ),
             ),
           ),
           tabPanel(
             "PCA-3D",
             sidebarLayout(
               sidebarPanel(
                 h3("PCA-3D"),
                 selectInput("input_RA_pca3d", label = "Selected input", choices = "No data selected! please load the data first"),
                 #selectInput("select_image_type_pca3d", label = "Output image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
                 actionButton("action_pca3d", "Submit"),
               ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel(
                     "PCA 3D Plot",
                     h3("Principal Component Analysis (PCA)"),
                     withSpinner(plotlyOutput("plot_pca3d", height = "700px", width= "800px")),
                  
                   ),
                   tabPanel(
                     "Summary Table",
                     fluidRow(
                       column(
                         withSpinner(dataTableOutput("pca3d_table")),
                    downloadButton(outputId = "download_result_pca3d", label = "Download as csv"), width = 6),
                     ),
                   ),
                 ),
               ),
             ),
           ),
           tabPanel(
             "t-SNE",
             sidebarLayout(
               sidebarPanel(
                 h3("t-SNE"),
                 selectInput("input_tsne", label = "Selected input", choices = "No data selected! please load the data first"),
                 selectInput("select_tsne_method", "Select method", choices = c("counts"="counts", "rclr"="rclr", "hellinger"="hellinger", "pa"="pa", "rank"="rank", "relabundance"="relabundance"), selected = "counts"),
                 selectInput("select_tsne_dimension", "Select dimension to display", choices = c("2" = 2, "3" = 3), selected = 2),
                 selectInput("select_image_type_tsne", label = "Output image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
                 actionButton("action_tsne", "Submit"),
               ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel(
                     "t-SNE Plot",
                     h3("t-distributed Stochastic Neighbor Embedding (t-SNE)"),
                     withSpinner(plotOutput("plot_tsne", width = "70%", height = "500px")),
                     fluidRow(
                       column(3, numericInput("tsne_plot_output_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("tsne_plot_output_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("tsne_plot_output_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px")),
                     ),
                     downloadButton(outputId = "download_plot_tsne", label = "Download plot"),
                   ),
                   tabPanel(
                     "Summary Table",
                     fluidRow(
                       column(
                         withSpinner(dataTableOutput("tsne_table")),
                     downloadButton(outputId = "download_result_tsne", label = "Download as csv"), width = 6),
                 ),
               ),
                 ),
               ),
             ),
           ),
           tabPanel(
             "UMAP",
             sidebarLayout(
               sidebarPanel(
                 h3("UMAP"),
                 selectInput("input_umap", label = "Selected input", choices = "No data selected! please load the data first"),
                 selectInput("select_umap_method", "Select method", choices = c("counts"="counts", "rclr"="rclr", "hellinger"="hellinger", "pa"="pa", "rank"="rank", "relabundance"="relabundance"), selected = "counts"),
                 selectInput("select_umap_kvalue", "Select k value (for graph construction)", choices = c("2" = 2, "3" = 3, "4" = 4, "5" = 5, "6" = 6, "7" = 7, "8" = 8, "9" = 9, "10" = 10, "11" = 11, "12" = 12, "13" = 13, "14" = 14, "15" = 15), selected = 2),
                 selectInput("select_image_type_umap", label = "Output image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
                 actionButton("action_umap", "Submit"),
               ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel(
                     "UMAP Plot",
                     h3("Uniform Manifold Approximation and Projection for Dimension Reduction (UMAP)"),
                     withSpinner(plotOutput("plot_umap", width = "100%", height = "500px")),
                     fluidRow(
                       column(3, numericInput("umap_plot_output_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("umap_plot_output_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("umap_plot_output_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px")),
                     ),
                     downloadButton(outputId = "download_plot_umap", label = "Download plot"),
                   ),
                   tabPanel(
                     "Summary Table based on condition",
                     fluidRow(
                       column(
                         withSpinner(dataTableOutput("umap_table")),
                     downloadButton(outputId = "download_result_umap", label = "Download as csv"), width = 6),
                 ),
               ),
                   tabPanel(
                     "Summary Table based on cluster",
                     fluidRow(
                       column(
                         withSpinner(dataTableOutput("umap_table1")),
                     downloadButton(outputId = "download_result_umap1", label = "Download as csv"), width = 6),
             ),
           ),
                 ),
               ),
             ),
           ),
),
#correlation
navbarMenu(
  "Correlation",
  tabPanel(
    "Taxa-based",
    sidebarLayout(
      sidebarPanel(
        h3("Compute correlation between taxa"), selectInput("input_taxa_based_correlation", label = "Selected input", choices = "No data selected! please load the data first"),
        radioButtons("select_taxa_based_correlation_method", "Correlation methods", choices = c("pearson" = "pearson","kendall" = "kendall", "spearman" = "spearman"), selected = "pearson"),
        numericInput("select_taxa_based_correlation_label_size", label = "Label size", value = 3),
        selectInput("select_image_type_taxa_based_correlation", label = "Output image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
        actionButton("action_taxa_based_correlation", "Submit"),
      ),
      mainPanel(
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Correlation plot",
            withSpinner(plotOutput("plot_taxa_based_correlation", width = "100%", height = "1000px")),
        h4("Taxa based correlation plot."),
        fluidRow(
          column(3, numericInput("taxa_based_correlation_output_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px")),
          column(3, numericInput("taxa_based_correlation_output_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px")),
          column(3, numericInput("taxa_based_correlation_output_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px")),
        ),
        downloadButton(outputId = "download_plot_taxa_based_correlation", label = "Download plot"),
      ),
      tabPanel(
        "Summary Table",
        fluidRow(
          column(
            withSpinner(dataTableOutput("taxa_based_correlation_table")),
        downloadButton(outputId = "download_result_taxa_based_correlation", label = "Download as csv"), width = 8),
        ),
      ),
     ),
    ),
  ),
  ),
  tabPanel(
    "Sample-based",
    sidebarLayout(
      sidebarPanel(
        h3("Compute correlation between samples"), selectInput("input_samples_based_correlation", label = "Selected input", choices = "No data selected! please load the data first"),
        radioButtons("select_samples_based_correlation_method", "Correlation methods", choices = c("pearson" = "pearson","kendall" = "kendall", "spearman" = "spearman"), selected = "pearson"),
        numericInput("select_samples_based_correlation_label_size", label = "Label size", value = 3),
        selectInput("select_image_type_samples_based_correlation", label = "Output image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
        actionButton("action_samples_based_correlation", "Submit"),
      ),
      mainPanel(
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Correlation plot",
            withSpinner(plotOutput("plot_samples_based_correlation", width = "100%", height = "1000px")),
            h4("Samples based correlation plot."),
            fluidRow(
              column(3, numericInput("samples_based_correlation_output_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px")),
              column(3, numericInput("samples_based_correlation_output_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px")),
              column(3, numericInput("samples_based_correlation_output_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px")),
            ),
            downloadButton(outputId = "download_plot_samples_based_correlation", label = "Download plot"),
          ),
          tabPanel(
            "Summary Table",
            fluidRow(
              column(
                withSpinner(dataTableOutput("samples_based_correlation_table")),
            downloadButton(outputId = "download_result_samples_based_correlation", label = "Download as csv"), width = 8),
        ),
      ),
        ),
      ),
    ),
  ),
),
#Heatmap
tabPanel(
  "Heatmap",
  sidebarLayout(
    sidebarPanel(
      h3("Heatmap - relative abundance"),
      selectInput("input_heatmap", label = "Selected input", choices = "No data selected! please load the data first"),
      selectInput("heatmap_row_names", "Show row names", choices = c("TRUE"="TRUE", "FALSE"="FALSE"), selected = TRUE),
      numericInput("heatmap_row_names_size", label = "Row name size", value = 7),
      selectInput("heatmap_column_names", "Show column names", choices = c("TRUE"="TRUE", "FALSE"="FALSE"), selected = TRUE),
      numericInput("heatmap_column_names_size", label = "Column name size", value = 7),
      selectInput("heatmap_row_dend", "Show row cladogram", choices = c("TRUE"="TRUE", "FALSE"="FALSE"), selected = TRUE),
      selectInput("heatmap_column_dend", "Show column cladogram", choices = c("TRUE"="TRUE", "FALSE"="FALSE"), selected = TRUE),
      
      selectInput("select_image_type_heatmap", label = "Output image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      actionButton("action_heatmap", "Submit"),
    ),
    mainPanel(
      withSpinner(plotOutput("plot_heatmap", width = "100%", height = "800px")),
      h4("Heatmap using relative abundance."),
      fluidRow(
        column(3, numericInput("heatmap_output_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px")),
        column(3, numericInput("heatmap_output_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px")),
        column(3, numericInput("heatmap_output_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px")),
      ),
      downloadButton(outputId = "download_heatmap", label = "Download Heatmap"),
    ),
  ),
),

navbarMenu("Differential abundance",
           "Two groups",
           tabPanel(
             "Wilcoxon Rank Sum test",
             sidebarLayout(
               sidebarPanel(
                 h3("Wilcoxon Rank Sum test"),
                 selectInput("input_wilcoxtest", label = "Selected input", choices = "No data selected! please load the data first"),
                 selectInput("group1_wilcoxtest", label = "Select condition1", choices = "Please upload metadata in upload page"),
                 selectInput("group2_wilcoxtest", label = "Select condition2", choices = "Please upload metadata in upload page"),
                 selectInput("select_wilcoxtest_pvalue", label = "Test correction", choices = list("Benjamini-Hochberg FDR" = "padj", "P-value" = "pvalue"), selected = "padj"),
                 numericInput("wilcoxtest_pvalue","FDR or Pvalue", value = "0.05"),
                 radioButtons("select_wilcoxtest_plot", "Types of plot", choices = c("Grouped box plot" = 1, "Individual box plot" = 2, "Volcano plot" = 3, "Heatmap" = 4), selected = 1),
                 selectInput("select_image_type_wilcoxtest", label = "Output image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
                 actionButton("action_wilcoxtest", "Submit"),
               ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel(
                     "Summary Table",
                     h3("Result - OTUs that were significantly different between two groups"),
                     withSpinner(verbatimTextOutput("text_wilcoxtest_level")),
                     fluidRow(
                       column(
                         withSpinner(dataTableOutput("wilcoxtest_table")), width = 12),
                 ),
                     hr(),
                     fluidRow(
                     column(3, h4("Download significant"),
                     downloadButton(outputId = "download_result_wilcoxtest_1", label = "Download as csv")),
                     column(3, h4("Download all"),
                     downloadButton(outputId = "download_result_wilcoxtest_2", label = "Download as csv")),
                     column(3,h4("Download relative frequency"),
                     downloadButton(outputId = "download_result_wilcoxtest_3", label = "Download as csv")),
                     column(3,h4("Total counts in each samples"),
                     downloadButton(outputId = "download_result_wilcoxtest_4", label = "Download as csv")),
                     ),
                   ),
                   tabPanel(
                     "Plot",
                     h3("Plot"),
                     withSpinner(plotOutput("boxplot_wilcoxtest", width = "100%", height = "800px")),
                     fluidRow(
                       column(3, numericInput("Boxplot_wilcoxtest_output_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("Boxplot_wilcoxtest_output_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("Boxplot_wilcoxtest_output_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px")),
                     ),
                     downloadButton(outputId = "download_Boxplot_wilcoxtest", label = "Download plot"),
                   ),
                 ),
               ),
             ),
           ),
           tabPanel(
             "t-test",
             sidebarLayout(
               sidebarPanel(
                 h3("t-test: Two sample t-test"),
                 selectInput("input_ttest", label = "Selected input", choices = "No data selected! please load the data first"),
                 selectInput("group1_ttest", label = "Select condition1", choices = "Please upload metadata in upload page"),
                 selectInput("group2_ttest", label = "Select condition2", choices = "Please upload metadata in upload page"),
                 selectInput("select_ttest_pvalue", label = "Test correction", choices = list("Benjamini-Hochberg FDR" = "padj", "P-value" = "pvalue"), selected = "padj"),
                 numericInput("ttest_pvalue","FDR or Pvalue", value = "0.05"),
                 radioButtons("select_ttest_plot", "Types of plot", choices = c("Grouped box plot" = 1, "Individual box plot" = 2, "Volcano plot" = 3, "Heatmap" = 4), selected = 1),
                 selectInput("select_image_type_ttest", label = "Output image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
                 actionButton("action_ttest", "Submit"),
               ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                    tabPanel(
                     "Summary Table",
                     h3("Result - OTUs that were significantly different between two groups"),
                     withSpinner(verbatimTextOutput("text_ttest_level")),
                      fluidRow(
                       column(
                         withSpinner(dataTableOutput("ttest_table")), width = 12),
                     ),
                     hr(),
                     fluidRow(
                     column(3, h4("Download significant"),
                     downloadButton(outputId = "download_result_ttest_1", label = "Download as csv")),
                     column(3, h4("Download all"),
                     downloadButton(outputId = "download_result_ttest_2", label = "Download as csv")),
                     column(3, h4("Download relative frequency"),
                     downloadButton(outputId = "download_result_ttest_3", label = "Download as csv")),
                     column(3, h4("Total counts in each samples"),
                     downloadButton(outputId = "download_result_ttest_4", label = "Download as csv")),
                     ),
                   ),
                   tabPanel(
                     "Plot",
                     h3("Plot"),
                     withSpinner(plotOutput("boxplot_ttest", width = "100%", height = "800px")),
                     fluidRow(
                       column(3, numericInput("Boxplot_ttest_output_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("Boxplot_ttest_output_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("Boxplot_ttest_output_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px")),
                     ),
                     downloadButton(outputId = "download_Boxplot_ttest", label = "Download plot"),
                   ),
                 ),
               ),
             ),
           ),
           tabPanel(
             "metagenomeSeq",
             sidebarLayout(
               sidebarPanel(
                 h3("metagenomeSeq"),
                 selectInput("input_RA_metagenomeseq", label = "Selected input", choices = "No data selected! please load the data first"),
                 selectInput("group1_metagenomeseq", label = "Select condition1", choices = "Please upload metadata in upload page"),
                 selectInput("group2_metagenomeseq", label = "Select condition2", choices = "Please upload metadata in upload page"),
                 selectInput("select_metagenomeseq_pvalue", label = "Test correction", choices = list("Benjamini-Hochberg FDR" = "padj", "P-value" = "pvalue"), selected = "padj"),
                 numericInput("metagenomeseq_pvalue","FDR or Pvalue", value = "0.05"),
                 radioButtons("select_metagenomeseq_plot", "Types of plot", choices = c("Grouped box plot" = 1, "Individual box plot" = 2, "Volcano plot" = 3, "Heatmap" = 4), selected = 1),
                 selectInput("select_image_type_metagenomeseq", label = "Output image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
                 actionButton("action_metagenomeseq", "Submit"),
               ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel(
                     "Summary Table",
                     h3("Result - OTUs that were significantly different between two groups"),
                     withSpinner(verbatimTextOutput("text_metagenomeseq_level")),
                     fluidRow(
                       column(
                         withSpinner(dataTableOutput("metagenomeseq_table")), width = 12),
                     ),
                     hr(),
                     fluidRow(
                     column(3, h4("Download significant"),
                     downloadButton(outputId = "download_result_metagenomeseq_1", label = "Download as csv")),
                     column(3, h4("Download all"),
                     downloadButton(outputId = "download_result_metagenomeseq_2", label = "Download as csv")),
                     column(3, h4("Total counts in each samples"),
                     downloadButton(outputId = "download_result_metagenomeseq_3", label = "Download as csv")),
                     ),
                   ),
                   tabPanel(
                     "Plot",
                     h3("Plot"),
                     withSpinner(plotOutput("boxplot_metagenomeseq", width = "100%", height = "800px")),
                     fluidRow(
                       column(3, numericInput("Boxplot_metagenomeseq_output_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("Boxplot_metagenomeseq_output_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("Boxplot_metagenomeseq_output_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px")),
                     ),
                     downloadButton(outputId = "download_Boxplot_metagenomeseq", label = "Download plot"),
                   ),
                   
                 ),
               ),
             ),
           ),
          
           tabPanel(
             "DESeq2",
             sidebarLayout(
               sidebarPanel(
                 h3("DESeq2"),
                 selectInput("input_RA_deseq2", label = "Selected input", choices = "No data selected! please load the data first"),
                 selectInput("group1_deseq2", label = "Select condition1", choices = "Please upload metadata in upload page"),
                 selectInput("group2_deseq2", label = "Select condition2", choices = "Please upload metadata in upload page"),
                 selectInput("select_deseq2_pvalue", label = "Test correction", choices = list("Benjamini-Hochberg FDR" = "padj", "P-value" = "pvalue"), selected = "padj"),
                 numericInput("deseq2_pvalue","FDR or Pvalue", value = "0.05"),
                 radioButtons("select_deseq2_plot", "Types of plot", choices = c("Grouped box plot" = 1, "Individual box plot" = 2, "Volcano plot" = 3, "Heatmap" = 4), selected = 1),
                 selectInput("select_image_type_deseq2", label = "Output image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
                 actionButton("action_deseq2", "Submit"),
               ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel(
                     "Summary Table",
                     h3("Result - OTUs that were significantly different between two groups"),
                     withSpinner(verbatimTextOutput("text_deseq2_level")),
                     fluidRow(
                       column(
                         withSpinner(dataTableOutput("deseq2_table")), width = 12),
                     ),
                     hr(),
                     fluidRow(
                     column(3, h4("Download significant"),
                     downloadButton(outputId = "download_result_deseq2_1", label = "Download as csv")),
                     column(3, h4("Download all"),
                     downloadButton(outputId = "download_result_deseq2_2", label = "Download as csv")),
                     column(3, h4("Download normalized count"),
                     downloadButton(outputId = "download_result_deseq2_3", label = "Download as csv")),
                     column(3, h4("Total counts in each samples"),
                     downloadButton(outputId = "download_result_deseq2_4", label = "Download as csv")),
                     ),
                   ),
                   tabPanel(
                     "Plot",
                     h3("Plot"),
                     withSpinner(plotOutput("boxplot_deseq2", width = "100%", height = "800px")),
                     fluidRow(
                       column(3, numericInput("Boxplot_deseq2_output_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("Boxplot_deseq2_output_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("Boxplot_deseq2_output_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px")),
                     ),
                     downloadButton(outputId = "download_Boxplot_deseq2", label = "Download plot"),
                   ),
                   
                 ),
               ),
             ),
           ),
           tabPanel(
             "Limma-Voom",
             sidebarLayout(
               sidebarPanel(
                 h3("Limma-Voom"),
                 selectInput("input_RA_limma", label = "Selected input", choices = "No data selected! please load the data first"),
                 selectInput("group1_limma", label = "Select condition1", choices = "Please upload metadata in upload page"),
                 selectInput("group2_limma", label = "Select condition2", choices = "Please upload metadata in upload page"),
                 selectInput("select_limma_pvalue", label = "Test correction", choices = list("Benjamini-Hochberg FDR" = "padj", "P-value" = "pvalue"), selected = "padj"),
                 numericInput("limma_pvalue","FDR or Pvalue", value = "0.05"),
                 radioButtons("select_limma_plot", "Types of plot", choices = c("Grouped box plot" = 1, "Individual box plot" = 2, "Volcano plot" = 3, "Heatmap" = 4), selected = 1),
                 selectInput("select_image_type_limma", label = "Output image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
                 actionButton("action_limma", "Submit"),
               ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel(
                     "Summary Table",
                     h3("Result - OTUs that were significantly different between two groups"),
                     withSpinner(verbatimTextOutput("text_limma_level")),
                     fluidRow(
                       column(
                         withSpinner(dataTableOutput("limma_table")), width = 12),
                     ),
                     hr(),
                     fluidRow(
                     column(3, h4("Download significant"),
                     downloadButton(outputId = "download_result_limma_1", label = "Download as csv")),
                     column(3, h4("Download all"),
                     downloadButton(outputId = "download_result_limma_2", label = "Download as csv")),
                     column(3, h4("Total counts in each samples"),
                     downloadButton(outputId = "download_result_limma_3", label = "Download as csv")),
                     ),
                    ),
                   tabPanel(
                     "Plot",
                     h3("Plot"),
                     withSpinner(plotOutput("boxplot_limma", width = "100%", height = "800px")),
                     fluidRow(
                       column(3, numericInput("Boxplot_limma_output_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("Boxplot_limma_output_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("Boxplot_limma_output_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px")),
                     ),
                     downloadButton(outputId = "download_Boxplot_limma", label = "Download plot"),
                   ),
                 ),
               ),
             ),
           ),
           
           tabPanel(
             "edgeR",
             sidebarLayout(
               sidebarPanel(
                 h3("edgeR"),
                 selectInput("input_RA_edger", label = "Selected input", choices = "No data selected! please load the data first"),
                 selectInput("group1_edger", label = "Select condition1", choices = "Please upload metadata in upload page"),
                 selectInput("group2_edger", label = "Select condition2", choices = "Please upload metadata in upload page"),
                 selectInput("select_edger_pvalue", label = "Test correction", choices = list("Benjamini-Hochberg FDR" = "padj", "P-value" = "pvalue"), selected = "padj"),
                 numericInput("edger_pvalue","FDR or Pvalue", value = "0.05"),
                 radioButtons("select_edger_plot", "Types of plot", choices = c("Grouped box plot" = 1, "Individual box plot" = 2, "Heatmap" = 3), selected = 1),
                 selectInput("select_image_type_edger", label = "Output image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
                 actionButton("action_edger", "Submit"),
               ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel(
                     "Summary Table",
                     h3("Result - OTUs that were significantly different between two groups"),
                     withSpinner(verbatimTextOutput("text_edger_level")),
                     fluidRow(
                       column(
                         withSpinner(dataTableOutput("edger_table")), width = 12),
                     ),
                     hr(),
                     fluidRow(
                     column(3, h4("Download significant"),
                     downloadButton(outputId = "download_result_edger_1", label = "Download as csv")),
                     column(3, h4("Download all"),
                     downloadButton(outputId = "download_result_edger_2", label = "Download as csv")),
                     column(3, h4("Total counts in each samples"),
                     downloadButton(outputId = "download_result_edger_3", label = "Download as csv")),
                     ),
                     ),
                   tabPanel(
                     "Plot",
                     h3("Plot"),
                     withSpinner(plotOutput("boxplot_edger", width = "100%", height = "800px")),
                     fluidRow(
                       column(3, numericInput("Boxplot_edger_output_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("Boxplot_edger_output_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("Boxplot_edger_output_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px")),
                     ),
                     downloadButton(outputId = "download_Boxplot_edger", label = "Download plot"),
                   ),
                 ),
               ),
             ),
           ),
           "----",
           "Multiple groups",
           tabPanel(
             "Kruskal-Wallis test",
             sidebarLayout(
               sidebarPanel(
                 h3("Kruskal-Wallis test"),
                 selectInput("input_kruskal_wallis_test", label = "Selected input", choices = "No data selected! please load the data first"),
                 selectInput("select_kruskal_wallis_test_pvalue", label = "Test correction", choices = list("Benjamini-Hochberg FDR" = "padj", "P-value" = "pvalue"), selected = "padj"),
                 numericInput("kruskal_wallis_test_pvalue","FDR or Pvalue", value = "0.05"),
                 radioButtons("kruskal_wallis_test_ad_hoc", "Post-hoc test", choices = c("Yes" = "Yes", "No" = "No"), selected = "No"),
                 radioButtons("kruskal_wallis_test_plot", "Types of plot", choices = c("Grouped box plot" = 1, "Individual box plot" = 2, "Heatmap" = 3), selected = 1),
                 selectInput("select_image_type_kruskal_wallis_test", label = "Output image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
                 actionButton("action_kruskal_wallis_test", "Submit"),
               ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                  tabPanel(
                     "Summary Table",
                     h3("Result - OTUs that were significantly different between multiple groups"),
                     withSpinner(verbatimTextOutput("text_kruskal_wallis_test_level")),
                     fluidRow(
                       column(
                         withSpinner(dataTableOutput("kruskal_wallis_test_table")), width = 12),
                     ),
                     hr(),
                     fluidRow(
                     column(3, h4("Download significant"),
                     downloadButton(outputId = "download_result_kruskal_wallis_test_1", label = "Download as csv")),
                     column(3, h4("Download all"),
                     downloadButton(outputId = "download_result_kruskal_wallis_test_2", label = "Download as csv")),
                     column(3, h4("Download relative frequency"),
                     downloadButton(outputId = "download_result_kruskal_wallis_test_3", label = "Download as csv")),
                     column(3, h4("Total counts in each samples"),
                     downloadButton(outputId = "download_result_kruskal_wallis_test_4", label = "Download as csv")),
                     ),
                   ),
                  tabPanel(
                    "Plot",
                    h3("Plot"),
                    withSpinner(plotOutput("boxplot_kruskal_wallis_test", width = "100%", height = "1500px")),
                    fluidRow(
                      column(3, numericInput("Boxplot_kruskal_wallis_test_output_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px")),
                      column(3, numericInput("Boxplot_kruskal_wallis_test_output_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px")),
                      column(3, numericInput("Boxplot_kruskal_wallis_test_output_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px")),
                    ),
                    downloadButton(outputId = "download_Boxplot_kruskal_wallis_test", label = "Download plot"),
                  ),
                 ),
               ),
             ),
           ),
           tabPanel(
             "ANOVA",
             sidebarLayout(
               sidebarPanel(
                 h3("Analysis of variance: ANOVA"),
                 selectInput("input_anova", label = "Selected input", choices = "No data selected! please load the data first"),
                 selectInput("select_anova_pvalue", label = "Test correction", choices = list("Benjamini-Hochberg FDR" = "padj", "P-value" = "pvalue"), selected = "padj"),
                 numericInput("anova_pvalue","FDR or Pvalue", value = "0.05"),
                 radioButtons("anova_ad_hoc", "Post-hoc test", choices = c("Yes" = "Yes", "No" = "No"), selected = "No"),
                 radioButtons("anova_plot", "Types of plot", choices = c("Grouped box plot" = 1, "Individual box plot" = 2, "Heatmap" = 3), selected = 1),
                 selectInput("select_image_type_anova", label = "Output image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
                 actionButton("action_anova", "Submit"),
               ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel(
                     "Summary Table",
                     h3("Result - OTUs that were significantly different between multiple groups"),
                     withSpinner(verbatimTextOutput("text_anova_level")),
                     fluidRow(
                       column(
                         withSpinner(dataTableOutput("anova_table")), width = 12),
                     ),
                     hr(),
                     fluidRow(
                     column(3, h4("Download significant"),
                     downloadButton(outputId = "download_result_anova_1", label = "Download as csv")),
                     column(3, h4("Download all"),
                     downloadButton(outputId = "download_result_anova_2", label = "Download as csv")),
                     column(3, h4("Download relative frequency"),
                     downloadButton(outputId = "download_result_anova_3", label = "Download as csv")),
                     column(3, h4("Total counts in each samples"),
                     downloadButton(outputId = "download_result_anova_4", label = "Download as csv")),
                     ),
                   ),
                   tabPanel(
                     "Plot",
                     h3("Plot"),
                     withSpinner(plotOutput("boxplot_anova", width = "100%", height = "1500px")),
                     fluidRow(
                       column(3, numericInput("Boxplot_anova_output_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("Boxplot_anova_output_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px")),
                       column(3, numericInput("Boxplot_anova_output_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px")),
                     ),
                     downloadButton(outputId = "download_Boxplot_anova", label = "Download plot"),
                   ),
                   
                 ),
               ),
             ),
           ),
),
)
)

