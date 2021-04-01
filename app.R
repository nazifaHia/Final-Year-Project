library(shiny)
library(shinydashboard)

library(scater)
library(plotly)
library(reshape2)
library(circlize)
library(sSeq)
library(shinycssloaders)
library(ComplexHeatmap)
library(tibble)
library(paletteer)
library(umap)
library(cluster)
library(dplyr)
library(tidyverse)
library(viridisLite)
library(HDF5Array)
library(formattable)
library(DT)
library(colourpicker)
library(shinyWidgets)
library(shinyjs)
library(ggpubr)
library(factoextra)
library(igraph)
library(scran)
library(dynamicTreeCut)
library(enrichR)
library(pheatmap)
library(Seurat)
library(rlist)



Sys.setenv(R_MAX_VSIZE = 16e9)

cdScFiltAnnot <- loadHDF5SummarizedExperiment(dir="cdScFiltAnnotHDF5", prefix="")

kmeans.res10 <- readRDS(file = "kmeans.res10.Rds")

list_cluster <- readRDS("Cluster.RDS")

cdScFiltAnnot$kmeansCluster <- kmeans.res10

c30 <- c("dodgerblue2",#1
         "#E31A1C", #2 red
         "green4", #3
         "#FF7F00", #4 orange
         "green1",#5
         "purple",#6
         "blue1",#7
         "deeppink1",#8
         "darkorange4",#9
         "black",#10
         "gold1",#11
         "darkturquoise",#12
         "#6A3D9A", #13 purple
         "orchid1",#14
         "gray70",#15
         "maroon",#16
         "palegreen2",#17
         "#333333",#18
         "#CAB2D6", #19 lt purple
         "#FDBF6F", #20 lt orange
         "khaki2",#21
         "skyblue2",#22
         "steelblue4",#23
         "green1",#24
         "yellow4",#25
         "yellow3",#26
         "#FB9A99", #27 lt pink
         "brown",#28
         "#000099",#29
         "#CC3300"#30
)



c_sample_col <- c30[c(1,3,23,19,30)]
c_clust_col <- c30[c(1,2,3,4,5,6,7,8,9,11,12,14,19,22,24,25)]

my.clusters <- cdScFiltAnnot$Clusters
#options(bitmapType='cairo')


df_shiny <<- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
df_shiny[,'Clusters'] <- as.factor(cdScFiltAnnot$Clusters)
df_shiny$Group <- "gray"
df_shiny$key <- row.names(df_shiny)


df_shiny_ForTwoClust <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
rownames(df_shiny_ForTwoClust) <- cdScFiltAnnot$Barcode
df_shiny_ForTwoClust$Group <- "gray88"
df_shiny_ForTwoClust$key <- row.names(df_shiny_ForTwoClust)


df_shiny_ForTwoSample <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
rownames(df_shiny_ForTwoSample) <- cdScFiltAnnot$Barcode
df_shiny_ForTwoSample$Group <- "gray88"
df_shiny_ForTwoSample$key <- row.names(df_shiny_ForTwoSample)

#cbPalette <- paletteer_c(package="pals", palette="glasbey", n=12)

colData(cdScFiltAnnot)[,'Clust_Sample'] <- paste0(colData(cdScFiltAnnot)$Clusters,'_',colData(cdScFiltAnnot)$Sample)

colData(cdScFiltAnnot)[,'Clust_SampleWithClust'] <- paste0(colData(cdScFiltAnnot)$Sample,'_',colData(cdScFiltAnnot)$Clusters)


GeneNameSorted <- sort(rownames(logcounts(cdScFiltAnnot)))
df = data.frame(Condition=colData(cdScFiltAnnot)$Sample)
#df <- df[order(df$Class),]
ha = HeatmapAnnotation(df = df)

flag <<- 0


ui <- dashboardPage(
  #skin = "purple",
  dashboardHeader(
    title = "Single cell RNA-Seq Clustering Analysis Tool"
  ),
  
  # Sidebar #############################
  dashboardSidebar(
    width = 230,
    
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Overview", tabName = "Overview", icon = icon("dashboard")),
     
      menuItem("Expression of Gene", tabName = "Gene_expressionAll", icon = icon("dna"),
               menuSubItem('Single gene expression', tabName = "Gene_expression")
               ),

      menuItem("Heatmap", tabName = "multiple_cluster_heatmap", icon = icon("columns"),
                 menuSubItem('Cluster heatmap', tabName = 'clusterHeatmap')
               ),
      
      menuItem('Differential Expression', tabName = 'DE_menus', icon = icon('bar-chart-o'), 
                 menuSubItem('DE between clusters', tabName = 'DE_between_clusters'),
                 menuSubItem('DE between clusters K-means', tabName = 'DE_between_clusters_Kmeans')
                ),
      menuItem('Cluster', tabName = 'Cluster', icon = icon('line-chart'),
               menuSubItem('K-means', tabName = 'KMEANS'),
               menuSubItem('Dynamic Tree Cut', tabName = 'Dynamic_Tree_Cut'),
               menuItem('Graph Based Cluster', tabName = 'Graph',
                        menuSubItem('KNN', tabName = 'KNN'),
                        menuSubItem('SNN', tabName = 'SNN'),
                        menuSubItem('Louvain Cluster', tabName = 'GraphBasedCluster'))
      ),
      menuItem("Marker genes", tabName = "MarkerGenes", icon = icon("hornbill")),
      menuItem("Enriched pathway", tabName = "Enriched_pathway", icon = icon("hubspot"))
      
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML('
        /* logo */
                              .skin-blue .main-header .logo {
                              background-color: #800080;
                              }
                              
                              /* logo when hovered */
                              .skin-blue .main-header .logo:hover {
                              background-color: #800080;
                              }
                              
                              /* navbar (rest of the header) */
                              .skin-blue .main-header .navbar {
                              background-color: #800080;
                              }        
                              
                              
                              
                              /* other links in the sidebarmenu when hovered */
                              .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                              background-color: #800080;
                              }
                              /* toggle button when hovered  */                    
                              .skin-blue .main-header .navbar .sidebar-toggle:hover{
                              background-color: #800080;
                              }
                              .box.box-solid.box-primary>.box-header {
                               color:#fff;
                               background:#800080
                              }
                              
                              .box.box-solid.box-primary{
                              
                              border-bottom-color:#080080;
                              border-left-color:#080080;
                              border-right-color:#080080;
                              border-top-color:#080080
                              }

                              '))),
   
    tabItems(
      tabItem(tabName = 'home',
              
              fluidRow(
                column(12,tags$h1('Load data')),
                column(12,
                  fileInput('fileInput', 'HDF5 file location', multiple = FALSE, accept = NULL,
                              width = NULL, buttonLabel = "Browse...",
                              placeholder = "No file selected")),
                
                br(),br(),br(),br(),br(),br(),br(),
                column(4, align="center",
                box(
                  title = tags$p(style = "font-size: 200%;font-weight: 900;", dim(cdScFiltAnnot)[2]), width=20,  status = "primary", "CELLS", solidHeader = TRUE
                )),
                column(4,align="center",
                box(
                  title = tags$p(style = "font-size: 200%;font-weight: 900;", nlevels(as.factor(cdScFiltAnnot$Sample))), width=20,  status = "primary","SAMPLES", solidHeader = TRUE
                  
                )),
                column(4, align="center",
                box(
                  title = tags$p(style = "font-size: 200%;font-weight: 900;", nlevels(as.factor(cdScFiltAnnot$Clusters))), width=20,status = "primary", "CLUSTERS", solidHeader = TRUE
                  
                ))
                    
                    #downloadButton("exportTsne", label = "Download t-SNE"),
                    #downloadButton("exportUmap", label = "Download UMAP")
              
            )
              
      ),
      tabItem(tabName = 'Overview',
              fluidRow(
                box(
                  title = p(tags$span("Overview", style="padding-right:8px;"), 
                           actionButton("projectionInfo", "help", 
                                        class = "btn-xs", title = "Additional information for this panel")
                  ), status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  column(width=6, selectInput("projection", "Overview",
                                              choices = c('tSNE','UMAP'),
                                              multiple = FALSE,
                                              selectize = FALSE)),
                  column(width = 6, selectInput("colorCellsBy", "Color cells by",
                                                choices = c('Sample','Cluster','nUMI','nGene','percent_mt','CellType'),
                                                multiple = FALSE,
                                                selectize = FALSE)),
                  column(width = 6, sliderInput("plotOverviewDotSize", "Dot size:", 0, 10, 3.5, 0.5)),
                  column(width = 6, sliderInput("plotOverviewDotOpacity", "Dot opacity:", 0, 1, 1, 0.1)),
                  column(width=12,plotlyOutput("tsnePlotCluster", width = "100%") %>% withSpinner(type = getOption("spinner.type", default = 8))),
                  column(width = 4, pickerInput(
                    inputId = "checkboxProjectionSample", 
                    label = "Select/deselect samples", 
                    choices = unique(levels(as.factor(cdScFiltAnnot$Sample))), 
                    selected = unique(levels(as.factor(cdScFiltAnnot$Sample))), 
                    options = list(
                      `actions-box` = TRUE, 
                      size = 10,
                      `selected-text-format` = "count > 20"
                    ), 
                    multiple = TRUE
                  )),
                  column(width = 4, pickerInput(
                    inputId = "checkboxProjectionCluster", 
                    label = "Select/deselect clusters", 
                    choices = unique(levels(as.factor(cdScFiltAnnot$Clusters))), 
                    selected = unique(levels(as.factor(cdScFiltAnnot$Clusters))),
                    options = list(
                      `actions-box` = TRUE, 
                      size = 10,
                      `selected-text-format` = "count > 20"
                    ), 
                    multiple = TRUE
                  )),
                  column(width = 4, pickerInput(
                    inputId = "checkboxProjectionCellType", 
                    label = "Select/deselect celltype", 
                    choices = unique(levels(as.factor(cdScFiltAnnot$cellType))), 
                    selected = unique(levels(as.factor(cdScFiltAnnot$cellType))),
                    options = list(
                      `actions-box` = TRUE, 
                      size = 10,
                      `selected-text-format` = "count > 20"
                    ), 
                    multiple = TRUE
                  ))
                )
              )
      ),
     
  
      tabItem(tabName = 'Gene_expression',
              fluidRow(
                box(
                  title = "Gene expression", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  column(width=6,selectInput("geneName", "Gene Name", selected = GeneNameSorted[1], 
                              choices = GeneNameSorted,
                              multiple = FALSE,
                              selectize = FALSE)),
                  column(width=6, selectInput("geneExprProjection", "Projection",
                              choices = c('tSNE','UMAP'),
                              multiple = FALSE,
                              selectize = FALSE)),
                  column(width=6, sliderInput("geneExpressionplotOverviewDotSize", "Dot size:", 0, 10, 0.5, 0.5)),
                  column(width=6,sliderInput("geneExpressionplotOverviewDotOpacity", "Dot opacity:", 0, 1, 1, 0.1)),
                  column(width=6, colourpicker::colourInput("colmaxgeneExp", "Select colour for maximum value", "firebrick1")),
                  column(width=6, colourpicker::colourInput("colmingeneExp", "Select colour for minimum", "gray88")),
                  column(width=12,plotOutput("PlotGeneExpr", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))),
                  column(width = 4, pickerInput(
                    inputId = "checkboxGeneExpressionSample", 
                    label = "Select/deselect samples", 
                    choices = unique(levels(as.factor(cdScFiltAnnot$Sample))), 
                    selected = unique(levels(as.factor(cdScFiltAnnot$Sample))), 
                    options = list(
                      `actions-box` = TRUE, 
                      size = 10,
                      `selected-text-format` = "count > 20"
                    ), 
                    multiple = TRUE
                  )),
                  column(width = 4, pickerInput(
                    inputId = "checkboxGeneExpressionCluster", 
                    label = "Select/deselect clusters", 
                    choices = unique(levels(as.factor(cdScFiltAnnot$Clusters))), 
                    selected = unique(levels(as.factor(cdScFiltAnnot$Clusters))),
                    options = list(
                      `actions-box` = TRUE, 
                      size = 10,
                      `selected-text-format` = "count > 20"
                    ), 
                    multiple = TRUE
                  )),
                  column(width = 4, pickerInput(
                    inputId = "checkboxGeneExpressionCellType", 
                    label = "Select/deselect celltype", 
                    choices = unique(levels(as.factor(cdScFiltAnnot$cellType))), 
                    selected = unique(levels(as.factor(cdScFiltAnnot$cellType))),
                    options = list(
                      `actions-box` = TRUE, 
                      size = 10,
                      `selected-text-format` = "count > 20"
                    ), 
                    multiple = TRUE
                  ))
                ),
                
              )
      ),
      
      
      tabItem(tabName = 'clusterHeatmap',
              fluidRow(
                box(title='Input for cluster heatmap', status = "primary", solidHeader = TRUE,
                    collapsible = TRUE, width = 4,
                    selectizeInput('ChooseClustersHeatmap', 'Select your clusters', 
                                   choices = as.factor(cdScFiltAnnot$Clusters[order(cdScFiltAnnot$Clusters)]), 
                                   selected = '', 
                                   multiple = TRUE),
                    textAreaInput('ClusterGeneList', 'Paste your Gene list', value = '', width = "180px", height = "250px"),
                    tags$form(
                      actionButton("buttonForClusterHeatmap", "Generate Heatmap", styleclass = "primary")
                    )
                    
                    #downloadButton("exportTsne", label = "Download t-SNE"),
                    #downloadButton("exportUmap", label = "Download UMAP")
                ),
                
                box(
                  title = "Heatmap", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 8,
                  plotOutput("plotClusterHeatmap", height = "600px"),
                  column(width=4, colourpicker::colourInput("colmaxClustHeatmap", "Maximum value", "red")),
                  column(width=4, colourpicker::colourInput("colmidClustHeatmap", "Middle value", "white")),
                  column(width=4, colourpicker::colourInput("colminClustHeatmap", "Minimum value", "blue"))
                )
              )
      ),
      tabItem(tabName = 'sampleHeatmap',
             fluidRow(
               box(title='Input for sample heatmap', status = "primary", solidHeader = TRUE,
                   collapsible = TRUE, width = 4,
                   selectizeInput('ChooseSampleHeatmap', 'Select your samples', 
                                  choices = as.factor(cdScFiltAnnot$Sample[order(cdScFiltAnnot$Sample)]), 
                                  selected = '', 
                                  multiple = TRUE),
                   textAreaInput('SampleGeneList', 'Paste your Gene list', value = '', width = "180px", height = "250px"),
                   tags$form(
                     actionButton("buttonForSampleHeatmap", "Generate Heatmap", styleclass = "primary")
                   )
                   
                   #downloadButton("exportTsne", label = "Download t-SNE"),
                   #downloadButton("exportUmap", label = "Download UMAP")
               ),
               
               box(
                 title = "Heatmap", status = "primary", solidHeader = TRUE,
                 collapsible = TRUE, width = 8,
                 plotOutput("plotSampleHeatmap", height = "600px"),
                 column(width=4, colourpicker::colourInput("colmaxSampleHeatmap", "Maximum value", "red")),
                 column(width=4, colourpicker::colourInput("colmidSampleHeatmap", "Middle value", "white")),
                 column(width=4, colourpicker::colourInput("colminSampleHeatmap", "Minimum value", "blue"))
               )
             )
      ),
      tabItem(tabName = 'sampleClusterHeatmap',
              fluidRow(
                box(title='Input for sample & cluster heatmap', status = "primary", solidHeader = TRUE,
                    collapsible = TRUE, width = 4,
                    selectizeInput('ChooseSampleClusterHeatmapSample', 'Select your samples', 
                                   choices = as.factor(cdScFiltAnnot$Sample[order(cdScFiltAnnot$Sample)]), 
                                   selected = '', 
                                   multiple = TRUE),
                    selectizeInput('ChooseSampleClusterHeatmapCluster', 'Select your clusters', 
                                   choices = as.factor(cdScFiltAnnot$Clusters[order(cdScFiltAnnot$Clusters)]), 
                                   selected = '', 
                                   multiple = TRUE),
                    textAreaInput('SampleClusterGeneList', 'Paste your Gene list', value = '', width = "180px", height = "250px"),
                    tags$form(
                      actionButton("buttonForSampleClusterHeatmap", "Generate Heatmap", styleclass = "primary")
                    )
                    
                    #downloadButton("exportTsne", label = "Download t-SNE"),
                    #downloadButton("exportUmap", label = "Download UMAP")
                ),
                
                box(
                  title = "Heatmap", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 8,
                  plotOutput("plotSampleClusterHeatmap", height = "600px"),
                  column(width=4, colourpicker::colourInput("colmaxSampleClusterHeatmap", "Maximum value", "red")),
                  column(width=4, colourpicker::colourInput("colmidSampleClusterHeatmap", "Middle value", "white")),
                  column(width=4, colourpicker::colourInput("colminSampleClusterHeatmap", "Minimum value", "blue"))
                )
              )
      ),
      
      tabItem(tabName = 'sampleBubblePlot',
              fluidRow(
                box(title='Input for sample bubbleplot', status = "primary", solidHeader = TRUE,
                    collapsible = TRUE, width = 4,
                    selectizeInput('ChooseSampleBubblePlot', 'Select your samples', 
                                   choices = as.factor(cdScFiltAnnot$Sample[order(cdScFiltAnnot$Sample)]), 
                                   selected = '', 
                                   multiple = TRUE),
                    textAreaInput('sampleGeneListSampleBubblePlot', 'Paste your Gene list', value = '', width = "180px", height = "250px"),
                    tags$form(
                      actionButton("buttonForSampleBubbleplot", "Generate Bubbleplot", styleclass = "primary")
                    )
                    
                    #downloadButton("exportTsne", label = "Download t-SNE"),
                    #downloadButton("exportUmap", label = "Download UMAP")
                ),
                
                box(
                  title = "Bubble plot", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 8,
                  plotOutput("plotSampleBubbleplot", height = "600px")%>% withSpinner(type = getOption("spinner.type", default = 8))
                )
              )
      ),
      tabItem(tabName = 'DE_menus',
              h2('Selected Sub-Item Four')
      ),
      tabItem(tabName = 'DE_between_clusters',
              fluidRow(
                box(title='Input for DE between clusters', status = "primary", solidHeader = TRUE,
                    collapsible = TRUE, width = 12,
                      column(width = 6,
                             selectizeInput("clust1", "Select Cluster/Clusters for Group1",
                                            choices = as.factor(cdScFiltAnnot$Clusters[order(cdScFiltAnnot$Clusters)]),
                                            multiple = TRUE)),
                      column(width = 6,
                             selectizeInput("clust2", "Select Cluster/Clusters for Group2",
                                            choices = as.factor(my.clusters[order(my.clusters)]),
                                            multiple = TRUE)),
                    
                      column(width = 4,
                        tags$form(
                        actionButton("buttonForTwoClustDE", "Run DE", styleclass = "primary")
                      )),
                    column(width = 12,
                           plotlyOutput("AllClustTwoClustComp") %>% withSpinner(type = getOption("spinner.type", default = 8))),
                    column(width = 12,
                           plotOutput("plotSelectClust") %>% withSpinner(type = getOption("spinner.type", default = 8)) )
                    
                    #downloadButton("exportTsne", label = "Download t-SNE"),
                    #downloadButton("exportUmap", label = "Download UMAP")
                ),
                
                box(
                  title = "DE results", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  DT::dataTableOutput("mytableTwoClust") %>% withSpinner(type = getOption("spinner.type", default = 8))
                )
              )
      ),
      tabItem(tabName = 'DE_between_clusters_Kmeans',
              fluidRow(
                box(title='Input for DE between clusters', status = "primary", solidHeader = TRUE,
                    collapsible = TRUE, width = 12,
                    column(width= 6,
                           selectizeInput("chooseClustType", "Select type of Cluster",
                                          choices = c('Kmeans','Graph Based Cluster','Dynamic Tree Cut','Default Cluster'),
                                          selected = (names(list_cluster)[1]),
                                          multiple = FALSE)),
                    column(width= 6,
                          selectizeInput("chooseClust", "Select Cluster resolution",
                                         choices = (names(list_cluster)),
                                         selected = (names(list_cluster)[1]),
                                         multiple = FALSE)),
                    column(width = 6,
                           selectizeInput("clust1forK", "Select Cluster",
                                          choices = as.factor(as.numeric(cdScFiltAnnot$kmeansCluster[order(cdScFiltAnnot$kmeansCluster)])),
                                          multiple = FALSE)),
                    column(width = 6,
                           selectizeInput("clust2forK", "Select Cluster",
                                          choices = c(1,2,3,4,5,6,7,8,9,10),
                                          multiple = FALSE)),
                    column(width = 4,
                           tags$form(
                             actionButton("buttonForTwoClustDEKmeans", "Run DE", styleclass = "primary")
                           )),
                    #column(width = 12,
                     #      plotlyOutput("AllClustTwoClustComp") %>% withSpinner(type = getOption("spinner.type", default = 8))),
                    #column(width = 12,
                     #      plotOutput("plotSelectClust") %>% withSpinner(type = getOption("spinner.type", default = 8)) )
                    
                    #downloadButton("exportTsne", label = "Download t-SNE"),
                    #downloadButton("exportUmap", label = "Download UMAP")
                ),
                
                box(
                  title = "DE results", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  DT::dataTableOutput("mytableTwoClustKmeans") %>% withSpinner(type = getOption("spinner.type", default = 8))
                )
              )
      ),
     tabItem(tabName = 'KMEANS',
          fluidRow(
             box( 
               title = "Kmeans With tSNE", status = "primary", solidHeader = TRUE,
               collapsible = TRUE, width = 12,
               
               column(width = 6,
                     numericInput('choice_K', 'select K', 12, min = 2, max = 50, step = 1)
               ),
               #column(width = 4,
                #    tags$form(
                 #     actionButton("buttonForKmeans", "Run Kmeans", styleclass = "primary")
                  #  )),
               column(width = 12,
               plotlyOutput("Cluster_KMEANSwithtSNE", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8)))            
             ),
            
               #box(
               #title = "Kmeans", status = "primary", solidHeader = TRUE,
               #collapsible = TRUE, width = 12,
               #plotlyOutput("Cluster_KMEANS", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))            
             #)
            )
             
     ),
     tabItem(tabName = 'KNN',
             box(
               title = "K Nearest Neighbor", status = "primary", solidHeader = TRUE,
               collapsible = TRUE, width = 12,
               column(width = 6,
                      numericInput('choice_K', 'select K', 12, min = 2, max = 50, step = 1)
               ),
               column(width = 12,
               plotOutput("KNN", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8)) 
               )
             )
     ),
     tabItem(tabName = 'SNN',
             box(
               title = "Shared Nearest Neighbor", status = "primary", solidHeader = TRUE,
               collapsible = TRUE, width = 12,
               column(width = 6,
                      numericInput('choice_K', 'select K', 12, min = 2, max = 50, step = 1)
               ),
               column(width = 12,
               plotOutput("SNN", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))     
               ))
     ),
     tabItem(tabName = 'GraphBasedCluster',
             fluidRow(
               box(
                 title = "Graph Based Cluster", status = "primary", solidHeader = TRUE,
                 collapsible = TRUE, width = 12,
                 column(width = 6,
                        numericInput('choice_K_for_Louvain', 'select K', 12, min = 2, max = 50, step = 1)
                 ),
                 column(width = 12,
                 plotlyOutput("GraphBasedLouvainCluster", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))            
                 )
                 ),
            
               #box(
                # title = "Graph View of Graph Based Cluster", status = "primary", solidHeader = TRUE,
                 #collapsible = TRUE, width = 12,
                 #plotOutput("GaphBased_Graph", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))            
               #)
             )
             
     ),
     tabItem(tabName = 'Dynamic_Tree_Cut',
             box(
               title = "Shared Nearest Neighbor", status = "primary", solidHeader = TRUE,
               collapsible = TRUE, width = 12,
               plotlyOutput("Dynamic_Tree_Cut_Plot", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))             )
     ),
     tabItem(tabName = 'MarkerGenes',
             fluidRow(
               box(
                 title = "Marker genes for clusters", status = "primary", solidHeader = TRUE,
                 collapsible = TRUE, width = 12,
                 selectInput('markerChooseCluster', 'Default Cluster', 
                             choices = levels(as.factor(my.clusters[order(my.clusters)])), 
                             multiple = FALSE,
                             selectize = FALSE),
                 DT::dataTableOutput("markerTableCluster") %>% withSpinner(type = getOption("spinner.type", default = 8))
               ),
               box(
                 title = "Marker genes for clusters", status = "primary", solidHeader = TRUE,
                 collapsible = TRUE, width = 12,
                 column(width= 6,
                        selectizeInput("chooseClustTypeforMarker", "Select type of Cluster",
                                       choices = c('Kmeans','Graph Based Cluster','Dynamic Tree Cut','Default Cluster'),
                                       selected = (names(list_cluster)[1]),
                                       multiple = FALSE)),
                 column(width= 6,
                        selectizeInput("ClusterRes", "Select Cluster resolution",
                                       choices = (names(list_cluster)),
                                       selected = (names(list_cluster)[1]),
                                       multiple = FALSE)),
                 column(width= 6,
                 selectizeInput("markerChoose_Kmeans_Cluster", "choose Cluster", 
                                choices = levels(as.factor(as.numeric(cdScFiltAnnot$kmeansCluster[order(cdScFiltAnnot$kmeansCluster)]))), 
                                selected = levels(as.factor(as.numeric(cdScFiltAnnot$kmeansCluster[order(cdScFiltAnnot$kmeansCluster)])[1])),
                                multiple = FALSE
                                )),
                 
                 DT::dataTableOutput("markerTableKmeansCluster") %>% withSpinner(type = getOption("spinner.type", default = 8))
               )
             )
     ),
     tabItem(tabName = 'Enriched_pathway',
             fluidRow(
              
               box(
                 title = "Pathway enrichment for clusters", status = "primary", solidHeader = TRUE,
                 collapsible = TRUE, width = 12,
                 column(width=6,selectInput('ENChooseCluster', 'Cluster', 
                                            choices = levels(as.factor(my.clusters[order(my.clusters)])), 
                                            multiple = FALSE,
                                            selectize = FALSE)),
                 # column(width=6, selectInput('ENChooseCluster', 'Cluster', 
                 #                             choices = levels(as.factor(my.clusters[order(my.clusters)])), 
                 #                             selected = '', multiple = FALSE,
                 #                             selectize = FALSE)),
                 column(width = 6, selectInput('ENChooseClusterGO', 'Cluster', 
                                               choices = levels(as.factor(my.clusters[order(my.clusters)])), 
                                               selected = '', multiple = FALSE,
                                               selectize = FALSE)),
                 column(width=12,DT::dataTableOutput("ENTableCluster") %>% withSpinner(type = getOption("spinner.type", default = 8)))
               ),
               box(
                 title = "Pathway enrichment for clusters", status = "primary", solidHeader = TRUE,
                 collapsible = TRUE, width = 12,
                 column(width=6,selectInput('ENChooseKmeansCluster', 'Cluster', 
                                            choices = levels(as.factor(as.numeric(cdScFiltAnnot$kmeansCluster[order(cdScFiltAnnot$kmeansCluster)]))), 
                                            multiple = FALSE,
                                            selectize = FALSE,
                                            selected = levels(as.factor(as.numeric(cdScFiltAnnot$kmeansCluster[order(cdScFiltAnnot$kmeansCluster)])))[3])),
                 # column(width=6, selectInput('ENChooseKmeansCluster', 'Cluster', 
                 #                             choices = levels(as.factor(my.clusters[order(my.clusters)])), 
                 #                             selected = '', multiple = FALSE,
                 #                             selectize = FALSE)),
                 column(width = 6, selectInput('ENChooseKmeansClusterGO', 'Cluster', 
                                               choices = levels(as.factor(as.numeric(cdScFiltAnnot$kmeansCluster[order(cdScFiltAnnot$kmeansCluster)]))), 
                                               selected = levels(as.factor(as.numeric(cdScFiltAnnot$kmeansCluster[order(cdScFiltAnnot$kmeansCluster)]))),
                                               multiple = FALSE,
                                               selectize = FALSE)),
                 column(width=12,DT::dataTableOutput("ENTableKmeansCluster") %>% withSpinner(type = getOption("spinner.type", default = 8)))
               ),
               
             )
     )
    
     
  )
 )
)


server <- function(input, output, session) { 
  
  #################################
  #
  # For overview panel
  #
  #################################
  
  
  
  
  output$violinPlot <- renderPlot({
    
    geneListToTest <- input$geneName
    #cdScFiltAnnotTmp <- cdScFiltAnnot
    dfViolin <- data.frame(Sample=colData(cdScFiltAnnot)$Sample, logCounts=logcounts(cdScFiltAnnot)[geneListToTest,], title=geneListToTest)
    p <- ggplot(dfViolin, aes(factor(Sample), logCounts)) +
      geom_violin(aes(fill=factor(Sample)), scale="width", trim=TRUE) +
      xlab('Sample') +
      ylab('Expression (logcounts)') +
      scale_fill_manual(values = c_sample_col) +
      #      scale_fill_manual(values = c("#E69F00", "#009E73", "#CC79A7")) + 
      theme_classic(base_size=14) +
      guides(fill=guide_legend(title="Samples")) + 
      facet_grid(. ~ title)
    p
    
    
  })
  
  
  output$violinPlotClusterOrig <- renderPlot({
    
    geneListToTest <- input$geneName
  
    
    dfViolin <- data.frame(Cluster=colData(cdScFiltAnnot)$Clusters, logCounts=logcounts(cdScFiltAnnot)[geneListToTest,], title=geneListToTest)
    p <- ggplot(dfViolin, aes(factor(Cluster), logCounts)) +
      geom_violin(aes(fill=factor(Cluster)), scale="width", trim=TRUE) +
      xlab('Cluster') +
      ylab('Expression (logcounts)') +
      scale_fill_manual(values = c_clust_col) +
      #      scale_fill_manual(values = c("#F8766D", "#DB8E00", "#AEA200", "#64B200", "#00BD5C", "#00C1A7", "#00BADE", "#00A6FF"))+
      theme_classic(base_size=14) +
      guides(fill=guide_legend(title="Clusters")) + 
      facet_grid(. ~ title)
    p
    
  })
  
  
  plot_tsnePlotCluster <- function(){ 	
    
    
    # df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    # df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    # df[,'Clusters'] <- as.factor(cdScFiltAnnot$Clusters)
    # as_tibble(df) %>%
    #    filter(Clusters %in% input$checkboxProjectionCluster) %>%
    #    filter(Sample %in% input$checkboxProjectionSample) %>%
    # ggplot(aes(x=V1, y=V2, Clusters=Clusters)) +
    #   geom_point(size=input$plotOverviewDotSize,alpha=input$plotOverviewDotOpacity, aes(colour=Clusters)) +
    #   #scale_colour_manual(values=cbPalette) +
    #   guides(colour = guide_legend(override.aes = list(size=4))) +
    #   xlab("") + ylab("") +
    #   ggtitle("t-SNE 2D coloured by cluster") +
    #   theme_classic(base_size=10)  +
    #   scale_color_manual(values = c_clust_col)
    
    
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    df[,'Cell']=as.factor(colData(cdScFiltAnnot)$Barcode)
    df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    df[,'Clusters'] <- as.factor(my.clusters)
    df[,'cellType'] <- as.factor(colData(cdScFiltAnnot)$cellType)
    df %>%
      filter(Clusters %in% input$checkboxProjectionCluster) %>%
      filter(Sample %in% input$checkboxProjectionSample) %>%
      filter(cellType %in% input$checkboxProjectionCellType) %>%
      plot_ly(color = ~Clusters, colors = c_clust_col[c(1:9)], type="scatter", mode="markers", hoverinfo = 'text',
                   marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
                   text = ~paste('</br> Cell: ', Cell,
                                 '</br> Clusters: ', Clusters,
                                 '</br> Samples: ', Sample,
                                 '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>%
      layout(title = 'tSNE with Clusters',showlegend = TRUE, legend = list(font = list(size = 10), itemsizing='constant'))
    
  }
  
  plot_tsnePlotCellType <- function(){ 	
    
    
    # df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    # df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    # df[,'Clusters'] <- as.factor(cdScFiltAnnot$Clusters)
    # as_tibble(df) %>%
    #    filter(Clusters %in% input$checkboxProjectionCluster) %>%
    #    filter(Sample %in% input$checkboxProjectionSample) %>%
    # ggplot(aes(x=V1, y=V2, Clusters=Clusters)) +
    #   geom_point(size=input$plotOverviewDotSize,alpha=input$plotOverviewDotOpacity, aes(colour=Clusters)) +
    #   #scale_colour_manual(values=cbPalette) +
    #   guides(colour = guide_legend(override.aes = list(size=4))) +
    #   xlab("") + ylab("") +
    #   ggtitle("t-SNE 2D coloured by cluster") +
    #   theme_classic(base_size=10)  +
    #   scale_color_manual(values = c_clust_col)
    
    
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    df[,'Cell']=as.factor(colData(cdScFiltAnnot)$Barcode)
    df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    df[,'Clusters'] <- as.factor(my.clusters)
    df[,'cellType'] <- as.factor(colData(cdScFiltAnnot)$cellType)
    df %>%
      filter(Clusters %in% input$checkboxProjectionCluster) %>%
      filter(Sample %in% input$checkboxProjectionSample) %>%
      filter(cellType %in% input$checkboxProjectionCellType) %>%
      plot_ly(color = ~cellType, colors = c30, type="scatter", mode="markers", hoverinfo = 'text',
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Sample,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>%
      layout(title = 'tSNE with Clusters',showlegend = TRUE, legend = list(font = list(size = 10), itemsizing='constant'))
    
  }
  
  plot_tsnePlotSample <- function(){ 	
    
    
    # df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    # df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    # df[,'Clusters'] <- as.factor(cdScFiltAnnot$Clusters)
    # as_tibble(df) %>%
    #   filter(Clusters %in% input$checkboxProjectionCluster) %>%
    #   filter(Sample %in% input$checkboxProjectionSample) %>%
    #   ggplot( aes(x=V1, y=V2, Sample=Sample)) +
    #   geom_point(size=input$plotOverviewDotSize,alpha=input$plotOverviewDotOpacity, aes(colour=Sample)) +
    #   #scale_colour_manual(values=cbPalette) +
    #   guides(colour = guide_legend(override.aes = list(size=4))) +
    #   xlab("") + ylab("") +
    #   ggtitle("t-SNE 2D coloured by sample") +
    #   theme_classic(base_size=10)  +
    #   scale_color_manual(values = c_sample_col)
    
    
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    df[,'Cell']=as.factor(colData(cdScFiltAnnot)$Barcode)
    df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    df[,'Clusters'] <- as.factor(my.clusters)
    df[,'cellType'] <- as.factor(colData(cdScFiltAnnot)$cellType)
    df %>%
      filter(Clusters %in% input$checkboxProjectionCluster) %>%
      filter(Sample %in% input$checkboxProjectionSample) %>%
      filter(cellType %in% input$checkboxProjectionCellType) %>%
      plot_ly(color = ~Sample, colors = c_sample_col[c(1:4)], type="scatter", mode="markers", hoverinfo = 'text',
                   marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
                   text = ~paste('</br> Cell: ', Cell,
                                 '</br> Clusters: ', Clusters,
                                 '</br> Samples: ', Sample,
                                 '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>%
      layout(title = 'tSNE with Samples',showlegend = TRUE, legend = list(font = list(size = 10), itemsizing='constant'))

  }
  plot_tsnePlotnUMI <- function(){ 	
    
    as.tibble(reducedDim(cdScFiltAnnot,'tSNE')) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Samples = colData(cdScFiltAnnot)$Sample) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Clusters = as.factor(my.clusters)) %>%
      mutate(log10_nUMI = log10(cdScFiltAnnot$total)) %>%
      mutate(cellType = colData(cdScFiltAnnot)$cellType) %>%
      plot_ly(color = ~log10_nUMI, type="scatter", mode="markers", hoverinfo = 'text', 
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Samples,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>% 
      layout(title = 'tSNE with total UMI in Log10',showlegend = FALSE)
    
  }
  plot_tsnePlotnGene <- function(){ 	
    
    as.tibble(reducedDim(cdScFiltAnnot,'tSNE')) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Samples = colData(cdScFiltAnnot)$Sample) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Clusters = as.factor(my.clusters)) %>%
      mutate(log10_nGenes = log10(cdScFiltAnnot$detected)) %>%
      mutate(cellType = colData(cdScFiltAnnot)$cellType) %>%
      plot_ly(color = ~log10_nGenes, type="scatter", mode="markers", hoverinfo = 'text' , 
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Samples,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>% 
      layout(title = 'tSNE with total gene in log10',showlegend = FALSE)
    
  }
  
  plot_tsnePlotPct_MT <- function(){ 	
    
    as.tibble(reducedDim(cdScFiltAnnot,'tSNE')) %>%
      mutate(Samples = colData(cdScFiltAnnot)$Sample) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Clusters = as.factor(my.clusters)) %>%
      mutate(pct_mt = cdScFiltAnnot$subsets_Mt_percent) %>%
      mutate(cellType = colData(cdScFiltAnnot)$cellType) %>%
      plot_ly(color = ~pct_mt, type="scatter", mode="markers", hoverinfo = 'text' , 
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Samples,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>% 
      layout(title = 'tSNE with pct_mt',showlegend = FALSE)
    
  }
  
  
  
  plot_UMAPPlotCluster <- function(){    
    
    # df <- as.data.frame(reducedDim(cdScFiltAnnot,'UMAP'))
    # df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    # df[,'Clusters'] <- as.factor(cdScFiltAnnot$Clusters)
    # as_tibble(df) %>%
    #   filter(Clusters %in% input$checkboxProjectionCluster) %>%
    #   filter(Sample %in% input$checkboxProjectionSample) %>%
    #   ggplot(aes(x=V1, y=V2, Clusters=Clusters)) +
    #   geom_point(size=input$plotOverviewDotSize,alpha=input$plotOverviewDotOpacity, aes(colour=Clusters)) +
    #   #scale_colour_manual(values=cbPalette) +
    #   guides(colour = guide_legend(override.aes = list(size=4))) +
    #   xlab("") + ylab("") +
    #   ggtitle("UMAP 2D coloured by cluster") +
    #   theme_classic(base_size=10)  +
    #   scale_color_manual(values = c_clust_col)
    
    as.tibble(reducedDim(cdScFiltAnnot,'UMAP')) %>%
      mutate(Sample = colData(cdScFiltAnnot)$Sample) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Clusters = as.factor(my.clusters)) %>%
      mutate(cellType = colData(cdScFiltAnnot)$cellType) %>%
      filter(Clusters %in% input$checkboxProjectionCluster) %>%
      filter(Sample %in% input$checkboxProjectionSample) %>%
      filter(cellType %in% input$checkboxProjectionCellType) %>%
      plot_ly(color = ~cellType, colors = c_clust_col[c(1:8)], type="scatter", mode="markers", hoverinfo = 'text',
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Sample,
                            '</br> CellType: ', cellType)) %>%
      add_trace(x=~V1,y=~V2) %>%
      layout(title = 'UMAP with Clusters',showlegend = TRUE, legend = list(title='Cluster',font = list(size = 10), itemsizing='constant'))
    
  }
  
  plot_UMAPPlotCellType <- function(){    
    
    # df <- as.data.frame(reducedDim(cdScFiltAnnot,'UMAP'))
    # df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    # df[,'Clusters'] <- as.factor(cdScFiltAnnot$Clusters)
    # as_tibble(df) %>%
    #   filter(Clusters %in% input$checkboxProjectionCluster) %>%
    #   filter(Sample %in% input$checkboxProjectionSample) %>%
    #   ggplot(aes(x=V1, y=V2, Clusters=Clusters)) +
    #   geom_point(size=input$plotOverviewDotSize,alpha=input$plotOverviewDotOpacity, aes(colour=Clusters)) +
    #   #scale_colour_manual(values=cbPalette) +
    #   guides(colour = guide_legend(override.aes = list(size=4))) +
    #   xlab("") + ylab("") +
    #   ggtitle("UMAP 2D coloured by cluster") +
    #   theme_classic(base_size=10)  +
    #   scale_color_manual(values = c_clust_col)
    
    as.tibble(reducedDim(cdScFiltAnnot,'UMAP')) %>%
      mutate(Sample = colData(cdScFiltAnnot)$Sample) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Clusters = as.factor(my.clusters)) %>%
      mutate(cellType = colData(cdScFiltAnnot)$cellType) %>%
      filter(Clusters %in% input$checkboxProjectionCluster) %>%
      filter(Sample %in% input$checkboxProjectionSample) %>%
      filter(cellType %in% input$checkboxProjectionCellType) %>%
      plot_ly(color = ~Clusters, colors = c_clust_col[c(1:8)], type="scatter", mode="markers", hoverinfo = 'text',
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Sample,
                            '</br> CellType: ', cellType)) %>%
      add_trace(x=~V1,y=~V2) %>%
      layout(title = 'UMAP with Clusters',showlegend = TRUE, legend = list(title='Cluster',font = list(size = 10), itemsizing='constant'))
    
  }
  
  plot_UMAPPlotSample <- function(){    
    
    
    # df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    # df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    # df[,'Clusters'] <- as.factor(cdScFiltAnnot$Clusters)
    # as_tibble(df) %>%
    #   filter(Clusters %in% input$checkboxProjectionCluster) %>%
    #   filter(Sample %in% input$checkboxProjectionSample) %>%
    #   ggplot(aes(x=V1, y=V2, Sample=Sample)) +
    #   geom_point(size=input$plotOverviewDotSize,alpha=input$plotOverviewDotOpacity, aes(colour=Sample)) +
    #   #scale_colour_manual(values=cbPalette) +
    #   guides(colour = guide_legend(override.aes = list(size=4))) +
    #   xlab("") + ylab("") +
    #   ggtitle("t-SNE 2D coloured by cluster") +
    #   theme_classic(base_size=10)  +
    #   scale_color_manual(values = c_clust_col)
    
    as.tibble(reducedDim(cdScFiltAnnot,'UMAP')) %>%
      mutate(Sample = colData(cdScFiltAnnot)$Sample) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Clusters = as.factor(my.clusters)) %>%
      mutate(cellType = colData(cdScFiltAnnot)$cellType) %>%
      filter(Clusters %in% input$checkboxProjectionCluster) %>%
      filter(Sample %in% input$checkboxProjectionSample) %>%
      filter(cellType %in% input$checkboxProjectionCellType) %>%
      plot_ly(color = ~Sample, colors = c_sample_col[c(1:4)], type="scatter", mode="markers", hoverinfo = 'text',
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Sample,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>%
      layout(title = 'UMAP with Samples',showlegend = TRUE, legend = list(title='Cluster',font = list(size = 10), itemsizing='constant'))
  }
  plot_UMAPPlotnUMI <- function(){    
    
    
      # as.tibble(reducedDim(cdScFiltAnnot,'UMAP')) %>%
      # mutate(Samples = colData(cdScFiltAnnot)$Sample) %>%
      # mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      # mutate(Clusters = as.factor(my.clusters)) %>%
      # mutate(log10_nUMI = log10(cdScFiltAnnot$total)) %>%
      # filter(Clusters %in% input$checkboxProjectionCluster) %>%
      # filter(Sample %in% input$checkboxProjectionSample) %>%
      # ggplot(aes(x=V1, y=V2, colour=log10_nUMI)) +
      # geom_point(size=input$plotOverviewDotSize,alpha=input$plotOverviewDotOpacity) +
      # #scale_colour_manual(values=cbPalette) +
      # guides(colour = guide_legend(override.aes = list(size=4))) +
      # xlab("") + ylab("") +
      # ggtitle("UMAP with total UMI in log10") +
      # theme_classic(base_size=10)
    
    
    as.tibble(reducedDim(cdScFiltAnnot,'UMAP')) %>%
      mutate(Samples = colData(cdScFiltAnnot)$Sample) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Clusters = as.factor(my.clusters)) %>%
      mutate(log10_nUMI = log10(cdScFiltAnnot$total)) %>%
      mutate(cellType = colData(cdScFiltAnnot)$cellType) %>%
      plot_ly(color = ~log10_nUMI, type="scatter", mode="markers", hoverinfo = 'text' ,
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Samples,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>%
      layout(title = 'UMAP with total UMI in log10',showlegend = FALSE)
  }
  
  plot_UMAPPlotnGene <- function(){ 	
    
    as.tibble(reducedDim(cdScFiltAnnot,'UMAP')) %>%
      mutate(Samples = colData(cdScFiltAnnot)$Sample) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Clusters = as.factor(my.clusters)) %>%
      mutate(log10_nGenes = log10(cdScFiltAnnot$detected)) %>%
      mutate(cellType = colData(cdScFiltAnnot)$cellType) %>%
      plot_ly(color = ~log10_nGenes, type="scatter", mode="markers", hoverinfo = 'text', 
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Samples,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>% 
      layout(title = 'UMAP with total gene in log10',showlegend = FALSE)
    
  }
  
  plot_UMAPPlotPercentMt <- function(){ 	
    
    as.tibble(reducedDim(cdScFiltAnnot,'UMAP')) %>%
      mutate(Samples = colData(cdScFiltAnnot)$Sample) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Clusters = as.factor(my.clusters)) %>%
      mutate(pct_mt = cdScFiltAnnot$subsets_Mt_percent) %>%
      mutate(cellType = colData(cdScFiltAnnot)$cellType) %>%
      plot_ly(color = ~pct_mt, type="scatter", mode="markers", hoverinfo = 'text', 
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Samples,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>% 
      layout(title = 'UMAP with pct_mt',showlegend = FALSE)
    
  }
  
  
  output$PlotGeneExpr <- renderPlot({
    if(input$geneExprProjection == 'tSNE')
      plot_tsnePlotGene()
    else if(input$geneExprProjection == 'UMAP')
      plot_UMAPPlotGene()
  })

  
  output$tsnePlotCluster <- renderPlotly({  
    if(input$projection == 'tSNE' & input$colorCellsBy == 'Cluster')
      plot_tsnePlotCluster()
    else if(input$projection == 'tSNE' & input$colorCellsBy == 'Sample')
      plot_tsnePlotSample()
    else if(input$projection == 'tSNE' & input$colorCellsBy == 'CellType')
      plot_tsnePlotCellType()
    else if(input$projection == 'tSNE' & input$colorCellsBy == 'nUMI')
      plot_tsnePlotnUMI()
    else if(input$projection == 'tSNE' & input$colorCellsBy == 'nGene')
      plot_tsnePlotnGene()
    else if(input$projection == 'tSNE' & input$colorCellsBy == 'percent_mt')
      plot_tsnePlotPct_MT()
    else if(input$projection == 'UMAP' & input$colorCellsBy == 'Cluster')
      plot_UMAPPlotCluster()
    else if(input$projection == 'UMAP' & input$colorCellsBy == 'Sample')
      plot_UMAPPlotSample()
    else if(input$projection == 'UMAP' & input$colorCellsBy == 'CellType')
      plot_UMAPPlotCellType()
    else if(input$projection == 'UMAP' & input$colorCellsBy == 'nUMI')
      plot_UMAPPlotnUMI()
    else if(input$projection == 'UMAP' & input$colorCellsBy == 'nGene')
      plot_UMAPPlotnGene()
    else if(input$projection == 'UMAP' & input$colorCellsBy == 'percent_mt')
      plot_UMAPPlotPercentMt()
  })
  
  output$tsnePlotSample <- renderPlotly({   
    if(input$projection == 'tSNE')
      plot_tsnePlotSample()
    else
      plot_UMAPPlotSample()
  })
  
  
  output$checkboxProjectionSample <- renderUI({
    choice <-  unique(levels(as.factor(cdScFiltAnnot$Sample)))
    pickerInput(
      inputId = "checkboxProjectionSample", 
      label = "Select/deselect samples for display", 
      choices = choice, 
      options = list(
        `actions-box` = TRUE, 
        size = 10,
        `selected-text-format` = "count > 3"
      ), 
      multiple = TRUE
    )
    
    
  })
  
  output$checkboxProjectionCluster <- renderUI({
    choice <-  unique(levels(cdScFiltAnnot$Clusters))
    pickerInput(
      inputId = "checkboxProjectionCluster", 
      label = "Select/deselect clusters for display", 
      choices = choice, 
      options = list(
        `actions-box` = TRUE, 
        size = 10,
        `selected-text-format` = "count > 3"
      ), 
      multiple = TRUE
    )
    
  })
  
  ############################################
  #
  # For Gene expression panel
  #
  ############################################
  
  plot_tsnePlotGene <- function(){
    
    as_tibble(reducedDim(cdScFiltAnnot,'tSNE')) %>%
      mutate(GeneExp = as.vector(logcounts(cdScFiltAnnot[input$geneName,]))) %>%
      mutate(Sample = as.factor(cdScFiltAnnot$Sample)) %>%
      mutate(Clusters = as.factor(cdScFiltAnnot$Clusters)) %>%
      mutate(cellType = as.factor(cdScFiltAnnot$cellType)) %>%
      filter(Clusters %in% input$checkboxGeneExpressionCluster) %>%
      filter(Sample %in% input$checkboxGeneExpressionSample) %>%
      filter(cellType %in% input$checkboxGeneExpressionCellType) %>%
      ggplot(aes(x=V1, y=V2, GeneExp = GeneExp)) +
      geom_point(size=input$geneExpressionplotOverviewDotSize,alpha=input$geneExpressionplotOverviewDotOpacity, aes(colour = GeneExp)) +
      #scale_colour_viridis_c()+
      scale_colour_gradientn(colours=c(input$colmingeneExp, input$colmaxgeneExp),
                             guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
      #scale_colour_gradientn(colours=c("grey90", "orangered", "orangered4"))+
      #guides(colour = guide_legend(override.aes = list(size=4))) +
      xlab("") + ylab("") +
      ggtitle(paste0('tSNE Gene Exp:',input$geneName))+
      theme_classic(base_size=14) 
      # theme(strip.background = element_blank(),
      #       strip.text.x     = element_blank(),
      #       axis.text.x      = element_blank(),
      #       axis.text.y      = element_blank(),
      #       axis.ticks       = element_blank(),
      #       axis.line        = element_blank(),
      #       panel.border     = element_blank())
    
  }
  
  
  plot_tsnePlotGeneFacetGene <- function(){
    
    GeneExp <- logcounts(cdScFiltAnnot)[input$geneName,]
    #GeneName = 'SPN'
    
    
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    GeneExp=logcounts(cdScFiltAnnot)[input$geneName,]
    Sample = as.factor(colData(cdScFiltAnnot)$Sample)
    
    as.tibble(df) %>%
      dplyr::mutate("GeneExp" = GeneExp) %>% 
      dplyr::mutate("Sample" = Sample) %>%
      ggplot(aes(x=V1, y=V2, GeneExp = GeneExp)) +
      geom_point(size=0.6,alpha=1, aes(colour = GeneExp)) +
      scale_colour_gradientn(colours=c("gray88", "red"))+
      #scale_colour_gradientn(colours=c("#FDBB84","#FC8D59","#E34A33","#B30000"))
      #scale_colour_viridis_c()+
      #scale_colour_gradientn(colours=c("grey90", "orangered", "orangered4"))+
      #guides(colour = guide_legend(override.aes = list(size=4))) +
      xlab("") + ylab("") +
      ggtitle(paste0('Gene Exp:',input$geneName)) +
      theme_classic(base_size=14) +  
      theme(#strip.background = element_blank(),
        #strip.text.x     = element_blank(),
        axis.text.x      = element_blank(),
        axis.text.y      = element_blank(),
        axis.ticks       = element_blank(),
        axis.line        = element_blank(),
        panel.border     = element_rect(colour = "gray50", fill=NA, size=0.25),
        panel.spacing.x=unit(0, "lines"), 
        panel.spacing.y=unit(0,"lines"),
        #legend.position = ("none")) +
        panel.background =  element_blank(), 
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank()) + 
      facet_wrap(~Sample, ncol=4,nrow=3)
    
  }
  
  
  plot_tsnePlotGeneFacetSample <- function(){
    
    GeneExp <- logcounts(cdScFiltAnnot)[input$geneName,]
    #GeneName = 'SPN'
    
    
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    GeneExp=logcounts(cdScFiltAnnot)[input$geneName,]
    Sample = as.factor(colData(cdScFiltAnnot)$Sample)
    
    as.tibble(df) %>%
      dplyr::mutate("GeneExp" = GeneExp) %>% 
      dplyr::mutate("Sample" = Sample) %>%
      ggplot(aes(x=V1, y=V2, GeneExp = GeneExp)) +
      geom_point(size=0.6,alpha=1, aes(colour = GeneExp)) +
      scale_colour_gradientn(colours=c("gray88", "red"))+
      #scale_colour_gradientn(colours=c("#FDBB84","#FC8D59","#E34A33","#B30000"))
      #scale_colour_viridis_c()+
      #scale_colour_gradientn(colours=c("grey90", "orangered", "orangered4"))+
      #guides(colour = guide_legend(override.aes = list(size=4))) +
      xlab("") + ylab("") +
      ggtitle(paste0('Gene Exp:',input$geneName)) +
      theme_classic(base_size=14) +  
      theme(#strip.background = element_blank(),
        #strip.text.x     = element_blank(),
        axis.text.x      = element_blank(),
        axis.text.y      = element_blank(),
        axis.ticks       = element_blank(),
        axis.line        = element_blank(),
        panel.border     = element_rect(colour = "gray50", fill=NA, size=0.25),
        panel.spacing.x=unit(0, "lines"), 
        panel.spacing.y=unit(0,"lines"),
        #legend.position = ("none")) +
        panel.background =  element_blank(), 
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank()) + 
      facet_wrap(~Sample, ncol=4,nrow=3)
    
  }
  
  plot_tsnePlotGeneFacetCluster <- function(){
    
    #GeneExp <- logcounts(cdScFiltAnnot)[input$geneName,]
    #GeneName = 'SPN'
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    #df[,'GeneExp']=logcounts(cdScFiltAnnot)[input$geneName,]
    #df[,'Clusters'] <- as.factor(my.clusters)
    GeneExp <- logcounts(cdScFiltAnnot)[input$geneName,]
    as.tibble(df) %>%
      dplyr::mutate("GeneExp" = GeneExp) %>% 
      dplyr::mutate("Clusters" = as.factor(my.clusters)) %>%
      ggplot(aes(x=V1, y=V2, GeneExp = GeneExp)) +
      geom_point(size=0.6,alpha=1, aes(colour = GeneExp)) +
      scale_colour_gradientn(colours=c("gray88", "red"))+
      #scale_colour_viridis_c()+
      #scale_colour_gradientn(colours=c("grey90", "orangered", "orangered4"))+
      #guides(colour = guide_legend(override.aes = list(size=4))) +
      xlab("") + ylab("") +
      ggtitle(paste0('Gene Exp:',input$geneName)) +
      theme_classic(base_size=14) +  
      theme(#strip.background = element_blank(),
        #strip.text.x     = element_blank(),
        axis.text.x      = element_blank(),
        axis.text.y      = element_blank(),
        axis.ticks       = element_blank(),
        axis.line        = element_blank(),
        panel.border     = element_rect(colour = "gray50", fill=NA, size=0.25),
        panel.spacing.x=unit(0, "lines"), 
        panel.spacing.y=unit(0,"lines"),
        #legend.position = ("none")) +
        panel.background =  element_blank(), 
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank()) + 
      facet_wrap(~Clusters, ncol=3,nrow=4)
    
  }
  
  
  plot_UMAPPlotGene <- function(){
    
    as_tibble(reducedDim(cdScFiltAnnot,'UMAP')) %>%
      mutate(GeneExp = as.vector(logcounts(cdScFiltAnnot[input$geneName,]))) %>%
      mutate(Sample = as.factor(cdScFiltAnnot$Sample)) %>%
      mutate(Clusters = as.factor(cdScFiltAnnot$Clusters)) %>%
      mutate(cellType = as.factor(cdScFiltAnnot$cellType)) %>%
      filter(Clusters %in% input$checkboxGeneExpressionCluster) %>%
      filter(Sample %in% input$checkboxGeneExpressionSample) %>%
      filter(cellType %in% input$checkboxGeneExpressionCellType) %>%
      ggplot(aes(x=V1, y=V2, GeneExp = GeneExp)) +
      geom_point(size=input$geneExpressionplotOverviewDotSize,alpha=input$geneExpressionplotOverviewDotOpacity, aes(colour = GeneExp)) +
      #scale_colour_viridis_c()+
      scale_colour_gradientn(colours=c(input$colmingeneExp, input$colmaxgeneExp))+
      #scale_colour_gradientn(colours=c("grey90", "orangered", "orangered4"))+
      #guides(colour = guide_legend(override.aes = list(size=4))) +
      xlab("") + ylab("") +
      ggtitle(paste0('UMAP Gene Exp:',input$geneName))+
      theme_classic(base_size=14) 
      # theme(strip.background = element_blank(),
      #       strip.text.x     = element_blank(),
      #       axis.text.x      = element_blank(),
      #       axis.text.y      = element_blank(),
      #       axis.ticks       = element_blank(),
      #       axis.line        = element_blank(),
      #       panel.border     = element_blank())
    
  }
  
  
  
  
  #####################################

  ####################
  
  
  
  
  
  
  
 
  
  
  #####################################
  
  ######################################
  #
  # Multiple cluster heatmap
  #
  ####################################
  
  HeatmapRendering <- eventReactive(input$buttonForClusterHeatmap, {
    
    # validate(
    #   need(input$ChooseClusters != "All Clusters", "Please select a single valid cluster")
    # )
    
    
    
    # cdScFiltAnnotTmp <- cdScFiltAnnot
    # cdScFiltAnnotTmp$Cluster <- as.numeric(cdScFiltAnnotTmp$Clusters)
    
    #id <- c("1","2","3","4","5","6","7" ,"8","9","10","11")
    #id <- c("1","2","4","10","3","7","11" ,"6","8","9","5")
    id <- c(input$ChooseClustersHeatmap)
    # cdScFiltAnnotTmp <-  cdScFiltAnnotTmp[, colData(cdScFiltAnnotTmp)$Cluster %in% id]
    # cdScFiltAnnotTmp <-  cdScFiltAnnotTmp[,order(match(colData(cdScFiltAnnotTmp)$Cluster, id))]
    # 
    # 
    # #cdScFiltAnnotTmp <- cdScFiltAnnotTmp[,order(colData(cdScFiltAnnotTmp)$Cluster)]
    geneListTmp <- strsplit(input$ClusterGeneList, split = '\n')
    geneListTmp <- geneListTmp[[1]]
    # 
    # #chosenGene.exprs <- logcounts(cdScFiltAnnotTmp)[geneListTmp,colData(cdScFiltAnnotTmp)$Cluster %in% c(input$ChooseClustersHeatmap)]
    # chosenGene.exprs <- logcounts(cdScFiltAnnotTmp)[geneListTmp,]
    clustTemp <- as.numeric(cdScFiltAnnot$Clusters[cdScFiltAnnot$Clusters %in% id])
    names(clustTemp) <- cdScFiltAnnot$Barcode[cdScFiltAnnot$Clusters %in% id]
    clustTemp <- clustTemp[order(clustTemp, decreasing = FALSE)]
    
    chosenGene.exprs <- logcounts(cdScFiltAnnot)[geneListTmp,names(clustTemp)]
    
    heat.vals = t(apply(as.matrix(chosenGene.exprs), 1, function(x) {
      q10 = quantile(x, 0.1)
      q97 = quantile(x, 0.97)
      x[x < q10] = q10
      x[x > q97] = q97
      scale(x)
    }))
    
    #df = data.frame(Cluster = colData(cdScFiltAnnot)[colData(cdScFiltAnnot)$Clusters %in% c(input$ChooseClustersHeatmap),'Clusters'])
    df = data.frame(Clusters = as.factor(clustTemp), Sample = as.factor(cdScFiltAnnot[,names(clustTemp)]$Sample))
    clusterClass <- c_clust_col
    names(clusterClass) <- c(1:length(c_clust_col))
    
    sampleClass <- c_sample_col[1:nlevels(as.factor(cdScFiltAnnot$Sample))]
    names(sampleClass) <- levels(as.factor(cdScFiltAnnot$Sample))
    
    #ha = HeatmapAnnotation(df = df)
    ha = HeatmapAnnotation(df = df, col = list(Clusters = clusterClass[as.numeric(id)], Sample = sampleClass), 
                           annotation_legend_param = list(
                             Clusters = list(nrow = 1),
                             Sample=list(nrow = 1)))
    # rm(cdScFiltAnnotTmp)
    # gc()
    
    ht_list <- Heatmap(heat.vals,
            col = colorRamp2(c(-1.5,0,1.5), c(input$colminClustHeatmap, input$colmidClustHeatmap, input$colmaxClustHeatmap)),
            heatmap_legend_param = list(
              color_bar = "continuous",
              title = "Scaled expr",
              direction = "horizontal"
            ),        top_annotation = ha, show_column_names=FALSE, cluster_rows = FALSE, show_row_dend = FALSE,
            cluster_columns  = FALSE, show_column_dend = FALSE,
            row_names_gp = gpar(fontsize = 10))
    draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

    
  })
  
  
  
  
  output$plotClusterHeatmap <- renderPlot({
    
    HeatmapRendering()
    
  })
  
  
  
  ######################################
  #
  # Multiple sample heatmap
  #
  ####################################
  
  HeatmapRenderingSample <- eventReactive(input$buttonForSampleHeatmap, {
    
    # validate(
    #   need(input$ChooseClusters != "All Clusters", "Please select a single valid cluster")
    # )
    
    
    
    # cdScFiltAnnotTmp <- cdScFiltAnnot
    # cdScFiltAnnotTmp$Cluster <- as.numeric(cdScFiltAnnotTmp$Clusters)
    
    #id <- c("1","2","3","4","5","6","7" ,"8","9","10","11")
    #id <- c("1","2","4","10","3","7","11" ,"6","8","9","5")
    id <- c(input$ChooseSampleHeatmap)
    # cdScFiltAnnotTmp <-  cdScFiltAnnotTmp[, colData(cdScFiltAnnotTmp)$Cluster %in% id]
    # cdScFiltAnnotTmp <-  cdScFiltAnnotTmp[,order(match(colData(cdScFiltAnnotTmp)$Cluster, id))]
    # 
    # 
    # #cdScFiltAnnotTmp <- cdScFiltAnnotTmp[,order(colData(cdScFiltAnnotTmp)$Cluster)]
    geneListTmp <- strsplit(input$SampleGeneList, split = '\n')
    geneListTmp <- geneListTmp[[1]]
    # 
    # #chosenGene.exprs <- logcounts(cdScFiltAnnotTmp)[geneListTmp,colData(cdScFiltAnnotTmp)$Cluster %in% c(input$ChooseClustersHeatmap)]
    # chosenGene.exprs <- logcounts(cdScFiltAnnotTmp)[geneListTmp,]
    clustTemp <- as.numeric(cdScFiltAnnot$Clusters[cdScFiltAnnot$Sample %in% id])
    names(clustTemp) <- cdScFiltAnnot$Barcode[cdScFiltAnnot$Sample %in% id]
    clustTemp <- clustTemp[order(clustTemp, decreasing = FALSE)]
    
    sampleTemp <- cdScFiltAnnot[,names(clustTemp)]$Sample
    names(sampleTemp) <- cdScFiltAnnot[,names(clustTemp)]$Barcode
    
    
    #sampleTemp <- sampleTemp[order(sampleTemp, decreasing = FALSE)]
    
    chosenGene.exprs <- logcounts(cdScFiltAnnot)[geneListTmp,names(sampleTemp)]
    
    heat.vals = t(apply(as.matrix(chosenGene.exprs), 1, function(x) {
      q10 = quantile(x, 0.1)
      q97 = quantile(x, 0.97)
      x[x < q10] = q10
      x[x > q97] = q97
      scale(x)
    }))
    
    #df = data.frame(Cluster = colData(cdScFiltAnnot)[colData(cdScFiltAnnot)$Clusters %in% c(input$ChooseClustersHeatmap),'Clusters'])
    df = data.frame(Sample = as.factor(cdScFiltAnnot[,names(sampleTemp)]$Sample), Clusters = as.factor(clustTemp))
    clusterClass <- c_clust_col
    names(clusterClass) <- c(1:length(c_clust_col))
    
    sampleClass <- c_sample_col[1:nlevels(as.factor(cdScFiltAnnot$Sample))]
    names(sampleClass) <- levels(as.factor(cdScFiltAnnot$Sample))
    
    #ha = HeatmapAnnotation(df = df)
    ha = HeatmapAnnotation(df = df, col = list(Clusters = clusterClass[as.numeric(1:nlevels(cdScFiltAnnot$Clusters))], Sample = sampleClass), 
                           annotation_legend_param = list(
                             Clusters = list(nrow = 1),
                             Sample=list(nrow = 1)))
    # rm(cdScFiltAnnotTmp)
    # gc()
    
    ht_list <- Heatmap(heat.vals,
                       col = colorRamp2(c(-1.5,0,1.5), c(input$colminSampleHeatmap, input$colmidSampleHeatmap, input$colmaxSampleHeatmap)),
                       heatmap_legend_param = list(
                         color_bar = "continuous",
                         title = "Scaled expr",
                         direction = "horizontal"
                       ),        top_annotation = ha, show_column_names=FALSE, cluster_rows = FALSE, show_row_dend = FALSE,
                       cluster_columns  = FALSE, show_column_dend = FALSE,
                       row_names_gp = gpar(fontsize = 10))
    draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    
    
  })
  
  
  
  
  output$plotSampleHeatmap <- renderPlot({
    
    HeatmapRenderingSample()
    
  })
  
  
  
  
  ######################################
  #
  # Multiple sample & cluster heatmap
  #
  ####################################
  
  HeatmapRenderingSampleCluster <- eventReactive(input$buttonForSampleClusterHeatmap, {
    
    # validate(
    #   need(input$ChooseClusters != "All Clusters", "Please select a single valid cluster")
    # )
    
    
    
    # cdScFiltAnnotTmp <- cdScFiltAnnot
    # cdScFiltAnnotTmp$Cluster <- as.numeric(cdScFiltAnnotTmp$Clusters)
    
    #id <- c("1","2","3","4","5","6","7" ,"8","9","10","11")
    #id <- c("1","2","4","10","3","7","11" ,"6","8","9","5")
    idSample <- c(input$ChooseSampleClusterHeatmapSample)
    idCluster <- c(input$ChooseSampleClusterHeatmapCluster)
    # cdScFiltAnnotTmp <-  cdScFiltAnnotTmp[, colData(cdScFiltAnnotTmp)$Cluster %in% id]
    # cdScFiltAnnotTmp <-  cdScFiltAnnotTmp[,order(match(colData(cdScFiltAnnotTmp)$Cluster, id))]
    # 
    # 
    # #cdScFiltAnnotTmp <- cdScFiltAnnotTmp[,order(colData(cdScFiltAnnotTmp)$Cluster)]
    geneListTmp <- strsplit(input$SampleClusterGeneList, split = '\n')
    geneListTmp <- geneListTmp[[1]]
    # 
    # #chosenGene.exprs <- logcounts(cdScFiltAnnotTmp)[geneListTmp,colData(cdScFiltAnnotTmp)$Cluster %in% c(input$ChooseClustersHeatmap)]
    # chosenGene.exprs <- logcounts(cdScFiltAnnotTmp)[geneListTmp,]
    clustTemp <- as.numeric(cdScFiltAnnot$Clusters[cdScFiltAnnot$Sample %in% idSample &
                                                     cdScFiltAnnot$Clusters %in% idCluster])
    names(clustTemp) <- cdScFiltAnnot$Barcode[cdScFiltAnnot$Sample %in% idSample &
                                                cdScFiltAnnot$Clusters %in% idCluster]
    clustTemp <- clustTemp[order(clustTemp, decreasing = FALSE)]
    
    sampleTemp <- cdScFiltAnnot[,names(clustTemp)]$Sample
    names(sampleTemp) <- cdScFiltAnnot[,names(clustTemp)]$Barcode
    
    
    #sampleTemp <- sampleTemp[order(sampleTemp, decreasing = FALSE)]
    
    chosenGene.exprs <- logcounts(cdScFiltAnnot)[geneListTmp,names(sampleTemp)]
    
    heat.vals = t(apply(as.matrix(chosenGene.exprs), 1, function(x) {
      q10 = quantile(x, 0.1)
      q97 = quantile(x, 0.97)
      x[x < q10] = q10
      x[x > q97] = q97
      scale(x)
    }))
    
    #df = data.frame(Cluster = colData(cdScFiltAnnot)[colData(cdScFiltAnnot)$Clusters %in% c(input$ChooseClustersHeatmap),'Clusters'])
    df = data.frame(Sample = as.factor(cdScFiltAnnot[,names(sampleTemp)]$Sample), Clusters = as.factor(clustTemp))
    clusterClass <- c_clust_col
    names(clusterClass) <- c(1:length(c_clust_col))
    
    sampleClass <- c_sample_col[1:nlevels(as.factor(cdScFiltAnnot$Sample))]
    names(sampleClass) <- levels(as.factor(cdScFiltAnnot$Sample))
    
    #ha = HeatmapAnnotation(df = df)
    ha = HeatmapAnnotation(df = df, col = list(Clusters = clusterClass[as.numeric(1:nlevels(cdScFiltAnnot$Clusters))], Sample = sampleClass), 
                           annotation_legend_param = list(
                             Clusters = list(nrow = 1),
                             Sample=list(nrow = 1)))
    # rm(cdScFiltAnnotTmp)
    # gc()
    
    ht_list <- Heatmap(heat.vals,
                       col = colorRamp2(c(-1.5,0,1.5), c(input$colminSampleClusterHeatmap, input$colmidSampleClusterHeatmap, input$colmaxSampleClusterHeatmap)),
                       heatmap_legend_param = list(
                         color_bar = "continuous",
                         title = "Scaled expr",
                         direction = "horizontal"
                       ),        top_annotation = ha, show_column_names=FALSE, cluster_rows = FALSE, show_row_dend = FALSE,
                       cluster_columns  = FALSE, show_column_dend = FALSE,
                       row_names_gp = gpar(fontsize = 10))
    draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    
    
  })
  
  
  
  
  output$plotSampleClusterHeatmap<- renderPlot({
    
    HeatmapRenderingSampleCluster()
    
  })
  
  ##################################
  
  
  
 
  
 
  
  
  #####################################
  #
  # DE between two clusters
  #
  ################################
  
  
  output$AllClustTwoClustComp <- renderPlotly({
    
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    df[,'Clusters'] <- as.factor(my.clusters)
    df %>% plot_ly(color = ~Clusters, colors = c_clust_col[c(1:nlevels(cdScFiltAnnot$Clusters))], type="scatter", mode="markers", hoverinfo = 'text',
                   text = ~paste('</br> Clusters: ', Clusters,
                                 '</br> Samples: ', Sample)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2,
                marker = list(
                  size = 3)) %>% 
      layout(title = 'tSNE with Clusters',showlegend = TRUE, legend = list(title='Cluster',font = list(size = 10), itemsizing='constant'))
    
  })
  
  
  output$plotSelectClust <- renderPlot({
    #df$Group <- "All Group"
    
    df_shiny_ForTwoClust[colData(cdScFiltAnnot)$Clusters %in% input$clust1, "Group"] <- "blue"
    df_shiny_ForTwoClust[colData(cdScFiltAnnot)$Clusters %in% input$clust2, "Group"] <- "red"
    
    
    
    p <- ggplot(df_shiny_ForTwoClust, aes(V1, V2, col = I(Group))) +
      geom_point(size=0.5)  +
      theme_classic(base_size=10) +
      theme(strip.background = element_blank(),
            strip.text.x     = element_blank(),
            axis.text.x      = element_blank(),
            axis.text.y      = element_blank(),
            axis.ticks       = element_blank(),
            axis.line        = element_blank(),
            panel.border     = element_blank())
    
    p
  })
  
  
  tableRenderingTwoClust <- eventReactive(input$buttonForTwoClustDE, {
    
    print(input$clust1)
    print(input$clust2)
    
    counts1 <- as.matrix(counts(cdScFiltAnnot)[,colData(cdScFiltAnnot)$Clusters %in% input$clust1])
    counts2 <- as.matrix(counts(cdScFiltAnnot)[,colData(cdScFiltAnnot)$Clusters %in% input$clust2])
    
    countsTable <- cbind(counts1,counts2)
    rownames(countsTable) <- make.names(rownames(countsTable), unique = TRUE)
    
    conds.Sel <- c(rep("1", dim(counts1)[2]), rep("2", dim(counts2)[2]))
    print(table(conds.Sel))
    
    res1 <- nbTestSH(countsTable, conds.Sel, "1", "2")
    
    normCounts1 <-  as.matrix(logcounts(cdScFiltAnnot)[,colData(cdScFiltAnnot)$Clusters %in% input$clust1])
    normCounts2 <-  as.matrix(logcounts(cdScFiltAnnot)[,colData(cdScFiltAnnot)$Clusters %in% input$clust2])
    res1[,'FDR'] <- p.adjust(res1[,7], method = "bonferroni")
    
    res1$Mean <- round(res1$Mean,3)
    res1$rawMeanA <- round(res1$rawMeanA,3)
    res1$rawMeanB <- round(res1$rawMeanB,3)
    res1$rawLog2FoldChange <- round(res1$rawLog2FoldChange,3)
    
    res1$pval <- formatC(res1$pval, format = "E", digits = 2)
    res1$FDR <- formatC(res1$FDR, format = "E", digits = 2)
    #res1[,'Z-score_1'] <- colMeans(scale(t(normCounts1)))
    #res1[,'Z-score_2'] <- colMeans(scale(t(normCounts2)))
    res1[,c('Mean','rawMeanA','rawMeanB','rawLog2FoldChange','pval','FDR')]
  }
  )
  
  
  output$mytableTwoClust = DT::renderDataTable(tableRenderingTwoClust(), server = FALSE, extensions = 'Buttons',
                                               
                                               options = list(dom = 'lBfrtip',
                                                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                                              pageLength = 5, autoWidth = TRUE)
  )
  
  
  ####################################
 
  observe({
    #choices = c('Kmeans','Graph Based Cluster','Dynamic Tree Cut','Default Cluster'),
    clustType <- input$chooseClustType
    if(clustType == 'Kmeans')
    {
      clust <- list_cluster
    }
    else if(clustType == 'Graph Based Cluster'){
      
      clust <- list_cluster
    }
    else if(clustType == 'Dynamic Tree Cut')
    {
      clust <- list_cluster
      
    }
    else{
      clust <- list_cluster
    }
    
    kmeansClusterRes <- input$chooseClust
    
    #Can also set the label and select items
    updateSelectizeInput(session, "clust1forK",
                      label = "Select Cluster",
                      choices = as.factor(sort(as.numeric(clust[[kmeansClusterRes]][order(list_cluster[[kmeansClusterRes]])]))),
                      )
    
    updateSelectizeInput(session, "clust2forK",
                         label = "Select Cluster",
                         choices = as.factor(sort(as.numeric(list_cluster[[kmeansClusterRes]][order(list_cluster[[kmeansClusterRes]])]))),
                         )
    
  })
  
  tableRenderingTwoClustKmeans <- eventReactive(input$buttonForTwoClustDEKmeans, {
    
    print(input$clust1forK)
    print(input$clust2forK)
    
    counts1 <- as.matrix(counts(cdScFiltAnnot)[,colData(cdScFiltAnnot)$kmeansCluster %in% input$clust1forK])
    counts2 <- as.matrix(counts(cdScFiltAnnot)[,colData(cdScFiltAnnot)$kmeansCluster %in% input$clust2forK])
    
    countsTable <- cbind(counts1,counts2)
    rownames(countsTable) <- make.names(rownames(countsTable), unique = TRUE)
    
    conds.Sel <- c(rep("1", dim(counts1)[2]), rep("2", dim(counts2)[2]))
    print(table(conds.Sel))
    
    res1 <- nbTestSH(countsTable, conds.Sel, "1", "2")
    
    #normCounts1 <-  as.matrix(logcounts(cdScFiltAnnot)[,colData(cdScFiltAnnot)$Clusters %in% input$clust1])
    #normCounts2 <-  as.matrix(logcounts(cdScFiltAnnot)[,colData(cdScFiltAnnot)$Clusters %in% input$clust2])
    res1[,'FDR'] <- p.adjust(res1[,7], method = "bonferroni")
    
    res1$Mean <- round(res1$Mean,3)
    res1$rawMeanA <- round(res1$rawMeanA,3)
    res1$rawMeanB <- round(res1$rawMeanB,3)
    res1$rawLog2FoldChange <- round(res1$rawLog2FoldChange,3)
    
    res1$pval <- formatC(res1$pval, format = "E", digits = 2)
    res1$FDR <- formatC(res1$FDR, format = "E", digits = 2)
    #res1[,'Z-score_1'] <- colMeans(scale(t(normCounts1)))
    #res1[,'Z-score_2'] <- colMeans(scale(t(normCounts2)))
    res1[,c('Mean','rawMeanA','rawMeanB','rawLog2FoldChange','pval','FDR')]
  }
  )
  
  
  output$mytableTwoClustKmeans = DT::renderDataTable(tableRenderingTwoClustKmeans(), server = FALSE, extensions = 'Buttons',
                                               
                                               options = list(dom = 'lBfrtip',
                                                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                                              pageLength = 5, autoWidth = TRUE)
  )
  
  
  ####################################
  
  
  
  

  
  
  ####################################
  
  
  
  ##########################################
  
  
  
  
  ####################################
  #
  # All info messages
  #
  ####################################
  
  observeEvent(input$projectionInfo, {

    showModal(modalDialog(
      title = "Projection",
      tags$div(
        tags$ul(
          tags$li("Projection of cells into 2D space. The projection can be toggled between tSNE and UMAP from the Projection dropdown menu"),
          tags$li("Cells can be coloured based on samples, clusters, number of UMIs each cell have and number of genes each cell express"),
          tags$li("Initial size of the dots are 1.5 which can be changed by moving the slide bar"),
          tags$li("Dot opacity controls the visual transperancy of the cells. It can be lowered by moving the slide bar so that cells 
                  underneath one eather can be can be viewed"),
          tags$li("Samples or clusters can be removed by clicking on each of the lables in the plot. One click will remove the cells associated with 
                  that feature and double click will remove all but the one that is clicked"),
          tags$li("Hovering mouse over each cell would show the cell barcode, cluster, sample and CellType"),
          tags$li("If you want to see a specific cell type, just double click on the legend of the image, so other cell
                  types would be deactivated from visualization."),
          tags$li("You can select/deselect samples or clusters from the dropdown list bottom of the figure. This would show only the selected cells.")
        )
      ),
      
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$sampleByclustersInfo, {
    
    showModal(modalDialog(
      title = "Samples into clusters",
      tags$div(
        "This table shows the number of samples distributed into different clusters. It also shows the total number of cells in a sample."
        ),
      
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  
  
  observeEvent(input$percentClusterInSampleInfo, {
    
    showModal(modalDialog(
      title = "Percent of cells from cluster going to sample",
      tags$div(
        "This figures shows the % of cells going into a cluster. All the cells in a specific sample are divided into different clusters and the percentage
        is calcualted based on cells only this sample."
      ),
      
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  
  observeEvent(input$umiinSampleInfo, {
    
    showModal(modalDialog(
      title = "UMI in sample",
      tags$div(
        "This figures shows the total number of Unique Molecular Identifier (UMI)s for each individul sample in log10 scale. Unique molecular identifiers (UMI) are molecular tags that are used to detect and 
         quantify unique mRNA transcripts. In this method, mRNA libraries are generated by fragmentation and reverse-transcribed to cDNA. 
         Oligo(dT) primers with specific sequencing linkers are added to the cDNA."
      ),
      
      easyClose = TRUE,
      footer = NULL
      ))
  })
  
  observeEvent(input$numberOfGenesExpressedInfo, {
    
    showModal(modalDialog(
      title = "Total genes expressed in sample",
      tags$div(
        "This violin plot shows the total number of genes expressed in each sample, i.e. genes having a UMI count > 0"
      ),
      
      easyClose = TRUE,
      footer = NULL
      ))
  })
  
  observeEvent(input$umiInClustersInfo, {
    
    showModal(modalDialog(
      title = "Total genes expressed in sample",
      tags$div(
        "This figures shows the total number of Unique Molecular Identifier (UMI)s for each individul cluster in log10 scale. Unique molecular identifiers (UMI) are molecular tags that are used to detect and 
         quantify unique mRNA transcripts. In this method, mRNA libraries are generated by fragmentation and reverse-transcribed to cDNA. 
        Oligo(dT) primers with specific sequencing linkers are added to the cDNA. "
      ),
      
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$numberOfGenesExpressedInClustersInfo, {
    
    showModal(modalDialog(
      title = "Total genes expressed in sample",
      tags$div(
        "This violin plot shows the total number of genes expressed in each cluster, i.e. genes having a UMI count > 0"
      ),
      
      easyClose = TRUE,
      footer = NULL
      ))
  })
  
  observeEvent(input$DE_between_sample_and_clustersSelectionInfo, {
    
    showModal(modalDialog(
      title = "Differential Expression between multiple sample and cluster",
      tags$div(
        "This input panel allows you to choose samples and clusters within a sample. DE will be calculated on cells that
        belong to the selected samples and clusters. This panel helps to select the cells that belong to a specific sample
        and cluster. For eg. cluster 1 is shared between Sample 1 and Sample 2. If you select cluster 1 and sample 1 then
        cells only in cluster 1 and sample 1 would be selected for the DE.",
          tags$ul(
            tags$li("The first selection is for the first group of cells."),
            tags$li("Second selection is for the second group of cells.")
          )
      ),
      
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  
  observeEvent(input$DE_between_sample_and_clustersProjectionInfo, {
    showModal(modalDialog(
      title = "Projection of cells",
      tags$div(
        "This panel projects the cells into a 2D dimension. Users can choose the projection of tSNE or UMAP. 
        They can also choose to colour the cells based on sample or cluster. This panel helps the users to identify
        the cells they want to choose for differential expression."
      ),
      
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$DE_in_manual_selection_input, {
    showModal(modalDialog(
      title = "Manual selection of cells for DE",
      tags$div(
        tags$ul(
          tags$li("This input panel configures the projection panel based on user selection of the inputs."),
          tags$li("In order to select cells manually on which the users want to run DE, they first need to select the first group of cells
                  by dragging with mouse pointer over the cells. The second group of cells are then selected by dragging again on a 
                  second group of cells. The first group of cells would then be coloured blue and the second group would be coloured blue."),
          tags$li("The user would then click DE between selected cells to run the Differential Expression analysis between the 
                  two groups of cells"),
          tags$li("In order to clear all the selection the uer needs to click Reset all selection")
          )
        
      ),
      
      easyClose = TRUE,
      footer = NULL
      ))
  })
  
  
  
  
  ###################################
  
  
  output$Cluster_KMEANS <- renderPlotly({
    
    gc()
    
    NormData<-assays(cdScFiltAnnot)$logcounts
    
    NormDataTransposed <- t(NormData)
    
    km.res10 <- kmeans(NormDataTransposed, algorithm = "Lloyd", 10, nstart = 15)

    cdScFiltAnnot$kmeansCluster = as.factor(km.res10$cluster)
    
    #saveHDF5SummarizedExperiment(cdScFiltAnnot, dir="C:/Users/hee yaa/Documents/cdScFiltAnnotHDF5", prefix="", replace = TRUE,
     #                            chunkdim=NULL, level=NULL, verbose=FALSE)
    
    p <- fviz_cluster(km.res10, data = NormDataTransposed , repel=TRUE,
                      ellipse.type = "convex") # save to access $data
    
    # save '$data'
    data <- p$data 
    
    # calculate the convex hull using chull(), for each cluster
    hull_data <-  data %>%
      group_by(cluster) %>%
      slice(chull(x, y))
    
     #saveRDS(hull_data, file = "dirName/kmRES10.Rds")
     #hull_data2 <- readRDS(file = "dirName/kmRES10.Rds")

   
      #View(data$name)
      # plot: you can now customize this by using ggplot sintax
      
      GGplot <- ggplot(data, aes(x, y)) + geom_point() +
        #geom_polygon(data = hull_data, alpha = 0.2,lwd=1, aes(color = cluster, linetype = cluster))
        geom_polygon(data = hull_data, alpha = 0.5, aes(fill=cluster, linetype=cluster))
      #GGplotmod <- ggplot(data, aes(x, y, text = name)) + geom_point() +
      # geom_polygon(data = hull_data, alpha = 0.5, aes(fill=cluster, linetype=cluster))
      
      #geom_point(aes(text = name), colour = "red", alpha = 1/2)
      #ggplotly(GGplot) %>% layout(dragmode = "lasso")
      ggplotly(GGplot) %>% layout(dragmode = "lasso")
   
    
    
  })
  
  RenderKmeansClusterPlot <-  reactive({
  
    #eventReactive(input$buttonForKmeans,
    #renderPlotly({
    # size of K
    
    k <- input$choice_K
    
    objectName <- paste0('kmeans.res',k)
    print(objectName)
    
    list_cluster <- readRDS("Cluster.RDS")

    if(is.null(list_cluster[[objectName]]))
    {
      gc()
      
      NormData<-assays(cdScFiltAnnot)$logcounts
      
      NormDataTransposed <- t(NormData)
      
      km.res <- kmeans(NormDataTransposed, algorithm = "Lloyd", k)
      
      list_cluster[[objectName]] <- km.res$cluster
      
      saveRDS(list_cluster, file = "Cluster.RDS")
      
      km.res <- list_cluster[[objectName]]
      
      }
    else{
      
      km.res <- list_cluster[[objectName]]
    }
    
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    df[,'Cell']=as.factor(colData(cdScFiltAnnot)$Barcode)
    df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    df[,'Clusters'] <- as.factor(km.res)
    df[,'cellType'] <- as.factor(colData(cdScFiltAnnot)$cellType)
    df %>%
      plot_ly(color = ~Clusters, colors = c_clust_col[c(1:9)], type="scatter", mode="markers", hoverinfo = 'text',
              #marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Sample,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>%
      layout(title = 'Kmeans with tSNE',showlegend = TRUE, legend = list(font = list(size = 10), itemsizing='constant'))

  
  })
  
  output$Cluster_KMEANSwithtSNE <- renderPlotly({
  
    RenderKmeansClusterPlot() 

    }) 
  
  
  output$GaphBased_Graph <- renderPlot({
   
    NormDataP<-assays(cdScFiltAnnot)$logcounts
    
    SNNGRAPH <- buildSNNGraph(NormDataP, k = 30)
    
    ModifiedGraph <- igraph::as_data_frame(SNNGRAPH, what="edges")
    
    GRAPHToLouvain <- graph_from_data_frame(ModifiedGraph, directed = FALSE)
    
    Louvain_GRAPH <- cluster_louvain(GRAPHToLouvain)
    
    membership(Louvain_GRAPH)
    communities(Louvain_GRAPH)
    
    plot(Louvain_GRAPH,GRAPHToLouvain)
    
  })
  
  
  output$KNN <- renderPlot({
    
    NormDataP<-assays(cdScFiltAnnot)$logcounts
    
    transposed <- t(NormDataP)
    
    KNNGraph <- FindNeighbors(transposed, dims = 1:30, k.param = 60, prune.SNN = 1/15)
    
    #pheatmap::pheatmap(KNNGraph$nn, border_color = "grey90", 
     #                  cluster_rows = F, cluster_cols = F, fontsize = 5)
    
    
    pheatmap(KNNGraph$nn, col = c("white", "black"), border_color = "grey90", 
             cluster_rows = F, cluster_cols = F, fontsize = 2)
  })
  
  
  output$SNN <- renderPlot({
    
    NormDataP<-assays(cdScFiltAnnot)$logcounts
    
    SNNGRAPH <- buildSNNGraph(NormDataP, k = 30)
    
    
    plot.igraph(SNNGRAPH)
    
  })
  
  
  
  
  output$GraphBasedLouvainCluster <- renderPlotly({
    
    k<- input$choice_K_for_Louvain
    
    NormDataP<-assays(cdScFiltAnnot)$logcounts
    SNNGRAPH <- buildSNNGraph(NormDataP, k = k)
    
    ModifiedGraph <- igraph::as_data_frame(SNNGRAPH, what="edges")
    
    GRAPHToLouvain <- graph_from_data_frame(ModifiedGraph, directed = FALSE)
    
    Louvain_GRAPH <- cluster_louvain(GRAPHToLouvain)
    
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    df[,'Cell']=as.factor(colData(cdScFiltAnnot)$Barcode)
    df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    df[,'Clusters'] <- as.factor(Louvain_GRAPH$membership)
    df[,'cellType'] <- as.factor(colData(cdScFiltAnnot)$cellType)
    df %>%
      plot_ly(color = ~Clusters, colors = c_clust_col[c(1:9)], type="scatter", mode="markers", hoverinfo = 'text',
              #marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Sample,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>%
      layout(title = 'tSNE with Graph Based Clusters',showlegend = TRUE, legend = list(font = list(size = 10), itemsizing='constant'))
    #fviz_silhouette(res.hc)
    #p <- plot.igraph(graphSNN10)
    #plotly(p)
  })
  
  output$Dynamic_Tree_Cut_Plot <- renderPlotly({
    
    NormData<-assays(cdScFiltAnnot)$logcounts
    
    NormDataTransposed <- t(NormData)
    ddRead <- readRDS(file = "Distance.Rds")
    
    if(is.null(ddRead)){
      
      ddRead <- dist(NormDataTransposed, method = "euclidean")
      saveRDS(ddRead, file = "Distance.Rds")
    }
    
    hcRead <- readRDS(file = "hclust.Rds")
    
    if(is.null(hcRead)){
      
      hcRead <- hclust(ddRead, method = "ward.D2")
      saveRDS(hcRead, file = "hclust.Rds")
    }
    
    Cluster.DTC <- unname(cutreeDynamic(hcRead, distM=as.matrix(ddRead), verbose=1))
    
    cdScFiltAnnot$Cluster.DTC <- as.factor(Cluster.DTC)
    
    #saveHDF5SummarizedExperiment(cdScFiltAnnot, dir="cdScFiltAnnotHDF5", prefix="", replace = TRUE,
     #                            chunkdim=NULL, level=NULL, verbose=FALSE)
    
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    df[,'Cell']=as.factor(colData(cdScFiltAnnot)$Barcode)
    df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    df[,'Clusters'] <- as.factor(Cluster.DTC)
    df[,'cellType'] <- as.factor(colData(cdScFiltAnnot)$cellType)
    df %>%
      plot_ly(color = ~Clusters, colors = c_clust_col[c(1:9)], type="scatter", mode="markers", hoverinfo = 'text',
              #marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Sample,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>%
      layout(title = 'tSNE with Dynamic Tree Cut',showlegend = TRUE, legend = list(font = list(size = 10), itemsizing='constant'))
    #fviz_silhouette(res.hc)
    #p <- plot.igraph(graphSNN10)
    #plotly(p)
  })
  
  
  ####################################
  #
  # Marker genes for Cluster
  #
  ####################################
  
  
  
  observe({
    
    Cluster_Res <- input$ClusterRes
    
    clustType <- input$chooseClustTypeforMarker
    
    if(clustType == 'Kmeans')
    {
      clust <- list_cluster
    }
    else if(clustType == 'Graph Based Cluster'){

      clust <- list_cluster
    }
    else if(clustType == 'Dynamic Tree Cut')
    {
      clust <- list_cluster

    }
    else{
      clust <- list_cluster
    }

    #kmeansClusterRes <- input$chooseClust
    
    #Can also set the label and select items
    updateSelectizeInput(session, "markerChoose_Kmeans_Cluster",
                         label = "Select Cluster",
                         choices = as.factor(sort(as.numeric(clust[[Cluster_Res]][order(clust[[Cluster_Res]])]))),
    )
  })
  
  
  
  tableRenderingmarkerTableCluster <- function(){
    
    
    clustID <- input$markerChooseCluster
    dfColName <- paste0('PercentClust',clustID)
    df <- data.frame("Cluster"=clustID,
                     "Gene"= rownames(cdScFiltAnnot),
                     "FDR" = metadata(cdScFiltAnnot)[['Cluster']][[1]][[clustID]][rownames(cdScFiltAnnot),'FDR'],
                     "PercentInClust" = rowData(cdScFiltAnnot)[,dfColName]/100)
    df <- as.data.frame(cbind(df, metadata(cdScFiltAnnot)[['Cluster']][[1]][[clustID]][rownames(cdScFiltAnnot),3:(max(as.numeric(cdScFiltAnnot$Clusters))+1)]))
    df <- df[order(df$FDR, decreasing = FALSE),]
    df$FDR <- formatC(df$FDR, format = "E", digits = 2)
    
    #df <- df[1:100,]
    #rownames(df) <- NULL
    df %>%
      datatable(colnames = c('% cells in this cluster' = 'PercentInClust'), 
                rownames = FALSE,
                caption = 'Table : Marker genes for cluster.') %>%
      #formatRound(columns = 'FDR', digits = 3) %>%
      formatRound(columns=c(colnames(df)[grep('logFC',colnames(df))]), digits=3) %>%
      formatPercentage(columns = '% cells in this cluster') %>% 
      formatStyle(names(df[,5:dim(df)[2]]),
                  background = styleColorBar(range(df[,5:dim(df)[2]]), 'lightblue'),
                  backgroundSize = '98% 58%',
                  backgroundRepeat = 'no-repeat',
                  backgroundPosition = 'center')
    
    
    
  }
  
  
  
  output$markerTableCluster = DT::renderDataTable(tableRenderingmarkerTableCluster (), server = FALSE, extensions = 'Buttons', 
                                                  options = list(dom = 'lBfrtip',
                                                                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
                                                                 pageLength = 5, autoWidth = TRUE))
  
  tableRenderingmarkerTableKCluster <- reactive({

    Cluster_Res <- input$ClusterRes
    
    objectName <- paste0('marker_genes_of_',Cluster_Res)
    
    list_marker_gene <- readRDS(file = "List_of_Marker_Genes.RDS")
      
    if(is.null(list_marker_gene[[objectName]])){
     
       marker.genes.cluster.up <- findMarkers(cdScFiltAnnot, groups= list_cluster[[Cluster_Res]], pval.type="all", lfc=0.5, direction="up")
     
       list_marker_gene[[objectName]] <- list(marker.genes.cluster.up) 
       
       cdScFiltAnnot[[objectName]] <- list(marker.genes.cluster.up)
       
       saveRDS(list_marker_gene, file = "List_of_Marker_Genes.RDS")
       
       cdScFiltAnnot[[Cluster_Res]]  <- list_cluster[[Cluster_Res]]
    }
    else{
      
      metadata(cdScFiltAnnot)[[objectName]] <- list_marker_gene[[objectName]]
      cdScFiltAnnot[[Cluster_Res]]  <- list_cluster[[Cluster_Res]]
      
    }

    
    clustID <- input$markerChoose_Kmeans_Cluster
    dfColName <- paste0('PercentClust',clustID)
    df <- data.frame("Cluster"=clustID,
                     "Gene"= rownames(cdScFiltAnnot),
                     "FDR" = metadata(cdScFiltAnnot)[[objectName]][[1]][[clustID]][rownames(cdScFiltAnnot),'FDR'],
                     "PercentInClust" = rowData(cdScFiltAnnot)[,dfColName]/100
                     )
    df <- as.data.frame(cbind(df, metadata(cdScFiltAnnot)[[objectName]][[1]][[clustID]][rownames(cdScFiltAnnot),3:(max(as.numeric(cdScFiltAnnot[[Cluster_Res]]))+1)]))
    df <- df[order(df$FDR, decreasing = FALSE),]
    df$FDR <- formatC(df$FDR, format = "E", digits = 2)
    
    #df <- df[1:100,]
    #rownames(df) <- NULL
    df %>%
      datatable(colnames = c('% cells in this cluster' = 'PercentInClust'), 
                rownames = FALSE,
                caption = 'Table : Marker genes for cluster.') %>%
      #formatRound(columns = 'FDR', digits = 3) %>%
      formatRound(columns=c(colnames(df)[grep('logFC',colnames(df))]), digits=3) %>%
      formatPercentage(columns = '% cells in this cluster') %>% 
      formatStyle(names(df[,5:dim(df)[2]]),
                  background = styleColorBar(range(df[,5:dim(df)[2]]), 'lightblue'),
                  backgroundSize = '98% 58%',
                  backgroundRepeat = 'no-repeat',
                  backgroundPosition = 'center')
    
  })
  
  
  
  output$markerTableKmeansCluster = DT::renderDataTable(tableRenderingmarkerTableKCluster (), server = FALSE, extensions = 'Buttons', 
                                                  options = list(dom = 'lBfrtip',
                                                                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
                                                                 pageLength = 5, autoWidth = TRUE))

  
  ####################################
  #
  # Enriched Pathway for cluster
  #
  ####################################
  
  observe({

    clusterID <- as.numeric(input$ENChooseCluster)
    
    # Can also set the label and select items
    updateSelectInput(session, "ENChooseClusterGO",
                      label = "Choose cluster enrichr",
                      choices = names(which(sapply(metadata(cdScFiltAnnot)[['enrichCluster']][[1]][[clusterID]], NROW) != 0)),
                      selected = names(which(sapply(metadata(cdScFiltAnnot)[['enrichCluster']][[1]][[clusterID]], NROW) != 0))[1])
    
    kmeansclusterID <- as.numeric(input$ENChooseKmeansCluster)
    
    #Can also set the label and select items
    Enrichment_Kmeans <- readRDS(file = "Enrichment_Kmeans.res.10.rds")
    updateSelectInput(session, "ENChooseKmeansClusterGO",
                     label = "Choose cluster enrichr",
                   choices = names(which(sapply(Enrichment_Kmeans[[1]][[kmeansclusterID]], NROW) != 0)),
                 selected = names(which(sapply(Enrichment_Kmeans[[1]][[kmeansclusterID]], NROW) != 0))[1])
     
  })
  
  
  tableRenderingENTableCluster <- function(){
    
    
    clusterID <- as.numeric(input$ENChooseCluster)
    clusterEN <- input$ENChooseClusterGO
    df <- metadata(cdScFiltAnnot)[['enrichCluster']][[1]][[clusterID]][[clusterEN]]  
    if(df$P.value <0.001)
      df$P.value <- formatC(df$P.value, format = "E", digits = 2)
    if(df$Adjusted.P.value < 0.01)
      df$Adjusted.P.value <- formatC(df$Adjusted.P.value, format = "E", digits = 2)
    df
  }
  
  
  
  output$ENTableCluster = DT::renderDataTable(tableRenderingENTableCluster(), server = FALSE, extensions = 'Buttons', 
                                              options = list(dom = 'lBfrtip',
                                                             buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
                                                             pageLength = 5, autoWidth = TRUE))
  
  
  tableRenderingENTableKCluster <- function(){
    
    #Enrichment_Kmeans.res.10 <- readRDS(file = "Enrichment_Kmeans.res.10.rds")
    check = TRUE;
    if(!check){
      
      dbs <- listEnrichrDbs()
      
      dbsSel <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018", "Enrichr_Submissions_TF-Gene_Coocurrence",
                  "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO", "KEGG_2015", "MGI_Mammalian_Phenotype_2013", "Human_Gene_Atlas","Mouse_Gene_Atlas")
      
      markerGene <- metadata(cdScFiltAnnot)[['kmeansCluster.res10']][[1]]
      
      resEnrichClust <- list()
      for(i in 1:max(as.numeric(cdScFiltAnnot$kmeansCluster))){
        
        resEnrich <- enrichr(rownames(markerGene[[i]][markerGene[[i]]$FDR<=0.01,]), databases = dbsSel)
        resEnrichTemp <- lapply(resEnrich, FDRsubsetting <- function(x){
          x <- x[x$Adjusted.P.value<0.1,c(1:4)]
          return(x)
        })
        resEnrichClust[[i]] <- resEnrichTemp
      }
      
      enrichment_Kmeans.res.10 <- list(resEnrichClust)
      saveRDS(enrichment_Kmeans.res.10, file = "Enrichment_Kmeans.res.10.rds")
     
      
      #metadata(cdScFiltAnnot)[['enrichKmeansCluster10']] <- list(resEnrichClust)
      
      #saveHDF5SummarizedExperiment(cdScFiltAnnot, dir="C:/Users/hee yaa/Documents/cdScFiltAnnotHDF5", prefix="", replace = TRUE,
       #                            chunkdim=NULL, level=NULL, verbose=FALSE)
      
    }
    
    Enrichment_Kmeans.res.10 <- readRDS(file = "Enrichment_Kmeans.res.10.rds")
    
    clusterID <- as.numeric(input$ENChooseKmeansCluster)
    clusterEN <- input$ENChooseKmeansClusterGO
    df <-Enrichment_Kmeans.res.10[[1]][[clusterID]][[clusterEN]]  
    if(df$P.value <0.001)
      df$P.value <- formatC(df$P.value, format = "E", digits = 2)
    if(df$Adjusted.P.value < 0.01)
      df$Adjusted.P.value <- formatC(df$Adjusted.P.value, format = "E", digits = 2)
    df
  }
  
  
  
  output$ENTableKmeansCluster = DT::renderDataTable(tableRenderingENTableKCluster(), server = FALSE, extensions = 'Buttons', 
                                              options = list(dom = 'lBfrtip',
                                                             buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
                                                             pageLength = 5, autoWidth = TRUE))
  
  
  
  ####################################
  
  
  }

shinyApp(ui, server)