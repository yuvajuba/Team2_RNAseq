#   Setting     ###################
# =============================== #
library(dplyr)
library(ggplot2)
library(tibble)
library(writexl)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ComplexHeatmap)
library(DT)
library(stringr)
library(shiny)
library(shinythemes)
library(shinyBS)
library(bslib)
library(circlize)
library(purrr)
library(ggvenn)
library(forcats)
library(ggupset)
library(openxlsx)

MyPalette <- c("#9933aa","#ffdd22","#aa4400","#ff0000","#337722","#00ff66","#005566","#002277",
               "#441144","#aa0077","#00bbff","#003333","#4422cc","#116611","#330077","#111111",
               "#667700","#ddaa00","#33ffff","#ff22ff","#ffff33","#00ff00","#0000ff","#444444")


#   ShinyApp options   ##################
# ====================================== #
options(shiny.maxRequestSize = 900*1024^2)


#  Fun      #######################
# =============================== #
Apply_filt <- function(obj, fc, p.adj, gene_type = NULL){
  #  -obj:        dataframe class object
  #  -fc:         num: min threshold to keep for log2FoldChange
  #  -padj:       num: max threshold to keep for p-value adjusted
  #  -gene_type:  str: if provided it would be the gene_type to keep
  
  mydata <- obj
  
  if(is.null(gene_type)){
    dat <- mydata %>% dplyr::filter(padj < p.adj,
                                    abs(log2FoldChange) > fc)
  } else {
    dat <- mydata %>% dplyr::filter(padj < p.adj,
                                    abs(log2FoldChange) > fc,
                                    Gene_type %in% gene_type)
  }
  
  return(dat)
}

#  Shiny ui    ####################
# ================================= #
ui <- fluidPage(
  theme = shinytheme("flatly"),
  
  
  ## I- Import data   ##########################
  h2("I- Importing the data",
     style = "color:gold ; font-weight:600 ; background-color:navy ; margin-top:10px ; margin-bottom:30px"),
  
  fluidRow(
    column(width = 4,
           fileInput(inputId = "Load",
                     label = "Load csv or xlsx file",
                     accept = c(".csv",".xlsx"),
                     width = "90%")),
    column(width = 8,
           dataTableOutput("tab_data"))
  ),
  
  
  
  ## II- Build marker list  ##################
  h2("II- Building marker list",
     style = "color:gold ; font-weight:600 ; background-color:navy ; margin-top:50px ; margin-bottom:20px"),
  
  fluidRow(
    column(width = 5,
           div(
             style = "margin-bottom:20px ; text-align:left",
             p("Filtering parammeters for the DEGs",
               style="color:darkgreen; font-size:120%; font-weight:600; margin-bottom:20px"),
             numericInput(inputId = "logfc_threshold",
                          label = "Log2FC threshold",
                          value = 1,
                          step = 0.5,
                          width = "30%"),
             numericInput(inputId = "padj_threshold",
                          label = "p-value threshold",
                          value = 0.05,
                          step = 0.01,
                          width = "30%"),
             selectInput(inputId = "gene_type_filt",
                         label = "Specify Gene_type",
                         multiple = T,
                         selectize = F,
                         choices = c(),
                         size = 8,
                         width = "50%")
           ),
           
           actionButton("ApplyFilter", "Apply filtering", width = "50%",
                        style = "background-color:darkgreen; color:white; font-weight:600")),
    column(width = 7,
           div(
             style = "text-align:left; margin-bottom:20px",
             p("Overview of the list :",
               style="font-size:120%; color:purple; font-weight:600; margin-bottom:20px"),
             tags$pre(
               style = "color: black;
                background-color: white;
                max-height: 500px;
                overflow-y: auto;
                font-size: 14px;
                border: .5px solid #fff;
                padding: 0px;",
               verbatimTextOutput("list_summary")
             )
             ))
  ),
  
  
  
  fluidRow(
    p("And below are the results after filtering :",
      style="font-size:130%; color:purple; font-weight:600; margin-bottom:20px; margin-top:30px"),
    dataTableOutput("tab_data_filt")
  ),
  
  
  
  
  ## III- PEA  ##################
  h2("III- Pathway Enrichment Analysis",
     style = "color:gold ; font-weight:600 ; background-color:navy ; margin-top:120px ; margin-bottom:20px"),
  
  ### 1- GO analysis      ###################
  h3("1- GO analysis",
     style = "color:darkred ; font-weight:600 ; background-color:gold ; margin-top:10px ; margin-bottom:20px"),
  
  fluidRow(
    column(width = 3,
           selectInput("go.list.genes", "Select gene list", 
                       choices = c("All DEGs","Up Genes","Down Genes"), 
                       width = "250px")),
    column(width = 3,
           selectInput("go.keytype", "Key Type", choices = c("SYMBOL","ENTREZID","ENSEMBL"), width = "250px")),
    column(width = 3,
           selectInput("go.ont", "Ontology", choices = c("BP","MF", "CC"), selected = "BP", width = "250px")),
    column(width = 3,
           sliderInput("go.maxgeneset", "Limit the geneset size", value = 200,
                       min = 100, max = 3000, step = 100, width = "250px"))
  ),
  
  
  fluidRow(
    div(
      style = "text-align:center; margin-top:30px; margin-bottom:30px",
      actionButton("run.go", "Run GO", width = "400px", 
                   style="background-color:darkgreen; color:white; font-weight:600; font-size:120%")
    )
  ),
  
  
  
  navbarPage(
    title = tags$span(
      style = "color:gold ; font-weight:600 ; font-size:120%",
      "GO results"
    ),
    collapsible = T,
    
    #### p1: Datatable     ##########################################
    tabPanel(
      tags$span(
        style = "color:white ; font-weight:600 ; font-size:120%",
        "Data table"
      ),
      
      fluidRow(
        dataTableOutput("tab.GO")
      )
    ),
    
    
    #### p2: Dotplot     ##########################################
    tabPanel(
      tags$span(
        style = "color:white ; font-weight:600 ; font-size:120%",
        "DotPlot"
      ),
      
      fluidRow(
        column(width = 3,
               div(
                 style = "margin-bottom:40px",
                 selectInput("dotplot.arrange.sel", 
                             "Arrange by:",
                             choices = c("Significance", "RichFactor", "BgGenes")),
                 selectInput("dotplot.show.cat", "Select terms", choices = c(), 
                             size = 20, multiple = T, selectize = F, width = "400px"),
                 verbatimTextOutput("dotp_selected_terms", placeholder = T)
               ),
               actionButton("go.dotplot", "Plot Dotplot", class = "btn-sm", width = "90%",
                            style = "font-size:18px; background-color:midnightblue; font-weight:600; margin-bottom:40px;
                            border-radius:10px; margin-top:20px; border-color:cadetblue")),
        column(width = 9,
               div(
                 style = "overflow-x:auto; width:100%; margin-top:10px",
                 plotOutput("go.dotplot", height = "auto")
               ))
      )
    ),
    
    
    
    #### p3: Cnet Plot        --------------------------------------------------
    tabPanel(
      tags$span(
        style = "color:white ; font-weight:600 ; font-size:120%",
        "CnetPlot"
      ),
      p("Now, to concider the potentially biological complexities in which a gene may belong 
        to multiple annotation categories and provide information of numeric changes if 
        available, the cnetplot() function will be used to extract complex associations. 
        It depicts the linkages of genes and biological concepts as a network.",
        style = "font-weight:600 ; color:darkgreen ; margin-top:20px ; margin-bottom:30px"),
      
      fluidRow(
        column(width = 4,
               div(
                 style = "margin-bottom:20px",
                 selectInput("cnetplot.arrange.sel", 
                             "Arrange by:",
                             choices = c("Significance", "RichFactor", "BgGenes")),
                 selectInput("cnetplot.show.cat", "Select terms", choices = c(), 
                             size = 15, multiple = T, selectize = F, width = "400px")
               )),
        column(width = 4,
               div(
                 style = "margin-bottom:20px",
                 p("Selected terms", style = "margin-bottom:5px ; font-weight:600"),
                 verbatimTextOutput("cnetp_selected_terms", placeholder = T)
               )),
        column(width = 4,
               div(
                 style = "margin-bottom:40px",
                 selectInput("cnet.layout", "Cnet layout",
                             choices = c("circle","kk","grid","fr"),
                             selected = "kk"),
                 numericInput("cnet.genesize", "Gene labels size", value = 0.8, step = 0.1, width = "150px"),
                 numericInput("cnet.catsize", "Term labels size", value = 1.4, step = 0.1, width = "150px")
               ),
               actionButton("go.cnetplot", "Plot Correlation network", class = "btn-sm", width = "90%",
                            style = "font-size:18px; background-color:midnightblue; font-weight:600; margin-bottom:40px;
                            border-radius:10px; margin-top:20px; border-color:cadetblue")),
        column(width = 10,
               offset = 1,
               div(
                 style = "overflow-x:auto; margin-top:20px; width:100%",
                 plotOutput("plt.cnetplot", height = "auto")
               ))
      )
    ),
    
    
    
    #### p4: Emap Plot        --------------------------------------------------
    tabPanel(
      
      tags$span(
        style = "color:white ; font-weight:600 ; font-size:120%",
        "EmapPlot"
      ),
      p("Enrichment map organizes enriched terms into a network with edges connecting 
        overlapping gene sets. In this way, mutually overlapping gene sets are tend to 
        cluster together, making it easy to identify functional module.",
        style = "font-weight:600 ; color:darkgreen ; margin-top:20px ; margin-bottom:30px"),
      
      fluidRow(
        column(width = 4,
               div(
                 style = "margin-bottom:20px",
                 selectInput("emapplot.arrange.sel", 
                             "Arrange by:",
                             choices = c("Significance", "RichFactor", "BgGenes")),
                 selectInput("emapplot.show.cat", "Select terms", choices = c(), 
                             size = 30, multiple = T, selectize = F, width = "400px")
               )),
        column(width = 4,
               div(
                 style = "margin-bottom:20px",
                 p("Selected terms", style = "margin-bottom:5px ; font-weight:600"),
                 verbatimTextOutput("emapp_selected_terms", placeholder = T)
               )),
        column(width = 4,
               div(
                 style = "margin-bottom:40px",
                 selectInput("emap.layout", "Select layout", 
                             choices = c("circle","kk","grid","fr"),
                             selected = "kk"),
                 numericInput("emap.cex.node", "Node size", value = 1.1, step = 0.1),
                 numericInput("emap.cex.label", "Label size", value = 1, step = 0.1),
                 numericInput("emap.cex.line", "edge width", value = .4, step = 0.1),
                 checkboxGroupInput("emap.displayparams", "Parameters to display",
                                    choices = c("Show edges", "Group into clusters", "Show cluster legend")),
                 numericInput("emap.edge.min", "Edge min.similarity", value = 0.4, min = 0, max = 1, step = 0.05),
                 numericInput("emap.clust.n", "Nb of clusters", value = 3, step = 1),
                 numericInput("emap.clust.labs", "Cluster n words", value = 3, step = 1)),
               actionButton("go.emapplot", "Plot Enrichment map", class = "btn-sm", width = "90%",
                            style = "font-size:18px; background-color:midnightblue; font-weight:600; margin-bottom:40px;
                            border-radius:10px; margin-top:20px; border-color:cadetblue"),
        ),
        column(width = 10,
               offset = 1,
               div(
                 style = "overflow-x:auto; margin-top:20px; width:100%",
                 plotOutput("plt.emapplot", height = "auto")
               ))
      )
      
  ),
  
  
  ### 2- GSEA      ###################
  h3("2- GSE analysis",
     style = "color:darkred ; font-weight:600 ; background-color:gold ; margin-top:80px ; margin-bottom:20px"),
  
  fluidRow(
    column(width = 3,
           div(
             style = "margin-bottom:10px",
             verbatimTextOutput("provided_ranking")
           ),
           selectInput("gsea_ont", "Ontology", choices = c("BP"), 
                       selected = "BP", width = "250px"),
           selectInput("gsea_keytype", "Keytype", choices = c("ENSEMBL"), 
                       selected = "ENSEMBL", width = "250px"),
           sliderInput("gsea_GsSize", "GeneSetSize", min = 10, max = 2000, step = 10,
                       value = c(10,500), width = "250px"),
           numericInput("gsea_npermsimple", "nPermSimple", value = 10000, width = "250px"),
           actionButton("run.gsea", "Run GSEA", style = "margin-top:20px", width = "250px",
                        style = "font-size:18px; background-color:midnightblue; font-weight:600; margin-bottom:40px;
                            border-radius:10px; margin-top:20px; border-color:cadetblue")
           ),
    column(width = 9,
           h4("- Datatable results:", style = "color:purple; font-weight:600"),
           div(
             style = "margin-top:20px; margin-bottom:30px",
             dataTableOutput("gseaTab")
           ),
           
           navbarPage(
             title = tags$span(
               style = "color:gold ; font-weight:600 ; font-size:120%",
               "GSEA results"
             ), 
             collapsible = T,
             
             #### p1: dotplot   ------------------------------------------------
             tabPanel(
               tags$span(
                 style = "color:white ; font-weight:600 ; font-size:120%",
                 "Dotplot"
               ),
               div(
                 style = "margin-bottom:30px",
                 numericInput("gsea_showcat", "ShowCategory", value = 20, step = 1),
                 radioButtons("gse.arrange", "Arrange by", 
                              choices = c("Significance","Ascending NES", "Descending NES"), inline = T, selected = ""),
                 plotOutput("gseaPlot1", height = "auto")
               )
             ),
             
             #### p2: gsea plot     --------------------------------------------
             tabPanel(
               tags$span(
                 style = "color:white ; font-weight:600 ; font-size:120%",
                 "Gseaplot"
               ),
               div(
                 style = "margin-bottom:20px",
                 column(width = 7,
                        selectInput("terms.gsea.sel", "Select terms", 
                                    choices = c(), 
                                    multiple = T, 
                                    selectize = F,
                                    size = 20,
                                    width = "500px")),
                 column(width = 5,
                        actionLink("conv_terms", "Convert the terms into GO IDs"),
                        verbatimTextOutput("GO_IDs"),
                        actionButton("gseaplot2", "Plot GSEAplot2", style = "margin-top:30px; background-color:navy"))
               ),
               plotOutput("gseaPlot2", height = "auto")
             )
           ))
  )
))



# Shiny server fun      #################
# ===================================== #
server <- function(session, input, output){
  
  # I- Import data      ------------------------------------
  Alldata <- eventReactive(input$Load, {
    req(input$Load)
    if(str_detect(input$Load$name, ".csv")){
      data <- read.csv(input$Load$datapath)
    } else if(str_detect(input$Load$name, ".xlsx")){
      data <- read_xlsx(input$Load$datapath)
    }
  })
  
  output$tab_data <- DT::renderDT({
    req(Alldata())
    Alldata <- Alldata()
    datatable(Alldata, options = list(scrollX = T))
  })
  
  
  # II- Build marker list     ----------------------------------
  MarkerList <- reactiveVal(list())
  AllGenes <- reactiveVal(character(0))
  GeneRanks <- reactiveVal(numeric(0))
  AllDEGs <- reactiveVal(character(0))
  AllFCs <- reactiveVal(numeric(0))
  UpGenes <- reactiveVal(character(0))
  DownGenes <- reactiveVal(character(0))
  
  ## 1- Filtering   -----------------------------------------------
  observeEvent(Alldata(), {
    req(Alldata())
    Alldata <- Alldata()
    message("==========================================")
    message(sprintf("%-25s: %s", "Total Nb of genes", nrow(Alldata)))
    
    gene_types <- unique(Alldata$Gene_type)
    updateSelectInput(session, "gene_type_filt", choices = gene_types)
    
    AllGenes(Alldata$ENSEMBL)
    GeneRanks(setNames(Alldata$Rank, Alldata$ENSEMBL))
  })
  
  Filtered_data <- eventReactive(input$ApplyFilter, {
    req(Alldata())
    Alldata <- Alldata()
    
    selected_fc <- input$logfc_threshold
    selected_p.adj <- input$padj_threshold
    selected_genetype <- ifelse(length(input$gene_type_filt)!=0,
                                input$gene_type_filt,
                                NULL)
    
    Apply_filt(obj = Alldata, 
               fc = selected_fc, 
               p.adj = selected_p.adj,
               gene_type = selected_genetype)
  })
  
  ## 2- Overview   -----------------------------------------------
  
  observeEvent(Filtered_data(), {
    req(Filtered_data(), AllGenes(), GeneRanks())
    Filtered_data <- Filtered_data()
    AllGenes <- AllGenes()
    GeneRanks <- GeneRanks()
    
    AllDEGs <- Filtered_data %>% dplyr::pull(Gene_name)
    AllFCs <- setNames(Filtered_data[["log2FoldChange"]], Filtered_data[["Gene_name"]])
    UpGenes <- Filtered_data %>% 
      dplyr::filter(log2FoldChange > 0) %>% 
      dplyr::arrange(desc(log2FoldChange)) %>% 
      dplyr::pull(Gene_name)
    DownGenes <- Filtered_data %>% 
      dplyr::filter(log2FoldChange < 0) %>% 
      dplyr::arrange(log2FoldChange) %>% 
      dplyr::pull(Gene_name)
    
    message("==========================================")
    message("=============  Apply Filtering ===========")
    message("==========================================")
    message(sprintf("%-25s: %s", "Total Nb of DEGs", nrow(Filtered_data)))
    message(sprintf("%-25s: %s", "Nb of Up genes", length(UpGenes)))
    message(sprintf("%-25s: %s", "Nb of Down genes", length(DownGenes)))
    message("==========================================")
    
    AllDEGs(AllDEGs)
    UpGenes(UpGenes)
    DownGenes(DownGenes)
    AllFCs(AllFCs)
    
    MarkerList <- list(Filtered_data = Filtered_data,
                       AllGenes = AllGenes,
                       GeneRanks = GeneRanks,
                       AllFCs = AllFCs(),
                       AllDEGs = AllDEGs(),
                       UpGenes = UpGenes(),
                       DownGenes = DownGenes())
    
    MarkerList(MarkerList)
  })
  
  output$list_summary <- renderPrint({
    req(MarkerList())
    MarkerList <- MarkerList()
    str(MarkerList)
  })
  
  output$tab_data_filt <- DT::renderDT({
    req(Filtered_data())
    Filtered_data <- Filtered_data()
    datatable(Filtered_data, options = list(scrollX = T))
  })
  
  
  
  # III- PEA     ----------------------------------
  
  ## 1- GO analysis     -----------------------------
  GOobject <- eventReactive(input$run.go, {
    req(MarkerList())
    List_markers <- MarkerList()
    
    MyGenes <- switch(input$go.list.genes,
                      "All DEGs" = List_markers[["AllDEGs"]],
                      "Up Genes" = List_markers[["UpGenes"]],
                      "Down Genes" = List_markers[["DownGenes"]])
    
    enrichGO(gene = MyGenes,
             OrgDb = org.Hs.eg.db,
             keyType = input$go.keytype,
             ont = input$go.ont,
             qvalueCutoff = 0.1,
             pvalueCutoff = 0.05,
             readable = T,
             maxGSSize = input$go.maxgeneset,
             minGSSize = 10)
  })
  
  observeEvent(input$run.go, {
    req(input$run.go)
    showNotification("GO analysis is running", duration = 10, type = "message")
  })
  
  observeEvent(GOobject(), {
    req(GOobject())
    GOobject <- GOobject()
    showNotification("GO analysis completed !", duration = 5, type = "message")
    message("> GO analysis completed successfully !!!")
    message("==========================================")
    message(sprintf("%-25s: %s", "Total Nb of terms found", nrow(GOobject@result)))
    
  })
  
  
  ### p1: Result table     ---------------------------------------------------------
  go.tab <- eventReactive(GOobject(), {
    req(GOobject())
    GOobject <- GOobject()
    GOobject@result %>% 
      dplyr::filter(p.adjust < 0.05) %>% 
      dplyr::mutate(RichFactor = round(Count / as.numeric(sub("/\\d+","",BgRatio)),5),
                    BgGenes = as.numeric(sub("/\\d+","",BgRatio))) %>% 
      dplyr::select(Description, RichFactor, p.adjust, GeneRatio, BgRatio, Count, BgGenes, geneID)
  })
  
  observeEvent(go.tab(), {
    req(go.tab())
    go.tab <- go.tab()
    message(sprintf("%-25s: %s", "Nb of significant term", nrow(go.tab)))
  })
  
  output$tab.GO <- DT::renderDT({
    req(go.tab())
    go.tab() %>% 
      dplyr::select(Description, GeneRatio, BgGenes, RichFactor, p.adjust) %>% 
      datatable(options = list(pageLength = 5, scrollX = TRUE))
  })
  
  ### p2: Updating selectinputs :   ----------------------------------------
  observeEvent(go.tab(), {
    req(go.tab())
    go.tab <- go.tab() 
    updateSelectInput(session, "dotplot.show.cat", choices = go.tab[["Description"]])
    updateSelectInput(session, "cnetplot.show.cat", choices = go.tab[["Description"]])
    updateSelectInput(session, "emapplot.show.cat", choices = go.tab[["Description"]])
  })
  
  
  ### p3: Dot plot      -----------------------------------------------------------
  
  observeEvent(input$dotplot.arrange.sel, {
    req(go.tab())
    go.tab <- go.tab()
    
    # Arrange by:
    go.tab <- switch(input$dotplot.arrange.sel,
                     "Significance" = arrange(go.tab, p.adjust),
                     "RichFactor" = arrange(go.tab, desc(RichFactor)),
                     "BgGenes" = arrange(go.tab, desc(BgGenes)))
    
    # Updating list order :
    updateSelectInput(session, "dotplot.show.cat", choices = go.tab[["Description"]])
  })
  
  # Selected terms:
  output$dotp_selected_terms <- renderPrint({
    selected_terms <- input$dotplot.show.cat
    cat(c(sprintf("%-15s: %s", "Nb of terms", nrow(go.tab())),
          sprintf("%-15s: %s", "Selected terms", length(selected_terms))), 
        sep = "\n")
  })
  
  # Plot Dot plot:
  go.dplot <- eventReactive(input$go.dotplot, {
    req(go.tab())
    go.tab() %>% 
      dplyr::filter(Description %in% input$dotplot.show.cat) %>% 
      dplyr::mutate(Description = ifelse(nchar(Description) <= 50,
                                         Description,
                                         paste0(substr(Description,1,46), "...."))) %>% 
      ggplot(aes(x= RichFactor, y= fct_reorder(Description, RichFactor)))+
      geom_segment(aes(xend= 0, yend= Description))+
      geom_point(aes(color= p.adjust, size= Count))+
      scale_color_viridis_c(guide = guide_colorbar(reverse = T))+
      scale_size_continuous(range = c(3,10))+
      theme_linedraw()+
      theme(panel.grid = element_blank(),
            panel.border = element_blank(),
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5, colour = "darkred", 
                                      margin = margin(b=0.2, unit = "in")),
            plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"),
            axis.title.x = element_text(size = 16, face = "bold", colour = "darkred", 
                                        margin = margin(t=0.2, unit = "in")),
            axis.title.y = element_blank(),
            axis.text = element_text(size = 14, face = "bold"),
            legend.title = element_text(size = 15, face = "bold", colour = "darkred",
                                        margin = margin(b=0.2, unit = "in")),
            legend.text = element_text(size = 13),
            legend.box.margin = margin(l=0.2, unit = "in"))+
      labs(title = "GO Enriched Terms")
  })
  
  output$go.dotplot <- renderPlot({
    req(go.dplot())
    go.dplot()
  }, width = 800, height = 700)
  
  
  ### p4: Cnet plot      -----------------------------------------------------------
  cnetplot_sel_terms <- reactiveVal(character(0))
  
  observeEvent(input$cnetplot.arrange.sel, {
    req(go.tab())
    go.tab <- go.tab()
    
    # Arrange by:
    go.tab <- switch(input$cnetplot.arrange.sel,
                     "Significance" = arrange(go.tab, p.adjust),
                     "RichFactor" = arrange(go.tab, desc(RichFactor)),
                     "BgGenes" = arrange(go.tab, desc(BgGenes)))
    
    # Updating list order :
    updateSelectInput(session, "cnetplot.show.cat", choices = go.tab[["Description"]])
  })
  
  observeEvent(input$cnetplot.show.cat, {
    selected_terms <- input$cnetplot.show.cat
    cnetplot_sel_terms(selected_terms)
  })
  
  # Selected terms:
  output$cnetp_selected_terms <- renderPrint({
    selected_terms <- input$cnetplot.show.cat
    cat(c(sprintf("%-15s: %s", "Nb of terms", nrow(go.tab())),
          sprintf("%-15s: %s", "Selected terms", length(selected_terms))), 
        sep = "\n")
  })
  
  # Construct the cnet plot:
  CNetPlot <- eventReactive(input$go.cnetplot, {
    req(GOobject(),  MarkerList(), Filtered_data())
    GOobject <- GOobject()
    Filtered_data <- Filtered_data()
    selected_terms <- cnetplot_sel_terms()
    List_markers <- MarkerList()
    fc <- List_markers[["AllFCs"]]
    
    gene_lbl_size <- input$cnet.genesize
    cat_lbl_size <- input$cnet.catsize
    selected_layout <- input$cnet.layout
    
    col_gradient <- switch(input$go.list.genes,
                           "All DEGs" = scale_color_gradientn(
                             colours = c("midnightblue", "white", "darkred"),
                             values = c(0, 0.5, 1),
                             limits = c(-4, 4)
                           ),
                           "Up Genes" = scale_color_gradientn(
                             colours = c("#ffffff", "#ffcc22", "#991111", "#500000"),
                             values = c(0, 0.3, 0.7, 1),
                             limits = c(0, 4)
                           ),
                           "Down Genes" = scale_color_gradientn(
                             colours = c("#000050", "#111199", "#22ccff", "#ffffff"),
                             values = c(0, 0.3, 0.7, 1),
                             limits = c(-4, 0)
                           ))
    
    cnetplot(x = GOobject,
             showCategory = selected_terms,
             layout = selected_layout,
             cex.params = list(gene_node = 0.7, 
                               gene_label = gene_lbl_size, 
                               category_node = 1.5,
                               category_label = cat_lbl_size),
             color.params = list(category = "#2277cc",
                                 gene = "#552299",
                                 edge = T,
                                 foldChange = fc))+
      labs(color = "logFC")+
      theme(legend.box.margin = margin(l=0.3, unit = "in"),
            legend.title = element_text(size = 14, face = "bold", colour = "darkred", 
                                        margin = margin(t=0.2,b=0.5, unit = "in")),
            legend.text = element_text(size = 10, face = "bold"))+
      col_gradient
  })
  
  # Plot cnet:
  output$plt.cnetplot <- renderPlot({
    req(CNetPlot())
    CNetPlot()
  }, height = 1000, width = 1000)
  
  
  
  ### p4: Emap plot      -----------------------------------------------------------
  emapplot_sel_terms <- reactiveVal(character(0))
  
  observeEvent(input$emapplot.arrange.sel, {
    req(go.tab())
    go.tab <- go.tab()
    
    # Arrange by:
    go.tab <- switch(input$emapplot.arrange.sel,
                     "Significance" = arrange(go.tab, p.adjust),
                     "RichFactor" = arrange(go.tab, desc(RichFactor)),
                     "BgGenes" = arrange(go.tab, desc(BgGenes)))
    
    # Updating list order :
    updateSelectInput(session, "emapplot.show.cat", choices = go.tab[["Description"]])
  })
  
  observeEvent(input$emapplot.show.cat, {
    selected_terms <- input$emapplot.show.cat
    emapplot_sel_terms(selected_terms)
  })
  
  # Selected terms:
  output$emapp_selected_terms <- renderPrint({
    selected_terms <- input$emapplot.show.cat
    cat(c(sprintf("%-15s: %s", "Nb of terms", nrow(go.tab())),
          sprintf("%-15s: %s", "Selected terms", length(selected_terms))), 
        sep = "\n")
  })
  
  # Apply pairwise term-sim to GOobj:
  GOobject_paired <- eventReactive(GOobject(), {
    req(GOobject())
    GOobject <- GOobject()
    pairwise_termsim(GOobject, method = "JC", showCategory = 250)
  })
  
  
  # Construct the emapplot:
  emapplot <- eventReactive(input$go.emapplot, {
    req(GOobject_paired(), MarkerList())
    GOobject_paired <- GOobject_paired()
    selected_terms <- emapplot_sel_terms()
    display_params <- input$emap.displayparams
    show_edges <- "Show edges" %in% display_params
    group_clusters <- "Group into clusters" %in% display_params
    show_legend <- "Show cluster legend" %in% display_params
    
    enrichplot::emapplot(GOobject_paired,
                         showCategory = selected_terms,
                         repel = T,
                         edge.params = list(show = show_edges, min = input$emap.edge.min),
                         cex.params = list(category_node = input$emap.cex.node, 
                                           category_label = input$emap.cex.label, 
                                           line = input$emap.cex.line,
                                           label_group = 1),
                         cluster.params = list(cluster = group_clusters, 
                                               method = stats::kmeans, 
                                               n = input$emap.clust.n, 
                                               legend = show_legend, 
                                               label_style = "shadowtext", 
                                               label_words_n = input$emap.clust.labs, 
                                               label_format = 30),
                         layout.params = list(layout = input$emap.layout),
                         node_label = "category")+
      theme(legend.box.margin = margin(l=0.2, unit = "in"),
            plot.margin = margin(t= 0.1, b=0.1, unit = "in"),
            legend.title = element_text(face = "bold", 
                                        colour = "darkred",
                                        size = 14,
                                        margin = margin(b=0.3, unit = "in")))+
      scale_fill_gradientn(colours = c("darkred","gold"),
                           values = c(0,1),
                           limits = c(0, 0.05))+
      labs(fill= "p.adjust")
  })
  
  # Plot emapplot:
  output$plt.emapplot <- renderPlot({
    req(emapplot())
    emapplot()
  }, height = 900, width = 900)
  
  
  
  
  
  
  ## 2- GSEA     -----------------------------
  
  output$provided_ranking <- renderPrint({
    req(Filtered_data(),MarkerList())
    MarkerList <- MarkerList()
    ranked_genes <- MarkerList[["GeneRanks"]]
    cat(paste0("Provided ranked genes : ",length(ranked_genes)),
        sep = "\n")
  })
  
  GSEAobj <- eventReactive(input$run.gsea, {
    req(MarkerList())
    List_markers <- MarkerList()
    ranked_genes <- List_markers[["GeneRanks"]]
    
    gseGO(geneList  = ranked_genes,
          OrgDb        = org.Hs.eg.db,
          ont          = "BP",
          keyType      = "ENSEMBL", 
          minGSSize    = input$gsea_GsSize[1],
          maxGSSize    = input$gsea_GsSize[2],
          pvalueCutoff = 0.05,
          verbose      = FALSE,
          nPermSimple = 10000) %>% 
      setReadable(OrgDb = 'org.Hs.eg.db', keyType = "ENSEMBL")
  })
  
  
  ### p1: Table        --------------------------------------------------------
  # GSEA table:
  gseaTab <- eventReactive(GSEAobj(), {
    req(GSEAobj())
    GSEAobj <- GSEAobj()
    GSEAobj@result %>% 
      dplyr::select(Description, setSize, NES, p.adjust, rank)
  })
  
  output$gseaTab <- DT::renderDT({
    req(gseaTab())
    gseaTab <- gseaTab()
    gseaTab %>% 
      datatable(options = list(pageLength = 5, scrollX = TRUE))
  })
  
  
  ### p2: Dotplot        --------------------------------------------------------
  
  # prepare dotplot:
  gseaDplot <- eventReactive(input$gse.arrange, {
    req(input$gse.arrange, GSEAobj())
    GSEAobj <- GSEAobj()
    
    GSEAobj <- switch(input$gse.arrange,
                      "Significance" = GSEAobj %>% dplyr::arrange(p.adjust),
                      "Ascending NES" = GSEAobj %>% dplyr::arrange(NES),
                      "Descending NES" = GSEAobj %>% dplyr::arrange(desc(NES)))
    
    dotplot(GSEAobj, 
            showCategory=input$gsea_showcat, 
            color = "p.adjust", 
            x="NES")+
      scale_size_continuous(range = c(3,8))+
      theme(axis.title.x = element_text(size = 14, face = "bold", colour = "darkred", margin = margin(t=10)),
            axis.text = element_text(face = "bold", size = 13),
            legend.title = element_text(size = 14, face = "bold", colour = "darkred", margin = margin(b=10)),
            legend.text = element_text(size = 12))
  })
  
  # GSEA plot1:
  output$gseaPlot1 <- renderPlot({
    req(gseaDplot())
    gseaDplot()
  }, width = 800, height = 900)
  
  
  ### p3: Gseaplot        --------------------------------------------------------
  # update sel terms:
  observeEvent(gseaTab(), {
    req(gseaTab())
    gseaTab <- gseaTab()
    updateSelectInput(session, "terms.gsea.sel", choices = gseaTab[["Description"]])
  })
  
  Selected.GO.IDs <- reactiveVal(character(0))
  
  observeEvent(input$conv_terms, {
    req(GSEAobj())
    GSEAobj <- GSEAobj()
    res <- GSEAobj@result
    res <- res %>% dplyr::filter(Description %in% input$terms.gsea.sel)
    selected_IDs <- rownames(res)
    
    Selected.GO.IDs(selected_IDs)
  })
  
  output$GO_IDs <- renderPrint({
    req(Selected.GO.IDs())
    Selected.GO.IDs <- Selected.GO.IDs()
    cat(Selected.GO.IDs, sep = "\n")
  })
  
  gseaPlot2 <- eventReactive(input$gseaplot2, {
    req(GSEAobj(), Selected.GO.IDs())
    GSEAobj <- GSEAobj()
    l <- list()
    for(i in rownames(GSEAobj@result)){
      l[[i]] <- GSEAobj@result[i,"Description"]
    }
    
    selected_IDs <- Selected.GO.IDs()
    gseaplot2(GSEAobj, 
              title = ifelse(length(selected_IDs) == 1,
                             l[[selected_IDs]],
                             ""),
              geneSetID = selected_IDs, 
              color = MyPalette[seq_along(selected_IDs)], 
              pvalue_table = F, 
              subplots = 1:2, 
              base_size = 14, 
              rel_heights = c(3,.3))
  })
  
  
  # GSEA plot2:
  output$gseaPlot2 <- renderPlot({
    req(gseaPlot2())
    gseaPlot2()
  }, width = 800, height = 600)
  
  
  
  
  
  
}





# Shiny build     ########################
# ====================================== #
shinyApp(ui, server)