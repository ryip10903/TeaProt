library(shiny)
library(shinydashboard)
library(fgsea)
library(data.table) 
library(ggplot2)
library(shinyWidgets)
library(markdown)
library(clusterProfiler)
library(dplyr)
library(DT)
library(shinycssloaders)
library(msigdb)
library(GSEABase)
library(vissE)
library(igraph)
library(ggpubr)
library(patchwork)
library(networkD3)
#library(pubgr)
#library(ggpubr)
#library(org.Hs.eg.db)
#library(org.Mm.eg.db)



# Define UI for app ----
ui <- dashboardPage(skin = 'black', title = "TeaProt", 
                    
  # Dashboardheader              
  dashboardHeader(title = tags$a(href='https://github.com/ryip10903/Protein_annotation_app', target = '_blank', tags$img(src=paste0("teaprot.svg"), height = "70%", width = "auto", align = "middle"))),
    
  # Sidebar panel to navigate tabs ----
  dashboardSidebar(
    sidebarMenu(id = "tabs",
                menuItem("START", tabName = "start", icon = icon("bookmark")),
                menuItem("urPTMdb", tabName = "urptmdb", icon = icon("database")),
                menuItem("Analysis", tabName = "ANALYSIS", icon = icon("bookmark"),
                         menuSubItem("View Data", tabName = "view", icon = icon("bookmark")),
                         menuSubItem("p-value & fold change", tabName = "PVALUE", icon= icon("info-circle")),
                         menuSubItem("annotations", tabName = "dugi", icon= icon("chart-bar")),
                         menuSubItem("enrichment", tabName = "genep",icon= icon("ellipsis-v")),
                         menuSubItem("fgsea", tabName = "fgseaa"))
                )
    ),
  
  # Main panel for displaying outputs ----
  dashboardBody(
    
    # The global options for the UI and meta tags are defined here, this includes the custom.css stylesheet into the body.
    tags$head(tags$meta(name = "description", content = "TeaProt is an easy to use and interactive tool to analyze proteomics data."),
              tags$meta(name = "keywords", content = "protein, proteomics, transcriptomics, analysis, visualization")),

    
  #UI layout for each tab
    tabItems(
      
      tabItem(tabName = "start",
              div(class = "jumbotron", style="background-image: url(dna-banner.svg); background-color:white; background-size: cover;", 
                  HTML("<center><h1>Welcome to TeaProt!</h1></center>"), 
                  HTML("<center><p>The online proteomics/transcriptomics analysis pipeline featuring novel underrepresented PTM genesets.</p></center>")),
              
              fluidRow(
                
                # README box
                box(status = 'primary', includeMarkdown("README.md")),
                
                # Demo data box
                box(title = "Demo data", status = 'primary', 
                    
                    selectizeInput("demoset", label = NULL, choices = list("" ,"Mouse tongue data" = "Parker"), selected = "Parker", options = list(placeholder = 'Select dataset')), 
                    htmlOutput('demoref'), 
                    uiOutput("dl_demoset_ui")), 
                
                # Inputs box
                box(title = "Inputs", status = 'primary', 
                    
                    # Input: Select a file ----
                    fileInput("file", "Choose CSV File", multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", '.xls', '.xlsx')),
                    
                    # Input: Select ID ----
                    selectizeInput(inputId = "col_id", label = "Choose ID column", choices = c(NULL), multiple = TRUE, selected = NULL, options = list(placeholder = "This is a placeholder", maxItems = 1)), 
                    
                    # Input: Select pval ----
                    selectizeInput(inputId = "col_pval", label = "Choose p-value column", choices = c(NULL), multiple = TRUE, selected = NULL, options = list(placeholder = "This is a placeholder", maxItems = 1)),  
                    
                    # Input: Select FC ----
                    selectizeInput(inputId = "col_fc", label = "Choose fold-change column", choices = c(NULL), multiple = TRUE, selected = NULL, options = list(placeholder = "This is a placeholder", maxItems = 1)), 
                    
                    # Input: Select a species ----
                    selectInput("species", "Choose Species", c("human", "mouse")),
                    
                    # Input: Select significance cutoff ----
                    sliderTextInput(inputId = "param_pval", label = "Choose a p-value cutoff:", choices = c(1, 0.1, 0.05, 0.01, 0.001), selected = 0.05, grid = TRUE),
                    
                    # Input: Select significance cutoff ----
                    sliderTextInput(inputId = "param_fc", label = "Choose a (log2) fold-change cutoff:", choices = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), selected = 0, grid = TRUE),
                    
                    # Display selected input parameters
                    verbatimTextOutput ("param_text"),
                    
                    # Start analysis button ----
                    actionButton("button", 'Start!', style='font-weight:600')
                    
                    ))),
      
      tabItem(tabName = "urptmdb", 
              
              div(class = "jumbotron", style="background-image: url(ptm-banner.svg); background-color:white; background-size: cover;", 
                  HTML("<center><h1 style='color:White;'>urPTMdb</h1></center>"), 
                  HTML("<center><p style='color:White;'>The underrepresented PTM gene-set database.</p></center>")),
              
              fluidRow(box(status = 'primary', includeMarkdown("README_urptmdb.md")),
                
                column(width = 6,
                box(title = "Download urPTMdb", status = 'primary', width = 12,

                    HTML('
    <table cellspacing=5>
    <tr><td style="padding-right: 10px">Number of studies:</td><td>58</td></tr>
    <tr><td style="padding-right: 10px">Number of PTMs:</td><td>18</td></tr>
    <tr><td style="padding-right: 10px">Number of gene-sets:</td><td>141</td></tr>
    <tr><td style="padding-right: 10px">Filesize:</td><td>1,188 KB</td></tr></table><br>'),
                    
                    HTML("<p>urPTMdb is generated by analyzing the genes reported by many studies to create novel PTM-related gene-sets. urPTMdb is provided in three formats, containing either the original identifier, or formats where genes from other species have been converted to the species of interest. It is recommended to download the database for the species you plan to analyze. In TeaProt, the database use is determined by the species selected at the start of the analysis.</p><br>
                         <ul>
                         <li>urPTMdb Original - Contains the gene identifiers as reported in the original studies</li>
                         <li>urPTMdb Human - All mouse genes have been converted to human homologs using homologene</li>
                         <li>urPTMdb Mouse - All human genes have been converted to mouse homologs using homologene</li>
                         </ul>"),
                    
                    downloadButton("downloadurptmdb", label = "Download urPTMdb - Original"),
                    downloadButton("downloadurptmdb_hs", label = "Download urPTMdb - Human"),
                    downloadButton("downloadurptmdb_mm", label = "Download urPTMdb - Mouse")), 
                
                box(title = "Browse geneset", status = 'primary', width = 12, 
                    selectizeInput(inputId = "geneset_selected", label = "Select a geneset", choices = c(read.delim(file = "database/urPTMdb/urptmdb_latest.gmt")[,1]), multiple = FALSE, selected = NULL),
                    uiOutput("geneset_info")),
                
                box(title = "Geneset jaccard network", status = 'primary', width = 12, 
                    uiOutput("network_d3"), HTML("Gene-set Jaccard index network of urPTMdb. The Jaccard index indicates the similarity between two gene-sets, where the connected nodes have an index > 0.15.")) 
                
                )  ) ),
      
      # View annotated table tab
      tabItem(tabName = "view",        
              fluidRow(box(title = 'About the Table',solidHeader = TRUE, status = 'primary',
              HTML("<p align='justify'> User-uploaded input data is annotated with information from various sources.
                   The annotated table contain information of: </p>
                   <ol>
                      <li> Drug Interaction </li>
                      <li> Cell ontology </li>
                      <li> Associated disease </li>
                  </ol>
                   <p> Export options are available at the bottom of the table </p>"))),
              fluidRow(
                box(title = 'Annotated table', DT::DTOutput("contents") %>% withSpinner(), width = 12) 
                )),
      
      tabItem(tabName = "ANALYSIS"),
      
      # P-value and fold-change distribution tab
      tabItem(tabName = "PVALUE",
        fluidRow( 
          box(title = 'About the Analysis',solidHeader = TRUE, status = 'primary',
                   HTML("<p align='justify'> Analysis are performed to analyze the p-values and fold-changes
                        of your data.</p>
                        <ol>
                            <li> Bar graphs that show the distribution of p-values and fold-changes in the data</li>
                            <li> Volcano plot that shows the fold-changes and corresponding p-values of each data point</li>
                            <ul> 
                              <li>(Hover onto each data point to view the exact values) </li>
                            </ul>
                        </ol>  "))),
        fluidRow(
          box( title = 'Distributions', plotOutput("histogram_pvalue") %>% withSpinner(), uiOutput("hist_download_ui")),
          box( title = 'Volcano plot', plotly::plotlyOutput("volcano_interactive") %>% withSpinner(), uiOutput("volcano_download_ui")  ))),
      
      # Drug interaction
      tabItem(tabName = "dugi", 
        fluidRow ( 
          box(title = 'About the Analysis',solidHeader = TRUE, status = 'primary',
              HTML("<p align='justify'> Your data is mapped with online databases to provide annotations.
              For each sets of graphs below, your data is mapped to a different database.
              The first graph of each section displays the number of genes that could be annotated by the mapped database.
              The second graph displays the annotation</p>
                   "))),
        fluidRow( box(title = 'Drug-gene interaction', width = 12, plotOutput("bargraph_drug") %>% withSpinner())),
        fluidRow( box(title = 'Subcellular localization', width = 12, plotOutput("bargraph_loc") %>% withSpinner())),
        fluidRow( box(title = 'IMPC procedure', width = 12, plotOutput("bargraph_impc") %>% withSpinner())),
        fluidRow( box(title = 'DisGeNet disease', width = 12, plotOutput("bargraph_disgenet") %>% withSpinner())),
        fluidRow( box(title = 'BRENDA enzymatic reactions', width = 12, plotOutput("bargraph_brenda") %>% withSpinner()))),
      
      # Gene pathway 
      tabItem(tabName = "genep",
        fluidRow(box(title = 'About the Analysis',solidHeader = TRUE, status = 'primary', 
        HTML("<p align='justify'> Analyses are performed to demonstrate the changes in gene expressions in relation to several annotations including (1) subcellular localization, (2) DisGeNet, (3) Drug-gene interactions and (4) International Mouse Phenotyping Consortium interactions. 
             A Pearson’s Chi-squared test based on protein annotations (subcellular localization) indicates whether specific annotations are primarily found in upregulated, downregulated or non-significant (NS) proteins. Only localizations with positive residuals in the upregulated group are shown. The data in the figure is colored by Pearson residuals, and sized by the absolute Pearson residuals.</p>"))),
        
        fluidRow( box(title = "Subcellular Localizations", plotOutput("contingency_loc") %>% withSpinner(), uiOutput("contingency_download_loc_ui"), width = 12)),
        fluidRow( box(title = "DisGeNet", plotOutput("contingency_disgenet") %>% withSpinner(), uiOutput("contingency_download_disgenet_ui"), width = 12)),
        fluidRow( box(title = "Drug-gene interactions", plotOutput("contingency_drug") %>% withSpinner(), uiOutput("contingency_download_drug_ui"), width = 12)),
        fluidRow( box(title = "IMPC genotype-phenotype Associations", plotOutput("contingency_impc") %>% withSpinner(), uiOutput("contingency_download_impc_ui"), width = 12))
        ),
              
      # Fgsea
      tabItem(tabName = "fgseaa", 
              
              # Row containing explanation and inputs
              fluidRow(box(title = 'About the Analysis', solidHeader = TRUE, status = 'primary', 
                           HTML("<p align='justify'>
                           This analysis is dependent on the fold-change values in your data. 
                           The graph displays the most enriched biological pathways that are associated with the differential expressions. 
                           In the input section, choose the geneset collection that you want the analysis to be based on. After running the analysis, the results will be displayed in the following tabs:
                                </p>
                                <ul>
                                <li><b>panel: </b>Image showing the top x positively and negatively enriched pathways</li>
                                <li><b>table: </b>Table showing all fgsea results</li>
                                <li><b>volcano: </b>Volcano plot showing the p-value and NES of each tested geneset</li>
                                <li><b>single: </b>Tab showing fgsea enrichment and coloured volcano plot for a single geneset of interest</li>
                                </ul>")),
                       box(title ="input", solidHeader = TRUE, status = 'primary',
                           selectInput('fgseadb', 'Choose gene-sets', choices = list(MSigDB = c(`h: hallmark gene sets` = 'h', `c1: positional gene sets` = 'c1', `c2: curated gene sets` = 'c2', `c3: regulatory target gene sets` = 'c3', `c4: computational gene sets` = 'c4',`c5: ontology gene sets` = 'c5',`c6: oncogenic signature gene sets` = 'c6',`c7: immunologic signature gene sets` = 'c7',`c8: cell type signature gene sets` = 'c8'), Transcription = c(`CHEA3 - ENCODE` = 'encode', `CHEA3 - REMAP` = 'remap', `CHEA3 - Literature` = 'literature'), urPTMdb = c(`underrepresented PTMs` = 'urptmdb')), selectize = FALSE),
                           actionButton('buttonfg', 'Start Analysis'))),
              
              # Row containing results
              fluidRow(tabBox(title = "FGSEA", id = "tabset1", width = 12, 
                        tabPanel("panel", height = "70vh", selectInput("fgnumber", "Choose Number of Pathways to display", c(10, 20, 30)), plotOutput("fgseaplot") %>% withSpinner()), 
                        tabPanel("table", height = "70vh", DT::dataTableOutput("fgseatable")), 
                        tabPanel("volcano", height = "70vh", plotly::plotlyOutput("fgseaplot_volcano") %>% withSpinner()),
                        tabPanel("single", height = "70vh", uiOutput("fgsea_select_ui"), fluidRow(column(5,plotOutput("fgseaplot_single") %>% withSpinner()), column(7,plotly::plotlyOutput("volcano_single") %>% withSpinner())), uiOutput("volcano_single_download_ui"))
      )))
      
    )

  )

)

# Define server logic to read selected file ----
server <- function(input, output, session) {
 
  # Create a reactiveValues object called mydata
  mydata <- reactiveValues()
  
  # Print the selected tab for testing purposes
  observe({print(input$tabs)})
  
  # Print the uploaded file for testing purposes
  observe({print(input$file)})
  
  # Generate the interactive network or urPTMdb genesets
  output$network_d3 <- renderUI({
    
    gph_d3 <- readRDS(file = "data/urptmdb_network.rds")

    networkD3::forceNetwork(Links = gph_d3$links, Nodes = gph_d3$nodes, Source = 'source', Target = 'target', NodeID = 'name', Group = 'group', zoom = TRUE, charge = -10, opacity = 0.9)
    
    })
  
  # Run tests once file is uploaded, save column names as input options, and reset all parameters and plots back to NULL
  observeEvent(input$file, {
    
    # Add code to set all reactiveValues back to NULL
    mydata$protdf <- NULL
    mydata$array <- NULL
    mydata$gene_array <- NULL
    mydata$visse1 <- NULL
    mydata$fgsea_res <- NULL
    mydata$fgsea_array <- NULL
    mydata$fgsea_pathways <- NULL
    mydata$fgsea_volcano <- NULL
    mydata$fgsea_descriptiontable <- NULL
    mydata$chisq_loc <- NULL
    mydata$chisq_impc <- NULL
    mydata$chisq_drug <- NULL
    mydata$chisq_disgenet <- NULL
    
    
    if(sub("^.*\\.","", input$file$datapath) %in% c("csv", "txt", "xls", "xlsx", "XLSX")){
    
    df <- load_data(input$file$datapath)
      
    col_names <- colnames(as.data.frame(df))
    
    if(length(col_names) > 2){
      
      updateSelectizeInput(session,"col_id", choices = col_names, selected = col_names[1])
      updateSelectizeInput(session,"col_pval", choices = col_names, selected = col_names[2])
      updateSelectizeInput(session,"col_fc", choices = col_names, selected = col_names[3])
      
    } else {
      
      updateSelectizeInput(session,"col_id", choices = col_names)
      updateSelectizeInput(session,"col_pval", choices = col_names)
      updateSelectizeInput(session,"col_fc", choices = col_names)
      
    }
    
    
    }
    
    
  })

  # Report user parameter selection
  output$param_text <- renderText({
    paste("Uploaded file: ", input$file[1], "\n", "ID column: ", input$col_id, "\n", "p-value column: ", input$col_pval, "\n", "Fold-change column: ", input$col_fc, "\n", "Selected species: ", input$species, "\n", "Selected p-value cutoff: ", input$param_pval, "\n", "Selected fold-change cutoff: ", input$param_fc, "\n" , sep="")
  })
  
  # Process data once the start button is pressed
  observeEvent(input$button, {
    
    req(input$file)

    sendSweetAlert(session = session, title = "Notification", 
                   text = "Data Mapping in progress", type = "warning",
                   btn_labels = NA,
                   closeOnClickOutside = FALSE, showCloseButton = FALSE)
    
    
    df <- load_data(input$file$datapath)
    
  #df is user's uploaded value
    df <- as.data.frame(df) %>% dplyr::select(input$col_id, input$col_pval, input$col_fc)

    if(validate_input(df) != TRUE){
      
      if(validate_size(df) != TRUE){
        
        sendSweetAlert(session = session, title = "data size error", text = "Please confirm that the 3 required columns (ID, pval, FC) have been selected.", type = "error", btn_labels = NA, closeOnClickOutside = TRUE)
        
      } else{
        
        sendSweetAlert(session = session, title = "Error", text = "Input data not valid", type = "error", btn_labels = NA, closeOnClickOutside = TRUE)
        
        print(validate_id(df[,1]))
        print(validate_pval(df[,2]))
        print(validate_fc(df[,3]))
        
      }
      
      

      
    }
    
    validate(need(validate_size(df) == TRUE, label = "error", message = "input data wrong size"))
    validate(need(validate_input(df) == TRUE, label = "error", message = "input data not validated"))
        

    
  #converting the uploaded column names to match our commands
    names(df)[1] <- "ID"
    names(df)[2] <- 'pvalue'
    names(df)[3] <- 'fold_change'

    df <- df %>% mutate(ID = sub("\\;.*", "", .$ID))
    df <- df %>% mutate(ID = sub("\\:.*", "", .$ID))
    
    
  #pre-made function
  #converting gene names to entrez id
    df <- cp_idconvert(df, cp_idtype(df$ID))
    
    if(validate_size(df) != TRUE){
      
      sendSweetAlert(session = session, title = "data size error", text = "Input data validated, but no ID's were converted. Please confirm that a valid ID type is present.", type = "error", btn_labels = NA, closeOnClickOutside = TRUE)
      
    }
    
    validate(need(validate_size(df) == TRUE, label = "error", message = "input data wrong size"))
    
    
    
  # joining the user's data with mouse gene id
    df <- left_join(df, entrezmapping(input$species))
    
  #joining user's data with drug interaction data
   
    df <- df %>% mutate(ID = tolower(ID)) 
    
    #Loading drug interaction data
    db_dgidb <- readRDS(file = "database/DGIdb_genename_drugname.Rds") %>% dplyr::select(drug_name, gene_name) %>% mutate(gene_name = tolower(gene_name)) %>% `colnames<-`(c("drug_name", "gene_name"))
    db_dgidb <- db_dgidb %>% group_by(gene_name) %>% summarise(drug_name = paste(drug_name, collapse = "|"))
    
    df <- left_join(df, db_dgidb, by = c("ID" = "gene_name"))
  
    
  
  # joining user's data with gene location data
    # Load Human Protein Atlas, but edit to include ancestor localization terms for more accurate overlap
    db_hpa <- readRDS(file = "database/db_hpa.Rds")
    ancestors <- readxl::read_excel(path = 'database/localization_ancestors.xlsx' ) %>% filter(!is.na(Ancestors)) %>% mutate(term_ancestor = paste0(Term, ";", Ancestors))
    
    db_hpa$CP_loc <- db_hpa$`IF main protein location`
    
    for(i in 1:nrow(ancestors)){
      db_hpa$CP_loc <- db_hpa$CP_loc %>% sub(ancestors$Term[i], ancestors$term_ancestor[i], .)
    }
    
    db_hpa <- db_hpa %>% mutate(CP_loc = sapply(strsplit(db_hpa$CP_loc, ";"), function(x) paste(unique(x), collapse = ";"))) %>% dplyr::select(-ENSG, -Uniprot, -`IF main protein location`, -`HyperLOPIT location`, -Reliability)
    
    df <- left_join (df, db_hpa, by = c('ID' = 'Gene'))
    
    
    # Annotate IMCP data
    db_impc_procedure <- read.csv("database/IMPC/impc_procedure.csv") %>% mutate(marker_symbol = tolower(marker_symbol)) %>% mutate(gene_id = as.character(gene_id)) 
    db_impc_parameter <- read.csv("database/IMPC/impc_parameter.csv") %>% mutate(marker_symbol = tolower(marker_symbol)) %>% mutate(gene_id = as.character(gene_id)) 
    
    if(input$species != "mouse"){
      df <- left_join(df, db_impc_procedure %>% dplyr::select(1,2) %>% distinct(), by = c('ID' = "marker_symbol")) %>% left_join(., db_impc_parameter %>% dplyr::select(1,2) %>% distinct(), by = c('ID' = "marker_symbol"))
    } else {
      df <- left_join(df, db_impc_procedure %>% dplyr::select(3,2) %>% distinct(), by = c('EntrezGeneID' = "gene_id")) %>% left_join(., db_impc_parameter %>% dplyr::select(3,2) %>% distinct(), by = c('EntrezGeneID' = "gene_id"))
    }
    
    # Annotate DisGeNet data
    db_disgenet <- read.csv("database/DisGeNet/disgenet_genesymbol.csv") %>% mutate(geneSymbol = tolower(geneSymbol)) %>% `colnames<-`(c("geneSymbol", "DisGeNet_disease"))
    df <- left_join(df, db_disgenet, by = c('ID' = "geneSymbol"))
    
    # Annotate Brenda data
    db_brenda <- read.csv("database/BRENDA/brenda_processed.csv") %>% dplyr::select(entrez, reaction) %>% mutate(entrez = as.character(entrez))
    df <- left_join(df, db_brenda, by = c('EntrezGeneID' = "entrez"))
    
    #Making columns used for contigency table
    df <- df %>% mutate(direction = case_when( (pvalue < input$param_pval & fold_change > input$param_fc) ~ "up", (pvalue < input$param_pval & fold_change < -(input$param_fc)) ~ "down", TRUE ~ "NS"))
    
  # Converting EntrezGeneID column into an array
  # So it can later be used for enrichment analysis
    array <- df %>% dplyr::select(3) %>% unlist()
    names(array) <- df %>% dplyr::select(EntrezGeneID) %>% unlist()
    array <- sort(array)
  
    gene_array <- df %>% dplyr::select(3) %>% unlist()
    names(gene_array) <- df %>% dplyr::select(ID) %>% unlist() %>% as.character() %>% toupper()
    gene_array <- sort(gene_array)
    
    
    mydata$protdf <- df
    mydata$array <- array 
    mydata$gene_array <-gene_array
    
    sendSweetAlert(session = session, title = "Notification", 
                   text = "Mapping has completed! You can start your Analysis! ", type = "success",
                   closeOnClickOutside = TRUE, showCloseButton = FALSE)
    

  })
  
  ## Demodata
  output$demoref <- renderUI({
    if(input$demoset == "Parker"){
      HTML('
    <table cellspacing=5>
    <tr><td style="padding-right: 10px">Publication:</td><td><a href="https://doi.org/10.3390/proteomes9020022" target="_blank">Western Diet Induced Remodelling of the Tongue Proteome</a></td></tr>
    <tr><td style="padding-right: 10px">Number of proteins:</td><td>5,296</td></tr>
    <tr><td style="padding-right: 10px">Filesize:</td><td>258 KB</td></tr></table><br>')
    } })
  
  # Download button
  # Downloadbutton demoset ----
  output$dl_demoset_ui <- renderUI({
    if(input$demoset != "")
      downloadButton("downloadData", label = "Download Demo Dataset")
  })
  
  
  
  
  # Download demo dataset handler ----
  output$downloadData <- downloadHandler(
    filename <- "tongue_de_mouse.xlsx",
    
    content =function(file) {message("Action: User downloaded the demo dataset")
      file.copy(paste0("data/tongue_de_mouse.xlsx"), file)},
    contentType = NA,
  )
  
  # Download demo dataset handler ----
  output$downloadurptmdb <- downloadHandler(
    filename <- "urptmdb_latest.gmt",
    
    content =function(file) {message("Action: User downloaded urPTMdb")
      file.copy(paste0("database/urPTMdb/urptmdb_latest.gmt"), file)},
    contentType = NA,
  )
  
  # Download demo dataset handler ----
  output$downloadurptmdb_mm <- downloadHandler(
    filename <- "urptmdb_latest_mm.gmt",
    
    content =function(file) {message("Action: User downloaded urPTMdb")
      file.copy(paste0("database/urPTMdb/urptmdb_latest_mm.gmt"), file)},
    contentType = NA,
  )
  
  # Download demo dataset handler ----
  output$downloadurptmdb_hs <- downloadHandler(
    filename <- "urptmdb_latest_hs.gmt",
    
    content =function(file) {message("Action: User downloaded urPTMdb")
      file.copy(paste0("database/urPTMdb/urptmdb_latest_hs.gmt"), file)},
    contentType = NA,
  )
  
  
  
  
  output$geneset_info <- renderUI({
    
    req(input$geneset_selected)
    print(input$geneset_selected)
    
    df <- read.delim(file = "database/urPTMdb/urptmdb_latest.gmt")
    description <- df %>% filter(V1 == input$geneset_selected) %>% .[,2] %>% unlist %>% unname
    genes <- df %>% filter(V1 == input$geneset_selected) %>% .[,3:ncol(df)] %>% unlist %>% unname
    genes <- genes[!is.na(genes)] %>% sort()
    
    if(stringr::str_extract(input$geneset_selected, "^....") == "GOBP"){
      
      HTML(paste('
    <table cellspacing=5>
    <tr><td style="padding-right: 10px">Source:</td><td>', stringr::str_extract(input$geneset_selected, "^....") ,'</td></tr>
    <tr><td style="padding-right: 10px">PTM:</td><td>', sub("^.....", "",input$geneset_selected )  ,'</td></tr>
    <tr><td style="padding-right: 10px">Description:</td><td>', description,'</td></tr></table><br>',
    '<br><h4>Included genes:</h4>',
    '<p align="justify"><code>', paste(genes, collapse = ", "),'</code></p>', sep = ""))
      
    } else {
      
      HTML(paste('
    <table cellspacing=5>
    <tr><td style="padding-right: 10px">Source:</td><td>', strsplit(input$geneset_selected, "_")[[1]][1],'</td></tr>
    <tr><td style="padding-right: 10px">Species:</td><td>', strsplit(input$geneset_selected, "_")[[1]][2] ,'</td></tr>
    <tr><td style="padding-right: 10px">PTM:</td><td>', strsplit(input$geneset_selected, "_")[[1]][3]  ,'</td></tr>
    <tr><td style="padding-right: 10px">Description:</td><td>', description,'</td></tr></table><br>',
    '<br><h4>Included genes:</h4>',
    '<p align="justify"><code>', paste(genes, collapse = ", "),'</code></p>', sep = ""))
      
    }
    
  })
  
  
    
  ## Fgsea enrichment analysis ##---------------------------------------
  observeEvent(input$buttonfg, {
    
    if(!exists('mydata$array')){sendSweetAlert(session = session, title = "Error", text = "Please upload your data first", type = "error")}
    
    
    req(mydata$array, input$fgnumber, input$fgseadb, input$species)
    
    sendSweetAlert(session = session, title = "Notification",
                   btn_labels = NA,
                   text = "Analysis in Progress", type = "warning",
                   closeOnClickOutside = FALSE , showCloseButton = FALSE)
    
    # Select gene-sets based on input$fgseadb ----
    if(input$fgseadb %in% c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8")){
      
      if(input$species == "mouse"){msigdb_mm =  msigdb::getMsigdb(org = 'mm', id = 'EZID', version = '7.2')} else {msigdb_mm =  msigdb::getMsigdb(org = 'hs', id = 'EZID', version = '7.2')}
      
      msigdb_mm = appendKEGG(msigdb_mm)
      msigdb_mm = subsetCollection(msigdb_mm, input$fgseadb)
      genedb <- geneIds(msigdb_mm)
      msig_db_gsea <- data.frame(gs_name = rep(names(genedb), lengths(genedb)), gene_symbol = unlist(genedb, use.names=TRUE))
      
      descriptiontable <- list()
      
      for(i in 1:length(msigdb_mm)){
      
        descriptiontable[i] <- description(msigdb_mm[[i]])
          
      }
      
      descriptiontable <- data.frame(pathway = names(msigdb_mm), description = unlist(descriptiontable))
      
      pathways = genedb
      
      }
    
    if(input$fgseadb %in% c("encode", "remap", "literature", "urptmdb")){
      
      if(input$fgseadb == "encode"){pathway <- read.delim(file = "database/CHEA3/ENCODE_ChIP-seq.gmt", header = FALSE); descriptiontable <- data.frame(pathway = read.delim(file = "database/CHEA3/ENCODE_ChIP-seq.gmt", header = FALSE)[,1], description = read.delim(file = "database/CHEA3/ENCODE_ChIP-seq.gmt", header = FALSE)[,1])}
      if(input$fgseadb == "remap"){pathway <- read.delim(file = "database/CHEA3/ReMap_ChIP-seq.gmt", header = FALSE); descriptiontable <- data.frame(pathway = read.delim(file = "database/CHEA3/ReMap_ChIP-seq.gmt", header = FALSE)[,1], description = read.delim(file = "database/CHEA3/ReMap_ChIP-seq.gmt", header = FALSE)[,1])}
      if(input$fgseadb == "literature"){pathway <- read.delim(file = "database/CHEA3/Literature_ChIP-seq.gmt", header = FALSE); descriptiontable <- data.frame(pathway = read.delim(file = "database/CHEA3/Literature_ChIP-seq.gmt", header = FALSE)[,1], description = read.delim(file = "database/CHEA3/Literature_ChIP-seq.gmt", header = FALSE)[,1])}
      if(input$fgseadb == "urptmdb" & input$species == "mouse"){pathway <- read.delim(file = "database/urPTMdb/urptmdb_latest_mm.gmt", header = TRUE) %>% dplyr::select(-V2); descriptiontable <- read.delim(file = "database/urPTMdb/urptmdb_latest_mm.gmt", header = TRUE) %>% dplyr::select(V1, V2) %>% `colnames<-`(c("pathway", "description"))}
      if(input$fgseadb == "urptmdb" & input$species != "mouse"){pathway <- read.delim(file = "database/urPTMdb/urptmdb_latest_hs.gmt", header = TRUE) %>% dplyr::select(-V2); descriptiontable <- read.delim(file = "database/urPTMdb/urptmdb_latest_hs.gmt", header = TRUE) %>% dplyr::select(V1, V2) %>% `colnames<-`(c("pathway", "description"))}
      
      pathway <- pathway %>% tidyr::unite(., col = "genes", -V1, na.rm = TRUE)
      pathway$genes <- sub("\\_\\_.*","", pathway$genes)
      pathway <- pathway %>% filter(genes != "")
      
      
      pathwaylist <- list()
      
      for(i in 1:nrow(pathway)){
        
        pathwaylist[pathway$V1[i]] <- strsplit(toupper(pathway$genes[i]), split = "_")
        
      }
      
      pathways = pathwaylist
      
    }
    
    if(input$fgseadb %in% c("encode", "remap", "literature", "urptmdb")){array <- mydata$gene_array} else {array <- mydata$array}
    
    fgseaRes <- fgsea(pathways = pathways,
          stats = array,
          minSize = 5,
          maxSize = 500)
    
    mydata$fgsea_res <- fgseaRes
    mydata$fgsea_array <- array
    mydata$fgsea_pathways <- pathways
    mydata$fgsea_descriptiontable <- descriptiontable
    
    mydata$fgsea_volcano <- plotly::ggplotly(ggplot(fgseaRes,aes(y=-log10(pval), x=NES, label = pathway)) + geom_point(col = "blue", alpha = 0.2) + theme_bw())
    
    sendSweetAlert(session = session, title = "Notification", 
                   text = "Analysis has completed!", type = "success",
                   closeOnClickOutside = TRUE, showCloseButton = TRUE)
    
    
  })
  
  output$fgseaplot_volcano <- plotly::renderPlotly({ 
    req(mydata$fgsea_volcano)
    return(mydata$fgsea_volcano)
  })
  
  output$fgseaplot <- renderPlot({
    
    req(mydata$fgsea_res, mydata$fgsea_pathways, mydata$fgsea_array, input$fgnumber)
    
    array <- isolate(mydata$fgsea_array)
    pathways <- isolate(mydata$fgsea_pathways)
    fgseaRes <- isolate(mydata$fgsea_res)
    
    topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=as.numeric(isolate(input$fgnumber))), pathway]
    topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=as.numeric(isolate(input$fgnumber))), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    
    p1 <- plotGseaTable(pathways[topPathways], array, fgseaRes, gseaParam=0.5, render = FALSE)
    plot(ggpubr::as_ggplot(p1))
    
  })
  
  output$fgsea_select_ui <- renderUI({
    
    req(mydata$fgsea_res)
    
    pathways <- isolate(mydata$fgsea_res)[,1]
    
    selectInput("fgsea_select", "Select gene-set to display", choices = pathways, width = '100%')
    
  })
  
  output$fgseaplot_single <- renderPlot({
    
    req(mydata$fgsea_pathways, mydata$fgsea_array, input$fgsea_select)
    
    array <- isolate(mydata$fgsea_array)
    pathways <- isolate(mydata$fgsea_pathways)

    
    p1 <- plotEnrichment(pathways[[input$fgsea_select]], array) + labs(title=input$fgsea_select)
    plot(p1)
    
  })
  
  
  output$volcano_single <- plotly::renderPlotly({ 
    req(mydata$protdf, input$fgsea_select, mydata$fgsea_pathways, input$fgseadb)
    
    pathways <- isolate(mydata$fgsea_pathways)
    
    if(input$fgseadb %in% c("encode", "remap", "literature", "urptmdb")){
      
      df <- isolate(mydata$protdf) %>% mutate(target = toupper(ID) %in% toupper(pathways[[input$fgsea_select]])) %>% arrange(target)
      
    } else {
      
      df <- isolate(mydata$protdf) %>% mutate(target = toupper(EntrezGeneID) %in% toupper(pathways[[input$fgsea_select]])) %>% arrange(target)
      
    }
    
    p1 <- plotly::ggplotly(ggplot(df,aes(y=-log10(pvalue), x=fold_change, label = ID, col = target)) + geom_point() + theme_bw() + scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "#AAAAEE15")))
    
    return(p1)
  })
  
  
  output$volcano_single_download_ui <- renderUI({
    req(mydata$protdf, input$fgsea_select, mydata$fgsea_pathways)
    tagList(
      downloadButton("dl_enrichment_single_image", label = "Download enrichment Image"),
      downloadButton("dl_volcano_single_image", label = "Download volcano Image")
      
    )
  })

  
  output$dl_enrichment_single_image <- downloadHandler(
    filename = function() {
      paste("enrichment_fgsea", Sys.Date(), ".svg", sep="")
    },
    content = function(file) {
      
      array <- isolate(mydata$fgsea_array)
      pathways <- isolate(mydata$fgsea_pathways)
      
      
      p <- plotEnrichment(pathways[[input$fgsea_select]], array) + labs(title=input$fgsea_select)

      svglite::svglite(filename = file, width = 6, height = 3)
      plot(p)
      dev.off()
    }
  )
  
    
  output$dl_volcano_single_image <- downloadHandler(
    filename = function() {
      paste("volcano_fgsea", Sys.Date(), ".svg", sep="")
    },
    content = function(file) {
      
      pathways <- isolate(mydata$fgsea_pathways)
      
      if(input$fgseadb %in% c("encode", "remap", "literature", "urptmdb")){
        
        df <- isolate(mydata$protdf) %>% mutate(target = toupper(ID) %in% toupper(pathways[[input$fgsea_select]])) %>% arrange(target)
        
      } else {
        
        df <- isolate(mydata$protdf) %>% mutate(target = toupper(EntrezGeneID) %in% toupper(pathways[[input$fgsea_select]])) %>% arrange(target)
        
      }
      
      p <- ggplot(df,aes(y=-log10(pvalue), x=fold_change, label = ID, col = target)) + geom_point() + theme_bw() + scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "#AAAAEE15")) + theme(legend.position = "null")
      
      svglite::svglite(filename = file, width = 2.5, height = 2.5)
      plot(p)
      dev.off()
    }
  )
  
  
  
  
  output$fgseatable <- DT::renderDataTable(server = FALSE, {
    
    req(mydata$fgsea_res, mydata$fgsea_descriptiontable)
    
    df <- mydata$fgsea_res %>% arrange(-padj) %>% left_join(., mydata$fgsea_descriptiontable) %>% dplyr::select(pathway, description, everything())
    
    df[,3] <- formatC(as.data.frame(df)[,3], digits = 2, format = 'e')
    df[,4] <- formatC(as.data.frame(df)[,4], digits = 2, format = 'e')
    df[,5] <- formatC(as.data.frame(df)[,5], digits = 3) 
    df[,7] <- formatC(as.data.frame(df)[,6], digits = 3) 
    df[,7] <- formatC(as.data.frame(df)[,7], digits = 3) 

    return(DT::datatable(df, extensions = 'Buttons', rownames = FALSE, options = list(dom = 'tpB', fixedColumns = TRUE, autoWidth = FALSE, pagingType = "numbers", scrollX = T, buttons = list(
      list(extend = 'csv', filename = paste("fgsea_", input$fgseadb, "_", Sys.Date(), sep="")),
      list(extend = 'excel', filename = paste("fgsea_", input$fgseadb, "_", Sys.Date(), sep="")) ))) %>%
             DT::formatStyle(columns = colnames(data), fontSize = '80%'))
    
  })
  
  
    
  # Creates an overview of our mapped data
  output$contents <- DT::renderDataTable(server = FALSE, {
    
    req(mydata$protdf)
    
    df <- mydata$protdf
    
    df[,2] <- formatC(as.data.frame(df)[,2], digits = 2, format = 'e')
  
    df[,3] <- formatC(as.data.frame(df)[,3], digits = 3) 
  
    return(DT::datatable(df, extensions = 'Buttons', rownames = FALSE, options = list(dom = 'tpB', fixedColumns = TRUE, autoWidth = FALSE, pagingType = "numbers", scrollX = T, buttons = list(
      list(extend = 'csv', filename = paste("TeaProt_annotated_data_", Sys.Date(), sep="")),
      list(extend = 'excel', filename = paste("TeaProt_annotated_data_", Sys.Date(), sep="")) ) )) %>%
      DT::formatStyle(columns = colnames(data), fontSize = '80%'))
    
  })
  


  

  
  #P value UI-1
  output$volcano_interactive <- plotly::renderPlotly({ 
    req(mydata$protdf, input$param_pval, input$param_fc)
    
    p1 <- plotly::ggplotly(ggplot(isolate(mydata$protdf),aes(y=-log10(pvalue), x=fold_change, label = ID, col = ((pvalue < isolate(input$param_pval)) & (abs(fold_change) > isolate(input$param_fc)) ))) + geom_point(alpha = 0.3) + scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "skyblue")) + theme_bw() + theme(legend.position = "none"))
    
    return(p1)
    })
  
  output$volcano_download_ui <- renderUI({
    req((mydata$protdf))
    downloadButton("dl_volcano_image", label = "Download volcano Image")
  })
  
  output$dl_volcano_image <- downloadHandler(
    filename = function() {
      paste("volcano_", Sys.Date(), ".svg", sep="")
    },
    content = function(file) {
      
      p <- ggplot(isolate(mydata$protdf),aes(y=-log10(pvalue), x=fold_change, label = ID, col = ((pvalue < isolate(input$param_pval)) & (abs(fold_change) > isolate(input$param_fc)) ))) + geom_point(alpha = 0.3) + scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "skyblue")) + theme_bw() + theme(legend.position = "none")
      
      svglite::svglite(filename = file, width = 2.5, height = 2.5)
      plot(p)
      dev.off()
    }
  )
  
  
  
  #P value UI-2
  output$histogram_pvalue <- renderPlot({ 
    req(mydata$protdf)
    
    # A Density plot for P value------------------------------
    
    p1 <- ggplot(isolate(mydata$protdf), aes(x = pvalue)) + geom_histogram(bins = 40, fill = "blue", alpha = 0.2) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme_bw()
    p2 <- ggplot(isolate(mydata$protdf), aes(x = fold_change)) + geom_histogram(bins = 100, fill = "blue", alpha = 0.2) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme_bw()
    
    p <- (p1 + p2)
    
    plot(p)
    })
  
  output$hist_download_ui <- renderUI({
    req((mydata$protdf))
      downloadButton("dl_hist_image", label = "Download histogram Image")
  })
  
  output$dl_hist_image <- downloadHandler(
    filename = function() {
      paste("distribution_histogram_", Sys.Date(), ".svg", sep="")
    },
    content = function(file) {
      
      p1 <- ggplot(isolate(mydata$protdf), aes(x = pvalue)) + geom_histogram(bins = 40, fill = "blue", alpha = 0.2) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme_bw()
      p2 <- ggplot(isolate(mydata$protdf), aes(x = fold_change)) + geom_histogram(bins = 100, fill = "blue", alpha = 0.2) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme_bw()
      
      p <- (p1 + p2)
      
      svglite::svglite(filename = file, width = 5, height = 2.5)
      plot(p)
      dev.off()
    }
  )


  
  #annotation dgi
  output$bargraph_drug <- renderPlot({
    req((mydata$protdf))
    
    df <- isolate(mydata$protdf) %>% tidyr::separate_rows(drug_name, sep = "\\|") %>% filter(drug_name != "") %>% 
      filter(drug_name != "") %>% group_by(drug_name) %>% summarise (frequency =n()) %>% slice_max(frequency, n = 40, with_ties = FALSE)
    
    p1 <- (isolate(mydata$protdf) %>% dplyr::select(ID, drug_name) %>% mutate(drug_name = !is.na(drug_name)) %>% group_by(ID) %>% summarise(n = sum(drug_name)) %>% mutate(n = (n != 0)) %>% ggplot(., aes(x = n)) + geom_bar(stat = "count", fill = "blue", alpha = 0.2, col = "black") + theme_bw())
    p2 <- ggplot(df, aes(x=frequency, y=reorder(drug_name, frequency))) + geom_bar(stat = "identity", fill = "blue", alpha = 0.2, col = "black") + theme_minimal() + labs(x = "Frequency", y = "")
    
    p <- (p1 + p2)
    
    plot(p)
  })
  
  #annotation loc
  output$bargraph_loc <- renderPlot({
    req((mydata$protdf))
    
    df <- isolate(mydata$protdf) %>% tidyr::separate_rows(CP_loc) %>% filter(CP_loc != "") %>% 
      filter(CP_loc != "") %>% group_by(CP_loc) %>% summarise (frequency =n())
    
    p1 <- (isolate(mydata$protdf) %>% dplyr::select(ID, CP_loc) %>% mutate(CP_loc = !is.na(CP_loc)) %>% group_by(ID) %>% summarise(n = sum(CP_loc)) %>% mutate(n = (n != 0)) %>% ggplot(., aes(x = n)) + geom_bar(stat = "count", fill = "blue", alpha = 0.2, col = "black") + theme_bw())
    p2 <- ggplot(df, aes(x=frequency, y=reorder(CP_loc, frequency))) + geom_bar(stat = "identity", fill = "blue", alpha = 0.2, col = "black") + theme_minimal() + labs(x = "Frequency", y = "")
    
    p <- (p1 + p2)
    
    plot(p)
  })
  
  #annotation impc
  output$bargraph_impc <- renderPlot({
    req((mydata$protdf))
    
    df <- isolate(mydata$protdf) %>% tidyr::separate_rows(impc_significant_procedure_name, sep = "\\|") %>% filter(impc_significant_procedure_name != "") %>% 
      filter(impc_significant_procedure_name != "") %>% group_by(impc_significant_procedure_name) %>% summarise (frequency =n()) %>% slice_max(frequency, n = 40, with_ties = FALSE)
    
    p1 <- (isolate(mydata$protdf) %>% dplyr::select(ID, impc_significant_procedure_name) %>% mutate(impc_significant_procedure_name = !is.na(impc_significant_procedure_name)) %>% group_by(ID) %>% summarise(n = sum(impc_significant_procedure_name)) %>% mutate(n = (n != 0)) %>% ggplot(., aes(x = n)) + geom_bar(stat = "count", fill = "blue", alpha = 0.2, col = "black") + theme_bw())
    p2 <- ggplot(df, aes(x=frequency, y=reorder(impc_significant_procedure_name, frequency))) + geom_bar(stat = "identity", fill = "blue", alpha = 0.2, col = "black") + theme_minimal() + labs(x = "Frequency", y = "")
    
    p <- (p1 + p2)
    
    plot(p)
  })
  
  #annotation disgenet
  output$bargraph_disgenet <- renderPlot({
    req((mydata$protdf))
    
    df <- isolate(mydata$protdf) %>% tidyr::separate_rows(DisGeNet_disease, sep = "\\|") %>% filter(DisGeNet_disease != "") %>% 
      filter(DisGeNet_disease != "") %>% group_by(DisGeNet_disease) %>% summarise (frequency =n()) %>% slice_max(frequency, n = 40, with_ties = FALSE)
    
    p1 <- (isolate(mydata$protdf) %>% dplyr::select(ID, DisGeNet_disease) %>% mutate(DisGeNet_disease = !is.na(DisGeNet_disease)) %>% group_by(ID) %>% summarise(n = sum(DisGeNet_disease)) %>% mutate(n = (n != 0)) %>% ggplot(., aes(x = n)) + geom_bar(stat = "count", fill = "blue", alpha = 0.2, col = "black") + theme_bw())
    p2 <- ggplot(df, aes(x=frequency, y=reorder(DisGeNet_disease, frequency))) + geom_bar(stat = "identity", fill = "blue", alpha = 0.2, col = "black") + theme_minimal() + labs(x = "Frequency", y = "")
    
    p <- (p1 + p2)
    
    plot(p)
  })
  
  #annotation disgenet
  output$bargraph_brenda <- renderPlot({
    req((mydata$protdf))
    
    df <- isolate(mydata$protdf) %>% tidyr::separate_rows(reaction, sep = "\\|") %>% filter(reaction != "") %>% 
      filter(reaction != "") %>% group_by(reaction) %>% summarise (frequency =n()) %>% slice_max(frequency, n = 40, with_ties = FALSE)
    
    p1 <- (isolate(mydata$protdf) %>% dplyr::select(ID, reaction) %>% mutate(reaction = !is.na(reaction)) %>% group_by(ID) %>% summarise(n = sum(reaction)) %>% mutate(n = (n != 0)) %>% ggplot(., aes(x = n)) + geom_bar(stat = "count", fill = "blue", alpha = 0.2, col = "black") + theme_bw())
    p2 <- ggplot(df, aes(x=frequency, y=reorder(reaction, frequency))) + geom_bar(stat = "identity", fill = "blue", alpha = 0.2, col = "black") + theme_minimal() + labs(x = "Frequency", y = "")
    
    p <- (p1 + p2)
    
    plot(p)
  })
  
  
  
  


  # Contingency
  output$contingency_loc <- renderPlot({
    req((mydata$protdf))
    
    # Localization analysis
    df <- isolate(mydata$protdf) %>% dplyr::select(ID, CP_loc, direction) %>% distinct()
    
    locs <- df$CP_loc %>% unique() %>% strsplit(., split = ";") %>% unlist() %>% unique()
    locs <- locs[!is.na(locs)]
    
    contingency <- list()
    
    for(i in 1:length(locs)){
      
      df_con <- df %>% mutate(group = grepl(locs[i], .$CP_loc))
      
      contingency[[locs[i]]] <- c(df_con %>% filter(direction == "up" & group == TRUE) %>% nrow(),
                                  df_con %>% filter(direction == "down" & group == TRUE) %>% nrow(),
                                  df_con %>% filter(direction == "NS" & group == TRUE) %>% nrow())
    }
    
    contingency[[length(locs) + 1]] <- c(df_con %>% filter(direction == "up") %>% nrow(),
                                         df_con %>% filter(direction == "down") %>% nrow(),
                                         df_con %>% filter(direction == "NS") %>% nrow())
    
    contingency <- do.call(rbind, contingency) %>% `rownames<-`(c(locs, "total")) %>% `colnames<-`(c("up", "down", "NS"))
    
    contingency <- contingency[rowSums(contingency) != 0,]
    
    chisq <- chisq.test(contingency)
    
    mydata$chisq_loc <- chisq
    
    corrplot::corrplot(mydata$chisq_loc$residuals %>% as.data.frame() %>% mutate(total = grepl("total", rownames(.))) %>% arrange(total, -up) %>% dplyr::select(-total) %>% t, is.cor = FALSE, title = "", mar=c(0,0,1,0), col.lim = c(floor(min(mydata$chisq_loc$residuals)), ceiling(max(mydata$chisq_loc$residuals))))
  })
  
  output$contingency_impc <- renderPlot({
    req((mydata$protdf))
    
    # IMPC analysis ----
    df <- isolate(mydata$protdf) %>% dplyr::select(ID, impc_significant_procedure_name, direction) %>% distinct()
    
    procedure <- df$impc_significant_procedure_name %>% unique() %>% strsplit(., split = "\\|") %>% unlist() %>% unique()
    procedure <- procedure[!is.na(procedure)]
    
    contingency <- list()
    
    for(i in 1:length(procedure)){
      
      df_con <- df %>% mutate(group = grepl(procedure[i], .$impc_significant_procedure_name))
      
      contingency[[procedure[i]]] <- c(df_con %>% filter(direction == "up" & group == TRUE) %>% nrow(),
                                       df_con %>% filter(direction == "down" & group == TRUE) %>% nrow(),
                                       df_con %>% filter(direction == "NS" & group == TRUE) %>% nrow())
    }
    
    contingency[[length(procedure) + 1]] <- c(df_con %>% filter(direction == "up") %>% nrow(),
                                              df_con %>% filter(direction == "down") %>% nrow(),
                                              df_con %>% filter(direction == "NS") %>% nrow())
    
    contingency <- do.call(rbind, contingency) %>% `rownames<-`(c(procedure, "total")) %>% `colnames<-`(c("up", "down", "NS"))
    
    contingency <- contingency[rowSums(contingency) != 0,]
    
    chisq <- chisq.test(contingency)
    
    mydata$chisq_impc <- chisq
    
    corrplot::corrplot(mydata$chisq_impc$residuals %>% as.data.frame() %>% mutate(total = grepl("total", rownames(.))) %>% arrange(total, -up) %>% dplyr::select(-total) %>% t, is.cor = FALSE, title = "", mar=c(0,0,1,0), col.lim = c(floor(min(mydata$chisq_impc$residuals)), ceiling(max(mydata$chisq_impc$residuals))))
  })
  
  
  output$contingency_drug <- renderPlot({
    req((mydata$protdf))
    
    # drug analysis ----
    df <- isolate(mydata$protdf) %>% dplyr::select(ID, drug_name, direction) %>% distinct()
    
    # Limit analysis to top 50 drugs (by number of annotations)
    procedure <- df$drug_name %>% unique %>% strsplit(., split = "\\|") %>% unlist %>% table %>% sort(decreasing = TRUE) %>% head(40) %>% names()
    procedure <- procedure[!is.na(procedure)]
    
    contingency <- list()
    
    for(i in 1:length(procedure)){
      
      df_con <- df %>% mutate(group = grepl(procedure[i], .$drug_name))
      
      contingency[[procedure[i]]] <- c(df_con %>% filter(direction == "up" & group == TRUE) %>% nrow(),
                                       df_con %>% filter(direction == "down" & group == TRUE) %>% nrow(),
                                       df_con %>% filter(direction == "NS" & group == TRUE) %>% nrow())
    }
    
    contingency[[length(procedure) + 1]] <- c(df_con %>% filter(direction == "up") %>% nrow(),
                                              df_con %>% filter(direction == "down") %>% nrow(),
                                              df_con %>% filter(direction == "NS") %>% nrow())
    
    contingency <- do.call(rbind, contingency) %>% `rownames<-`(c(procedure, "total")) %>% `colnames<-`(c("up", "down", "NS"))
    
    contingency <- contingency[rowSums(contingency) != 0,]
    
    chisq <- chisq.test(contingency)
    
    mydata$chisq_drug <- chisq
    
    corrplot::corrplot(mydata$chisq_drug$residuals  %>% as.data.frame() %>% mutate(total = grepl("total", rownames(.))) %>% arrange(total, -up) %>% dplyr::select(-total) %>% t, is.cor = FALSE, title = "", mar=c(0,0,1,0), col.lim = c(floor(min(mydata$chisq_drug$residuals)), ceiling(max(mydata$chisq_drug$residuals))))
  })
  
  
  output$contingency_download_drug_ui <- renderUI({
    req((mydata$chisq_drug))
    tagList(
      downloadButton("dl_contingency_drug", label = "Download contingency - Drug-gene interactions Table"),
      downloadButton("dl_contingency_drug_image", label = "Download contingency - Drug-gene interactions Image")
    )
  })
  
  output$dl_contingency_drug <- downloadHandler(
    filename = function() {
      paste("contingency_drug_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(cbind(mydata$chisq_drug$observed %>% `colnames<-`(paste(colnames(.), "observed", sep = "_")),
                      round(mydata$chisq_drug$expected,1) %>% `colnames<-`(paste(colnames(.), "expected", sep = "_")),
                      round(mydata$chisq_drug$residuals,1) %>% `colnames<-`(paste(colnames(.), "residuals", sep = "_"))), file)
    }
  )
  
  output$dl_contingency_drug_image <- downloadHandler(
    filename = function() {
      paste("contingency_drug_", Sys.Date(), ".svg", sep="")
    },
    content = function(file) {
      svglite::svglite(filename = file, width = 15, height = 8)
      corrplot::corrplot(mydata$chisq_drug$residuals  %>% as.data.frame() %>% mutate(total = grepl("total", rownames(.))) %>% arrange(total, -up) %>% dplyr::select(-total) %>% t, is.cor = FALSE, title = "", mar=c(0,0,1,0), col.lim = c(floor(min(mydata$chisq_drug$residuals)), ceiling(max(mydata$chisq_drug$residuals))))
      dev.off()
    }
  )
  
  
  
  
  
  
  output$contingency_disgenet <- renderPlot({
    req((mydata$protdf))
    
    # disgenet analysis ----
    df <- isolate(mydata$protdf) %>% dplyr::select(ID, DisGeNet_disease, direction) %>% distinct()
    
    # Limit analysis to top 50 drugs (by number of annotations)
    procedure <- df$DisGeNet_disease %>% unique %>% strsplit(., split = "\\|") %>% unlist %>% table %>% sort(decreasing = TRUE) %>% head(40) %>% names()
    procedure <- procedure[!is.na(procedure)]
    
    contingency <- list()
    
    for(i in 1:length(procedure)){
      
      df_con <- df %>% mutate(group = grepl(procedure[i], .$DisGeNet_disease))
      
      contingency[[procedure[i]]] <- c(df_con %>% filter(direction == "up" & group == TRUE) %>% nrow(),
                                       df_con %>% filter(direction == "down" & group == TRUE) %>% nrow(),
                                       df_con %>% filter(direction == "NS" & group == TRUE) %>% nrow())
    }
    
    contingency[[length(procedure) + 1]] <- c(df_con %>% filter(direction == "up") %>% nrow(),
                                              df_con %>% filter(direction == "down") %>% nrow(),
                                              df_con %>% filter(direction == "NS") %>% nrow())
    
    contingency <- do.call(rbind, contingency) %>% `rownames<-`(c(procedure, "total")) %>% `colnames<-`(c("up", "down", "NS"))
    
    contingency <- contingency[rowSums(contingency) != 0,]
    
    chisq <- chisq.test(contingency)
    
    mydata$chisq_disgenet <- chisq
    
    corrplot::corrplot(mydata$chisq_disgenet$residuals  %>% as.data.frame() %>% mutate(total = grepl("total", rownames(.))) %>% arrange(total, -up) %>% dplyr::select(-total) %>% t, is.cor = FALSE, title = "", mar=c(0,0,1,0), col.lim = c(floor(min(mydata$chisq_disgenet$residuals)), ceiling(max(mydata$chisq_disgenet$residuals))))
  })
  
  
  output$contingency_download_disgenet_ui <- renderUI({
    req((mydata$chisq_disgenet))
    tagList(
      downloadButton("dl_contingency_disgenet", label = "Download contingency - DisGeNet Table"),
      downloadButton("dl_contingency_disgenet_image", label = "Download contingency - DisGeNet Image")
    )
  })
  
  output$dl_contingency_disgenet <- downloadHandler(
    filename = function() {
      paste("contingency_disgenet_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(cbind(mydata$chisq_disgenet$observed %>% `colnames<-`(paste(colnames(.), "observed", sep = "_")),
                      round(mydata$chisq_disgenet$expected,1) %>% `colnames<-`(paste(colnames(.), "expected", sep = "_")),
                      round(mydata$chisq_disgenet$residuals,1) %>% `colnames<-`(paste(colnames(.), "residuals", sep = "_"))), file)
    }
  )
  
  output$dl_contingency_disgenet_image <- downloadHandler(
    filename = function() {
      paste("contingency_disgenet_", Sys.Date(), ".svg", sep="")
    },
    content = function(file) {
      svglite::svglite(filename = file, width = 15, height = 8)
      corrplot::corrplot(mydata$chisq_disgenet$residuals %>% as.data.frame() %>% mutate(total = grepl("total", rownames(.))) %>% arrange(total, -up) %>% dplyr::select(-total) %>% t, is.cor = FALSE, title = "", mar=c(0,0,1,0), col.lim = c(floor(min(mydata$chisq_disgenet$residuals)), ceiling(max(mydata$chisq_disgenet$residuals))))
      dev.off()
    }
  )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  output$contingency_download_loc_ui <- renderUI({
    req((mydata$chisq_loc))
    tagList(
      downloadButton("dl_contingency_loc", label = "Download contingency - Localization Table"),
      downloadButton("dl_contingency_loc_image", label = "Download contingency - Localization Image")
    )
  })
  
  output$dl_contingency_loc <- downloadHandler(
    filename = function() {
      paste("contingency_localization_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(cbind(mydata$chisq_loc$observed %>% `colnames<-`(paste(colnames(.), "observed", sep = "_")),
                      round(mydata$chisq_loc$expected,1) %>% `colnames<-`(paste(colnames(.), "expected", sep = "_")),
                      round(mydata$chisq_loc$residuals,1) %>% `colnames<-`(paste(colnames(.), "residuals", sep = "_"))), file)
    }
  )
  
  output$dl_contingency_loc_image <- downloadHandler(
    filename = function() {
      paste("contingency_localization_", Sys.Date(), ".svg", sep="")
    },
    content = function(file) {
      svglite::svglite(filename = file, width = 15, height = 8)
      corrplot::corrplot(mydata$chisq_loc$residuals %>% as.data.frame() %>% mutate(total = grepl("total", rownames(.))) %>% arrange(total, -up) %>% dplyr::select(-total) %>% t, is.cor = FALSE, title = "", mar=c(0,0,1,0), col.lim = c(floor(min(mydata$chisq_loc$residuals)), ceiling(max(mydata$chisq_loc$residuals))))
      dev.off()
    }
  )
  
  
  
  
  output$contingency_download_impc_ui <- renderUI({
    req((mydata$chisq_loc))
    tagList(
      downloadButton("dl_contingency_impc", label = "Download contingency - IMPC Table"),
      downloadButton("dl_contingency_impc_image", label = "Download contingency - IMPC Image")
    )
  })
  
  output$dl_contingency_impc <- downloadHandler(
    filename = function() {
      paste("contingency_impc_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(cbind(mydata$chisq_impc$observed %>% `colnames<-`(paste(colnames(.), "observed", sep = "_")),
                      round(mydata$chisq_impc$expected,1) %>% `colnames<-`(paste(colnames(.), "expected", sep = "_")),
                      round(mydata$chisq_impc$residuals,1) %>% `colnames<-`(paste(colnames(.), "residuals", sep = "_"))), file)
    }
  )
  
  output$dl_contingency_impc_image <- downloadHandler(
    filename = function() {
      paste("contingency_impc_", Sys.Date(), ".svg", sep="")
    },
    content = function(file) {
      svglite::svglite(filename = file, width = 15, height = 8)
      corrplot::corrplot(mydata$chisq_impc$residuals  %>% as.data.frame() %>% mutate(total = grepl("total", rownames(.))) %>% arrange(total, -up) %>% dplyr::select(-total) %>% t, is.cor = FALSE, title = "", mar=c(0,0,1,0), col.lim = c(floor(min(mydata$chisq_impc$residuals)), ceiling(max(mydata$chisq_impc$residuals))))
      dev.off()
    }
  )
  

  

}


# Create Shiny app ----
shinyApp(ui, server)

