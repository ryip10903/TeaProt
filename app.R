library(shiny)
library(shinydashboard)
library(fgsea)
library(data.table) 
library(ggplot2)
library(shinyWidgets)
library(markdown)
library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(DT)
library(shinycssloaders)
library(msigdb)
library(GSEABase)
library(vissE)
library(igraph)
library(ggpubr)
library(patchwork)
#library(pubgr)
#library(org.Hs.eg.db)
#library(org.Mm.eg.db)
#library(org.Dr.eg.db)
#library(org.Rn.eg.db)
#library(org.Dm.eg.db)



# Define UI for data upload app ----
ui <- dashboardPage( skin = 'black',
                     
    # App title and header ----
    title = "TeaProt",
  
    dashboardHeader(title = tags$a(href='https://github.com/ryip10903/Protein_annotation_app', target = '_blank',
                                      tags$img(src=paste0("teaprot.svg"), height = "70%", width = "auto", align = "middle"))),
  
  
  # Sidebar panel for inputs ----
  dashboardSidebar(
    sidebarMenu(id = "tabs",
      # Input: Select a file ----
      fileInput("file", "Choose CSV File",
                multiple = FALSE,
                accept = c("text/csv", 
                           "text/comma-separated-values,text/plain",
                           ".csv", '.xls', '.xlsx')),
      
      
      selectInput("species", "Choose Species",
                  c("human", "mouse", "rat", "drosophila", "zebrafish")),
      
    #Start analysis button
    actionButton("button", 'Start Here!', style='font-weight:600' ),
      
    #sidebar Items
      menuItem("START", tabName = "start", icon = icon("bookmark")),
      
      menuItem("View Data", tabName = "view", icon = icon("bookmark")),
      
      menuItem("Analysis", tabName = "ANALYSIS", icon = icon("bookmark"),
              menuSubItem("p-value & fold change", tabName = "PVALUE", icon= icon("info-circle")),
              menuSubItem("annotations", tabName = "dugi"),
              menuSubItem("Gene Pathway", tabName = "genep"),
              menuSubItem("Fgsea analysis", tabName = "fgseaa"),
              menuSubItem('vissE analysis', tabName = 'vvss')
              )


    ),
    sidebarMenuOutput("menu")
  ),
  
  # Main panel for displaying outputs ----
  dashboardBody(
    
    # The global options for the UI and meta tags are defined here
    # This includes the custom.css stylesheet into the body.
    # The meta tags affect the search enginge results when looking for CoffeeProt
    tags$head(tags$meta(name = "description", content = "TeaProt is an easy to use and interactive tool to analyze proteomics data."),
              tags$meta(name = "keywords", content = "protein, proteomics, transcriptomics, analysis, visualization")),
    
    
  #UI layout for each tab
    tabItems(
      
      tabItem(tabName = "start",
              div(class = "jumbotron", style="background-image: url(dna-banner.svg); 
              background-color:white;
                  background-size: cover;", HTML("<center><h1>Welcome to TeaProt!</h1></center>"),
                  HTML("<center><p>For the annotation of Protein/transcript data.</p></center>")),
              fluidRow(
                box(status = 'primary', includeMarkdown("README.md")))),
      
      tabItem(tabName = "view",        
              fluidRow(box(title = 'About the Table',solidHeader = TRUE, status = 'primary',
              HTML("<p align='justify'> User-uploaded input data is annotated with information from various sources, including DGIdb and HPA.</p>"))),
              fluidRow(
                box(title = 'Annotated table', DT::DTOutput("contents") %>% withSpinner()) 
                )),
      
      tabItem(tabName = "ANALYSIS"),
      
      ## P value
      tabItem(tabName = "PVALUE",
        fluidRow( 
          box(title = 'About the Analysis',solidHeader = TRUE, status = 'primary',
                   HTML("<p align='justify'> Analysis are performed to analyze the p-values and fold-changes
                        of your data.</p>"))),
        fluidRow(
          box( title = 'Distributions', plotOutput("histogram_pvalue") %>% withSpinner()),
          box( title = 'Volcano plot', plotly::plotlyOutput("density_pvalue") %>% withSpinner() 
               ))
        ),
      
      ##Drug interaction
      tabItem(tabName = "dugi", 
        fluidRow ( 
          box(title = 'About the Analysis',solidHeader = TRUE, status = 'primary',
              HTML("<p align='justify'> Analysis are performed to identify the number of genes
                        that have known drug interactions.</p>"))),
        fluidRow(
          box(title = 'Drug-gene interaction - Annotations', plotOutput("bargraph_drug") %>% withSpinner())),
        fluidRow(
          box(title = 'Subcellular localization - Annotations', width = 12, plotOutput("bargraph_loc") %>% withSpinner()))),
      
    #Gene pathway 
      tabItem(
        actionButton('buttonpa', 'Start Analysis'),
        selectInput('pgoinput', 'Choose P_value cutoff',
                    c(0.001, 0.01, 0.05)),
        tabName = "genep",
              h2("Gene Pathway Analysis"), 
              fluidRow(column(12, box(plotOutput('goinput')%>% withSpinner(), width =12)),
              column(12, box(tableOutput('gotable'), width = 12))),
        ),
              
      #Fgsea
      tabItem(
        actionButton('buttonfg', 'Start Analysis'), 
        selectInput("fgnumber", "Choose Number of Pahways to be seen",
                    c(10, 20, 30)),
        tabName = "fgseaa", h2("Fgsea Analysis"), 
        fluidRow(column(12,
          box(width = 12, plotOutput("fgseainput") %>% withSpinner() )))),
      

      
    #Visse analysis
      tabItem(
        actionButton('buttonv', 'Start Analysis'),
        selectInput('DBselect', 'Options', c('c1','c2','c3','c4','c5','c6','c7', 'c8',
                                             'h'), 
                    selected = c('h', 'c2', 'c5'),
                    multiple=TRUE, selectize=TRUE),
        tabName = 'vvss', h2('Visse Analysis'),
          fluidRow(column(12, box(title = 'Visse Plot', plotOutput('visseinput')
                                 %>% withSpinner(), width = 12)))
      )
      
    ),
    
   
    
    
    
    
  ))

# Define server logic to read selected file ----
server <- function(input, output, session) {
 
  
  # Create a reactiveValues object called mydata
  mydata <- reactiveValues()
  forout_reactive <- reactiveValues()
  CLP <- reactiveValues()
  
  # The following code runs when a file is uploaded
  # We want to load data, and annotate it with our databases here
  # We first load the data based on its extension (.csv, .txt, .xlsx)
  # Then we load our databases and join the user data with out provided databases
  # Finally we save the df as a reactive object
  
  observe({
    print(input$tabs)
  })
  
  observeEvent(input$button, {
    
    
    sendSweetAlert(session = session, title = "Notification", 
                   text = "Data Mapping in progress :) ", type = "warning",
                   btn_labels = NA,
                   closeOnClickOutside = FALSE, showCloseButton = FALSE)
    
    
  #Loading drug interaction data
    db_dgidb <<- readRDS(file = "database/DGIdb_genename_drugname.Rds") %>% dplyr::select(drug_name, gene_name, interaction_claim_source, interaction_types) %>% mutate(gene_name = tolower(gene_name)) %>% `colnames<-`(c("drug_name", "gene_name", "dgi_interaction_claim_source", "dgi_interaction_types"))
   
    # Load Human Protein Atlas, but edit to include ancestor localization terms for more accurate overlap
    db_hpa <- readRDS(file = "database/db_hpa.Rds")
    ancestors <- readxl::read_excel(path = 'database/localization_ancestors.xlsx' ) %>% 
      filter(!is.na(Ancestors)) %>% 
      mutate(term_ancestor = paste0(Term, ";", Ancestors))
    
    db_hpa$CP_loc <- db_hpa$`IF main protein location`
    
    for(i in 1:nrow(ancestors)){
      db_hpa$CP_loc <- db_hpa$CP_loc %>% sub(ancestors$Term[i], ancestors$term_ancestor[i], .)
    }
    
    db_hpa <<- db_hpa %>% mutate(CP_loc = sapply(strsplit(db_hpa$CP_loc, ";"), function(x) paste(unique(x), collapse = ";"))) %>% dplyr::select(-ENSG, -Uniprot, -`HyperLOPIT location`, -Reliability) %>% dplyr::rename("HPA_IF_protein_location" = `IF main protein location`)
    
   
    print("file loaded")
    
    if(sub("^.*\\.","", input$file$datapath) == "csv"){
      tryCatch(
        {
          df <- read.csv(input$file$datapath, sep = ",", na.strings = c("NA", "Na", "NaN", "NAN", "na", "nan"))
        },
        error = function(e) {
          stop(safeError(e))
        }
      )
      
    } else if(sub("^.*\\.","", input$file$datapath) == "txt" | sub("^.*\\.","", input$file$datapath) == "tsv") {
      tryCatch(
        {
          df <- read.csv(input$file$datapath, sep = "\t", na.strings = c("NA", "Na", "NaN", "NAN", "na", "nan"))
        },
        error = function(e) {
          stop(safeError(e))
        }
      )
    } else if(sub("^.*\\.","", input$file$datapath) == "xls" | 
              sub("^.*\\.","", input$file$datapath) == "xlsx") {
      tryCatch(
        {
          df <- readxl::read_excel(path = input$file$datapath, guess_max = 21474836, na = c("NA", "Na", "NaN", "NAN", "na", "nan"))
        },
        error = function(e) {
          stop(safeError(e))
        }
      )
    }
  #df is user's uploaded value
    df <- as.data.frame(df)
    x1 <<- df
    
    
  #converting the uploaded column names to match our commands
    names(df)[1] <- "ID"
    names(df)[2] <- 'pvalue'
    names(df)[3] <- 'fold_change'
    x2 <<- df
    
    df <- df %>% mutate(ID = sub("\\;.*", "", .$ID))
    df <- df %>% mutate(ID = sub("\\:.*", "", .$ID))
  #pre-made function
  #converting gene names to entrez id
    df <- cp_idconvert(df, cp_idtype(df$ID))
    
    x3 <<- df
    
  # joining the user's data with mouse gene id
    df <- left_join(df, entrezmapping(input$species))
    
  #joining user's data with drug interaction data
   
    df <- df %>% mutate (ID = tolower(ID)) 
    df <- left_join(df, db_dgidb, by = c("ID" = "gene_name"))
  
    x4 <<- df
    

    
   
    
  
  # joining user's data with gene location data
    df <- left_join (df, db_hpa, by = c('ID' = 'Gene'))
    x5 <<- df
    

    
    
  # Converting EntrezGeneID column into an array
  # So it can later be used for enrichment analysis
    array <- df %>% dplyr::select(3) %>% unlist()
    names(array) <- df %>% dplyr::select(EntrezGeneID) %>% unlist()
    array <- sort(array)
    
    x6 <<- array
    
    go_array <- df %>% dplyr::select(3) %>% unlist()
    names(go_array) <- df %>% dplyr::select(EntrezGeneID) %>% unlist() %>% as.character()
    go_array <- sort(go_array)
    
    x7 <<- go_array
   
#attaching df to "mydata" reactive value
  
    mydata$CP_summary  <- df %>% tidyr:: separate_rows(CP_loc) %>% filter(CP_loc != "") %>% 
      filter(CP_loc != "") %>% group_by(CP_loc) %>% summarise (frequency =n())
      
    
    
    z1 <<- mydata$CP_summary

    mydata$protdf <- df
    mydata$array <- array 
    mydata$go_array <-array
    
    
    sendSweetAlert(session = session, title = "Notification", 
                   text = "Mapping has completed! You can start your Analysis! ", type = "success",
                   closeOnClickOutside = TRUE, showCloseButton = FALSE)
    
  })
  


    
  ## visse analysis --------------------------------------------
  
  observeEvent(input$buttonv, {
    
    if(!exists('mydata$protdf')){sendSweetAlert(session = session, title = "Error", text = "Please upload your data first", type = "error")}
    
    req(mydata$protdf)
    sendSweetAlert(session = session, title = "Notification",
                   btn_labels = NA,
                   text = "Analysis in Progress, this will take approximately 10 minutes....", type = "warning",
                   closeOnClickOutside = FALSE , showCloseButton = FALSE)
    
    
    
      
    # Create input data 
    d <- mydata$protdf[,c(1,3)] %>% filter(!is.nan(ID) & !is.na(fold_change)) %>% dplyr::rename(gene = 1, fc = 2) %>% group_by(gene) %>% summarize(fc = mean(fc))
    
    geneList <- d[,2] %>% unlist
    names(geneList) <- as.character(d[,1] %>% unlist %>% stringr::str_to_title())
    geneList <- base::sort(geneList, decreasing = TRUE)
    
    
    #load the MSigDB from the msigdb package
    msigdb_mm =  msigdb::msigdb.v7.2.mm.SYM()
    #append KEGG gene-sets
    msigdb_mm = appendKEGG(msigdb_mm)
    #dplyr::select h, c2, and c5 collections (recommended)
    msigdb_mm = subsetCollection(msigdb_mm, input$DBselect)
    
    # Create an input gene set for the GSEA function
    genedb <- geneIds(msigdb_mm)
    msig_db_gsea <- data.frame(gs_name = rep(names(genedb), lengths(genedb)), gene_symbol = unlist(genedb, use.names=TRUE))
    
    # Run GSEA
    set.seed(36)
    
    em <- GSEA(geneList, TERM2GENE = msig_db_gsea, pvalueCutoff = 0.05)
    geneset_res = em@result$Description
    
    if(nrow(em@result) < 1) {
      sendSweetAlert(session = session, title = "Notification",
        text = "Error,No enrichment found with GSEA", type = "error",
       closeOnClickOutside = TRUE , showCloseButton = TRUE)}
    validate(need(nrow(em@result)!=0, "Error, dataset contains 0 rows after filtering. "))
    
    z2 <<- em
    
    sendSweetAlert(session = session, title = "Notification",
                   btn_labels = NA,
                   text = "step1/4 has been completed", type = "warning",
                   closeOnClickOutside = FALSE , showCloseButton = FALSE)
    
    #create a GeneSetCollection using the gene-set analysis results
    geneset_gsc = msigdb_mm[geneset_res]
    
    #compute gene-set overlap
    gs_ovlap = computeMsigOverlap(geneset_gsc, thresh = 0.25)
    #create an overlap network
    gs_ovnet = computeMsigNetwork(gs_ovlap, msigdb_mm)
    #plot the network
    set.seed(36) #set seed for reproducible layout
    
    sendSweetAlert(session = session, title = "Notification",
                   btn_labels = NA,
                   text = "step 1/4 has been completed", type = "warning",
                   closeOnClickOutside = FALSE , showCloseButton = FALSE)
    
    
 
    #simulate gene-set statistics
    geneset_stats = em@result$NES
    names(geneset_stats) = geneset_res
    head(geneset_stats)
    
    #plot the network and overlay gene-set statistics
    set.seed(36) #set seed for reproducible layout
    
    
    #identify clusters
    grps = cluster_walktrap(gs_ovnet)
    #extract clustering results
    grps = groups(grps)
    #sort by cluster size
    grps = grps[order(sapply(grps, length), decreasing = TRUE)]
    #plot the top 12 clusters
    set.seed(36) #set seed for reproducible layout
    
    sendSweetAlert(session = session, title = "Notification",
                   btn_labels = NA,
                   text = "step 2/4 has been done", type = "warning",
                   closeOnClickOutside = FALSE , showCloseButton = FALSE)
    
  
    genes = names(geneList)
    gene_stats = unname(geneList)
    names(gene_stats) = genes
    
   
    #create independent plots
    set.seed(36) #set seed for reproducible layout
    p1 = plotMsigWordcloud(msigdb_mm, grps[1:6], type = 'Name')
    p2 = plotMsigNetwork(gs_ovnet, markGroups = grps[1:6], genesetStat = geneset_stats)
    p3 = plotGeneStats(gene_stats, msigdb_mm, grps[1:6]) +
      geom_hline(yintercept = 0, colour = 2, lty = 2)
    
    sendSweetAlert(session = session, title = "Notification",
                   btn_labels = NA,
                   text = "step 3/4 has been done", type = "warning",
                   closeOnClickOutside = FALSE , showCloseButton = FALSE)
    
    #combine using functions from ggpubr
    mydata$visse1 <- ggarrange(p1, p2, p3, ncol = 3, common.legend = TRUE, legend = 'bottom')
    
    sendSweetAlert(session = session, title = "Notification",
                   text = "Analysis is done!", type = "success",
                   closeOnClickOutside = TRUE , showCloseButton = TRUE)
  })

  #visse UI
  output$visseinput <- renderPlot({
    req(mydata$visse1)
    plot(mydata$visse1)
  })
    
  ## Fgsea enrichment analysis ##---------------------------------------
  observeEvent(input$buttonfg, {
    
    if(!exists('mydata$array')){sendSweetAlert(session = session, title = "Error", text = "Please upload your data first", type = "error")}
    
    
    req(mydata$array, input$fgnumber)
    
    sendSweetAlert(session = session, title = "Notification",
                   btn_labels = NA,
                   text = "Analysis in Progress", type = "warning",
                   closeOnClickOutside = FALSE , showCloseButton = FALSE)
    
    fgseaRes <- fgsea(pathways =examplePathways,
          stats = mydata$array,
          minSize = 15,
          maxSize = 500)
    
    topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=as.numeric(input$fgnumber)), pathway]
    topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=as.numeric(input$fgnumber)), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    
    
    mydata$fgseaplot <- (plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, 
                  gseaParam=0.5, render = FALSE))
    
    f1 <<- mydata$fgseaplot
    
    
    sendSweetAlert(session = session, title = "Notification", 
                   text = "Analysis has completed!", type = "success",
                   closeOnClickOutside = TRUE, showCloseButton = TRUE)
    
    
  })
  
  output$fgseainput <- renderPlot({
    
    req(mydata$fgseaplot)
    plot(mydata$fgseaplot)
  })
  
  
  #Gene Ontology Analysis -----------------------------------------------
  observeEvent(input$buttonpa, {
    
    if(!exists('mydata$go_array')){sendSweetAlert(session = session, title = "Error", text = "Please upload your data first", type = "error")}
    
    req(mydata$go_array, input$pgoinput)
    
    sendSweetAlert(session = session, title = "Notification",
                   btn_labels = NA,
                   text = "Analysis in Progress", type = "warning",
                   closeOnClickOutside = FALSE , showCloseButton = FALSE)
    
    gene_List <- mydata$go_array
    gene <-  names(gene_List)[abs(gene_List) > 1]
    
    
    
    CLP$ego <- enrichGO(gene      = gene,
                    universe      = names(gene_List),
                    OrgDb         = org.Mm.eg.db::org.Mm.eg.db,
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = as.numeric(input$pgoinput),
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
    
    mydata$CLPG <-goplot(CLP$ego)
    mydata$CLPT <- (CLP$ego@result)
  
    sendSweetAlert(session = session, title = "Notification", 
                   text = "Analysis has completed!", type = "success",
                   closeOnClickOutside = TRUE, showCloseButton = FALSE)
    
  })
    
    #Gene pathway UI-1
  
  output$goinput <- renderPlot({ 
    
    req(mydata$CLPG)
    plot(mydata$CLPG)
    })
  
  #Gene pathway  UI-2
  output$gotable <- renderTable({
    
    req(mydata$CLPT)
  
    plot((mydata$CLPT))
  })
    
  # Creates an overview of our mapped data
  output$contents <- DT::renderDataTable(server = FALSE, {
    
    req(mydata$protdf)
    
    df <- mydata$protdf
    
    df[,2] <- formatC(as.data.frame(df)[,2], digits = 2, format = 'e')
  
    df[,3] <- formatC(as.data.frame(df)[,3], digits = 3) 
  
    return(DT::datatable(df, extensions = 'Buttons', rownames = FALSE, options = list(dom = 'tpB', fixedColumns = TRUE, autoWidth = FALSE, pagingType = "numbers", scrollX = T, buttons = c('copy', 'csv', 'excel','pdf'))) %>%
      DT::formatStyle(columns = colnames(data), fontSize = '80%'))
    
  })
  

  
  # Render a histogram of pvalues----------------------------
  observeEvent(input$tabs, {
    
    if (req(input$tabs) == "PVALUE") {
      
      if(!exists('mydata$protdf')){sendSweetAlert(session = session, title = "Error", text = "Please upload your data first", type = "error")}
      
      req(isolate(mydata$protdf))
      
      sendSweetAlert(session = session, title = "Notification", btn_labels = NA, text = "Analysis in Progress", type = "warning", closeOnClickOutside = FALSE , showCloseButton = FALSE)
      
      
      # A Density plot for P value------------------------------
      mydata$pplot1 <- plotly::ggplotly(ggplot(isolate(mydata$protdf),aes(y=-log10(pvalue), x=fold_change, label = ID)) + geom_point(col = "blue", alpha = 0.2) + theme_bw())
  
      p1 <- ggplot(isolate(mydata$protdf), aes(x = pvalue)) + geom_histogram(bins = 40, fill = "blue", alpha = 0.2) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme_bw()
      p2 <- ggplot(isolate(mydata$protdf), aes(x = fold_change)) + geom_histogram(bins = 100, fill = "blue", alpha = 0.2) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme_bw()
      
      mydata$pplot2 <- (p1 + p2)
      
      sendSweetAlert(session = session, title = "Notification", text = "Analysis has completed!", type = "success", closeOnClickOutside = TRUE, showCloseButton = FALSE)
      
      } 
  
  })
  
  #P value UI-1
  output$density_pvalue <- plotly::renderPlotly({ 
    req(mydata$pplot1)
    return(mydata$pplot1)
    })
  
  #P value UI-2
  output$histogram_pvalue <- renderPlot({ 
    req(mydata$pplot2)
    plot(mydata$pplot2)
    })
  
  

  

  


  
  # Render a bargraph for drug interaction data----------------------------
  observeEvent(input$tabs, {
    
    if (req(input$tabs) == "dugi") {
      
      if(!exists('mydata$protdf')){sendSweetAlert(session = session, title = "Error", text = "Please upload your data first", type = "error")}
      
      req(isolate(mydata$protdf), isolate(mydata$CP_summary))
      
      sendSweetAlert(session = session, title = "Notification", btn_labels = NA, text = "Analysis in Progress", type = "warning", closeOnClickOutside = FALSE , showCloseButton = FALSE)
      
      mydata$plot_anno_dgi <- (mydata$protdf %>% dplyr::select(ID, drug_name) %>% mutate(drug_name = !is.na(drug_name)) %>% group_by(ID) %>% summarise(n = sum(drug_name)) %>% mutate(n = (n != 0)) %>% ggplot(., aes(x = n)) + geom_bar(stat = "count", fill = "blue", alpha = 0.2, col = "black") + theme_bw())
      
      p1 <- (mydata$protdf %>% dplyr::select(ID, CP_loc) %>% mutate(CP_loc = !is.na(CP_loc)) %>% group_by(ID) %>% summarise(n = sum(CP_loc)) %>% mutate(n = (n != 0)) %>% ggplot(., aes(x = n)) + geom_bar(stat = "count", fill = "blue", alpha = 0.2, col = "black") + theme_bw())
      p2 <- ggplot(mydata$CP_summary, aes(x=frequency, y=reorder(CP_loc, frequency))) + geom_bar(stat = "identity", fill = "blue", alpha = 0.2, col = "black") + theme_minimal() + labs(x = "Frequency", y = "")
      
      mydata$plot_anno_loc <- (p1 + p2)
      
      sendSweetAlert(session = session, title = "Notification", text = "Analysis has completed!", type = "success", closeOnClickOutside = TRUE, showCloseButton = FALSE)
      
    } 
})
  
  #annotation dgi
  output$bargraph_drug <- renderPlot({
    req(isolate(mydata$protdf))
    plot(mydata$plot_anno_dgi)
  })
  
  #annotation dgi
  output$bargraph_loc <- renderPlot({
    req(isolate(mydata$protdf))
    plot(mydata$plot_anno_loc)
  })
  

  
  
  
  
  
  #Making columns used for contigency tabele
 # x <- x %>% mutate(direction = case_when( (p < 0.05 & fold_change > 0) ~ "up", (q < 0.05 & cor < 0) ~ "down", TRUE ~ "NS"  ))
#
 
  # Case_when function
  #x5 %>% mutate(significant = case_when( (pvalue < 0.05 & for_change > 1 ) ~ "significant", 
                                           #TRUE ~ "not significant" ))
}


# Create Shiny app ----
shinyApp(ui, server)

