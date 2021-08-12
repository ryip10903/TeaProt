library(shiny)
library(shinydashboard)
library(fgsea)
library(data.table) 
library(ggplot2)
library(org.Hs.eg.db)
library(shinyWidgets)
library(markdown)
library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(DT)
library(shinycssloaders)
library(msigdbr)
library(GSEABase)


# Define UI foror data upload app ----
ui <- dashboardPage( skin = 'black',
  
  # App title ----
  dashboardHeader(title = span('Uploading Files', style = "font-weight: bold")),
  
  
  
  # Sidebar panel for inputs ----
  dashboardSidebar(
    sidebarMenu(
      # Input: Select a file ----
      fileInput("file", "Choose CSV File",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      
      selectInput("species", "Choose Species",
                  c("human", "drosophila", "mouse")),
      
    #Start analysis button
    actionButton("button", 'Start Here!', style='font-weight:600' ),
      
    #sidebar Items
      menuItem("START", tabName = "start", icon = icon("bookmark")),
      
      menuItem("View Data", tabName = "view", icon = icon("bookmark")),
      
      
      menuItem("Analysis", tabName = "ANALYSIS", icon = icon("bookmark"),
              menuSubItem("Pvalue", tabName = "PVALUE", icon= icon("info-circle")),
              menuSubItem("Fchange", tabName = "forchange"),
              menuSubItem("Drug Interaction", tabName = "dugi"),
              menuSubItem("Gene Location", tabName = "genel"),
              menuSubItem("Gene Pathway", tabName = "genep"),
              menuSubItem("Fgsea analysis", tabName = "fgseaa"),
              menuSubItem('vissE analysis', tabName = 'vvss')
              ),
      
     menuItem("youtube", icon= icon("file-code-o"),
              href= "https://www.youtube.com/")


    ),
    sidebarMenuOutput("menu")
  ),
  
  # Main panel for displaying outputs ----
  dashboardBody(
    
    
    
  #UI layout for each tab
    tabItems(
      
      
      
      tabItem(tabName = "start",
              div(class = "jumbotron", style="background-image: url(dna-banner.svg); 
                  background-size: cover;", HTML("<center><h1>Welcome to TeaProt!</h1></center>"),
                  HTML("<center><p>For the annotation of Protein/transcript data.</p></center>")),
                  includeMarkdown("README.md")),
      
      tabItem(tabName = "view", h2("View Your Data Here"), DT::DTOutput("contents")),
      
      tabItem(
        actionButton('buttonp', 'Start Analysis!'),
              tabName = "PVALUE",
              h2("P-value Analysis"),  
        fluidRow(
          box( title = 'Histogram', plotOutput("histogram_pvalue") %>% withSpinner()),
          box( title = 'Density Plot', plotOutput("density_pvalue") %>% withSpinner() 
               ))),
              
    
      tabItem(tabName = "ANALYSIS"),
      
      tabItem(
        actionButton('buttond', 'Start Analysis'),#style = 'background-color:MistyRose;'),
        tabName = "dugi", 
              h2("Drug Interaction Analysis"), plotOutput("bargraph_drug") %>% withSpinner()),
      
      tabItem(
        actionButton('buttong', 'Start Analysis'),
        tabName = "genel", 
              h2("Gene Location Analysis"), 
              fluidRow( box(width =12, plotOutput('CP_location') %>% withSpinner()))
                       ),
    
      tabItem(
        actionButton('buttonpa', 'Start Analysis'),
        selectInput('pgoinput', 'Choose P_value cutoff',
                    c(0.001, 0.01, 0.05)),
        tabName = "genep",
              h2("Gene Pathway Analysis"), 
              fluidRow(column(12, box(plotOutput('goinput')%>% withSpinner(), width =12)),
              column(12, box(tableOutput('gotable'), width = 12))),
        ),
              
      tabItem(
        actionButton('buttonfg', 'Start Analysis'), 
        selectInput("fgnumber", "Choose Number of Pahways to be seen",
                    c(10, 20, 30)),
        tabName = "fgseaa", h2("Fgsea Analysis"), 
        fluidRow(column(12,
          box(width = 12, plotOutput("fgseainput") %>% withSpinner() )))),
      
      tabItem(
        actionButton('buttonf', 'Start Analysis'),
        tabName = "forchange", h2("Foldchange Analysis"), 
              fluidRow(column(6, box(title="Boxplot",plotOutput("boxplot_fc") %>% withSpinner(),width = 12)), 
                       column(6,plotOutput('zplot_fc',height = 500) %>% withSpinner())),
              downloadButton("dl_table", "Download your table here!")),
      
      tabItem(
        actionButton('buttonv', 'Start Analysis'),
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
    
    x1 <<- df
  
  #converting the uploaded column names to match our commands
    names(df)[1] <- "ID"
    names(df)[2] <- 'pvalue'
    names(df)[3] <- 'fold_change'
    x2 <<- df
  #pre-made function
  #converting gene names to entrez id
    df <- cp_idconvert(df, cp_idtype(df$ID))
    
    x3 <<- df
    
  # joining the user's data with mouse gene id
    df <- left_join(df, entrezmapping("mouse"))
    
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
  
  #Genelocation analysis
  observeEvent(input$buttong, {
    
    #sendSweetAlert(session = session, title = "Notification", 
                   #text = "Analysis is now in progress ", type = "waning",
                   #btn_labels = NA,
                   #closeOnClickOutside = FALSE, showCloseButton = FALSE)
    
    output$CP_location <- renderPlot ({
    
    req(mydata$CP_summary)
    
    return(ggplot(mydata$CP_summary, aes(x=frequency, y=reorder(CP_loc, frequency)))
           + geom_bar(stat = "identity") + theme_minimal() + labs(x = "Frequency", y = ""))
      
      #sendSweetAlert(session = session, title = "Notification", 
                     #text = "Your analysis has completed!", type = "success",
                     #closeOnClickOutside = TRUE, showCloseButton = FALSE)
  })})

    
  ## visse analysis ##
  
  observeEvent(input$buttonv, {
    
    output$visseinput <- renderPlot({
    
    # Create input data 
    d <- mydata$protdf[,c(1,3)] %>% filter(!is.nan(ID) & !is.na(fold_change)) %>% dplyr::rename(gene = 1, fc = 2) %>% group_by(gene) %>% summarize(fc = mean(fc))
    
    geneList <- d[,2] %>% unlist
    names(geneList) <- as.character(d[,1] %>% unlist %>% stringr::str_to_title())
    geneList <- base::sort(geneList, decreasing = TRUE)
    
    
    #load the MSigDB from the msigdb package
    msigdb_mm = msigdb.v7.2.mm.SYM()
    #append KEGG gene-sets
    msigdb_mm = appendKEGG(msigdb_mm)
    #dplyr::select h, c2, and c5 collections (recommended)
    msigdb_mm = subsetCollection(msigdb_mm, c('h', 'c2', 'c5'))
    
    # Create an input gene set for the GSEA function
    genedb <- geneIds(msigdb_mm)
    msig_db_gsea <- data.frame(gs_name = rep(names(genedb), lengths(genedb)), gene_symbol = unlist(genedb, use.names=TRUE))
    
    # Run GSEA
    set.seed(36)
    
    em <- GSEA(geneList, TERM2GENE = msig_db_gsea, pvalueCutoff = 0.05)
    geneset_res = em@result$Description
    
    #create a GeneSetCollection using the gene-set analysis results
    geneset_gsc = msigdb_mm[geneset_res]
    
    #compute gene-set overlap
    gs_ovlap = computeMsigOverlap(geneset_gsc, thresh = 0.25)
    #create an overlap network
    gs_ovnet = computeMsigNetwork(gs_ovlap, msigdb_mm)
    #plot the network
    set.seed(36) #set seed for reproducible layout
    plotMsigNetwork(gs_ovnet)
    
    ## -----------------------------------------------------------------------------
    #simulate gene-set statistics
    geneset_stats = em@result$NES
    names(geneset_stats) = geneset_res
    head(geneset_stats)
    
    #plot the network and overlay gene-set statistics
    set.seed(36) #set seed for reproducible layout
    plotMsigNetwork(gs_ovnet, genesetStat = geneset_stats)
    
    #identify clusters
    grps = cluster_walktrap(gs_ovnet)
    #extract clustering results
    grps = groups(grps)
    #sort by cluster size
    grps = grps[order(sapply(grps, length), decreasing = TRUE)]
    #plot the top 12 clusters
    set.seed(36) #set seed for reproducible layout
    plotMsigNetwork(gs_ovnet, markGroups = grps[1:6], genesetStat = geneset_stats)
    
    ## -----------------------------------------------------------------------------
    #compute and plot the results of text-mining
    #using gene-set Names
    plotMsigWordcloud(msigdb_mm, grps[1:6], type = 'Name')
    #using gene-set Short descriptions
    plotMsigWordcloud(msigdb_mm, grps[1:6], type = 'Short')
    
    set.seed(36)
    
    genes = names(geneList)
    gene_stats = unname(geneList)
    names(gene_stats) = genes
    
    #plot the gene-level statistics
    plotGeneStats(gene_stats, msigdb_mm, grps[1:6]) +
      geom_hline(yintercept = 0, colour = 2, lty = 2)
    
    #create independent plots
    set.seed(36) #set seed for reproducible layout
    p1 = plotMsigWordcloud(msigdb_mm, grps[1:6], type = 'Name')
    p2 = plotMsigNetwork(gs_ovnet, markGroups = grps[1:6], genesetStat = geneset_stats)
    p3 = plotGeneStats(gene_stats, msigdb_mm, grps[1:6]) +
      geom_hline(yintercept = 0, colour = 2, lty = 2)
    
    #combine using functions from ggpubr
    return(ggarrange(p1, p2, p3, ncol = 3, common.legend = TRUE, legend = 'bottom'))
  })})

  
  ## Fgsea enrichment analysis ##
  observeEvent(input$buttonfg, {
    
    sendSweetAlert(session = session, title = "Notification", 
                   text = "Analysis in Progress", type = "warning",
                   closeOnClickOutside = TRUE, showCloseButton = FALSE)
    
    output$fgseainput <- renderPlot({
    
    req(mydata$array, input$fgnumber)
    
    fgseaRes <- fgsea(pathways =examplePathways,
          stats = mydata$array,
          minSize = 15,
          maxSize = 500)
    
    topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=as.numeric(input$fgnumber)), pathway]
    topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=as.numeric(input$fgnumber)), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    
    sendSweetAlert(session = session, title = "Notification", 
                   text = "Your analysis has completed!", type = "success",
                   closeOnClickOutside = TRUE, showCloseButton = FALSE)
    
    return(plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, 
                  gseaParam=0.5))
    
    
    sendSweetAlert(session = session, title = "Notification", 
                   text = "Analysis has completed!", type = "success",
                   closeOnClickOutside = TRUE, showCloseButton = FALSE)
    
  })})
  
  #Gene Ontology Analysis
  observeEvent(input$buttonpa, {
    
    sendSweetAlert(session = session, title = "Notification", 
                   text = "Analysis in progress!", type = "warning",
                   closeOnClickOutside = TRUE, showCloseButton = FALSE)
    
    output$goinput <- renderPlot({
    
    req(mydata$go_array, input$pgoinput)
    #add message "starting analysis", disable the function of closing
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
    
    return(goplot(CLP$ego))
  })
    
    
  #Show the Gene Ontology result table
    output$gotable <- renderTable({
      
      req(CLP$ego@result)
    z3 <<- CLP$ego@result
    return((CLP$ego@result))
    
    sendSweetAlert(session = session, title = "Notification", 
                   text = "Analysis has completed!", type = "success",
                   closeOnClickOutside = TRUE, showCloseButton = FALSE)
    
  })})
  
  
  
    
  # Creates an overview of our mapped data
  output$contents <- DT::renderDT({
    
    req(mydata$protdf)
    
    df <- mydata$protdf
    
    df[,2] <- formatC(df[,2], digits = 2, format = 'e')
    df[,3] <- formatC(df[,3], digits = 3) 
  
    return(DT::datatable(df, options = list(scrollX=TRUE)))
    
  })
  

  
  # Render a histogram of pvalues
  observeEvent(input$buttonp, {
    
    #sendSweetAlert(session = session, title = "Notification", 
                  #button = FALSE,
                  #text = "Analysis in Progress ", type = "warning",
                  #closeOnClickOutside = FALSE, showCloseButton = FALSE)
    
    output$histogram_pvalue <- renderPlot({
    
    req(mydata$protdf)
    
    return(hist(mydata$protdf[,2]))
      
      
  })
  
  # A Density plot for P value
  output$density_pvalue <- renderPlot({
    
    req(mydata$protdf)
    
    return(ggplot(mydata$protdf,aes(y=-log10(pvalue), x=fold_change)) + geom_point() +
             ggtitle("Foldchange VS P_value") + theme(plot.title = element_text(
               hjust = 0.5, size = 15, face = 'bold', color = 'blue')) )

  })})
  # Render a boxplot of fold-change
  observeEvent(input$buttonf, {
    
   
    
   output$boxplot_fc <- renderPlot({
    
    req(mydata$protdf)
    
    return(boxplot(mydata$protdf[,3]))    
  })
  
  # Render a cool plot of fold-change
  output$zplot_fc <- renderPlot({
    
    req(mydata$protdf)
    
    return(ggplot(mydata$protdf, aes(y=fold_change, x=pvalue)) + geom_polygon()
                  )
    #sendSweetAlert(session = session, title = "Notification", 
                   #button = FALSE,
                   #text = "Analysis in Progress ", type = "success",
                   #closeOnClickOutside = FALSE, showCloseButton = FALSE)
  })})
  
  # Render a bargraph for drug interaction data
  
  observeEvent(input$buttond, {
    
    #sendSweetAlert(session = session, title = "Notification", 
                   #button = FALSE,
                   #text = "Analysis in Progress ", type = "warning",
                   #closeOnClickOutside = FALSE, showCloseButton = FALSE)
    
    output$bargraph_drug <- renderPlot({
    
    req(mydata$protdf)
    
    return(ggplot(mydata$protdf, aes(x= !is.na(drug_name))) + geom_bar())
      
      #sendSweetAlert(session = session, title = "Notification", 
                    # button = FALSE,
                     #text = "Analysis has Completed ", type = "success",
                     #closeOnClickOutside = FALSE, showCloseButton = FALSE)
  })})
  

 
  # Downloader for the table
  output$dl_table <- cp_dl_table_csv(mydata$protdf, "annotated_data.csv")
  
  # Case_when function
  #x5 %>% mutate(significant = case_when( (pvalue < 0.05 & for_change > 1 ) ~ "significant", 
                                           #TRUE ~ "not significant" ))
}


# Create Shiny app ----
shinyApp(ui, server)

