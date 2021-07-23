library(shiny)
library(shinydashboard)
library(fgsea)
library(data.table) 
library(ggplot2)
library(org.Hs.eg.db)
library(shinyWidgets)
library(markdown)
library(clusterProfiler)
library(dplyr)


# Define UI foror data upload app ----
ui <- dashboardPage( skin = 'purple',
  
  # App title ----
  dashboardHeader(title = "Uploading Files"),
  
  
  
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
    actionButton("button", "Start your analysis!"),
      
    #sidebar Items
      menuItem("START", tabName = "start", icon = icon("bookmark")),
      
      menuItem("View Data", tabName = "view", icon = icon("bookmark")),
      
      
      menuItem("Analysis", tabName = "ANALYSIS", icon = icon("bookmark"),
              menuSubItem("Pvalue", tabName = "PVALUE", icon= icon("info-circle")),
              menuSubItem("Foldchange", tabName = "forchange"),
              menuSubItem("Drug Interaction", tabName = "dugi"),
              menuSubItem("Gene Location", tabName = "genel"),
              menuSubItem("Gene Pathway", tabName = "genep"),
              menuSubItem("Fgsea analysis", tabName = "fgseaa")),
      
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
      
      tabItem(tabName = "view", h2("View Your Data Here"), tableOutput("contents")),
      
      tabItem(tabName = "PVALUE",
              h2("P-value Analysis",  plotOutput("histogram_pvalue"),
                 plotOutput("density_pvalue"))
              ),
      tabItem(tabName = "ANALYSIS"),
      tabItem(tabName = "dugi", 
              h2("Drug Interaction Analysis"), plotOutput("bargraph_drug")),
      tabItem(tabName = "genel", 
              h2("Gene Location Analysis"), plotOutput('CP_location')),
      tabItem(tabName = "genep",
              h2("Gene Pathway Analysis"), plotOutput('goinput'),
              tableOutput('gotable')),
      tabItem(tabName = "fgseaa", h2("Fgsea Analysis"), plotOutput("fgseainput")),
      tabItem(tabName = "forchange", h2("Foldchange Analysis"), 
              fluidRow(column(6, box(title="Boxplot", status = 'warning',
                                     plotOutput("boxplot_fc"),
                                     width = 12)), column(6,plotOutput('zplot_fc'
                                                                       ,height = 500))),
              downloadButton("dl_table", "Download your table here!"))
      
    ),
    
   
    
    
    
    
  ))

# Define server logic to read selected file ----
server <- function(input, output, session) {
  output$menu <- renderMenu({
    sidebarMenu(
      menuItem("Menu item", icon = icon("calender"))
    )
  })
  
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
    
    #db_hpa_locoverlap <<- db_hpa %>% select(3) %>% dplyr::mutate(loc2 = CP_loc) %>% `colnames<-`(c("loc1", "loc2")) %>% distinct %>% tidyr::expand(loc1 = .$loc1, loc2 = .$loc2) %>% mutate(overlap.loc = mapply(function(x, y) paste(intersect(x, y), collapse=";"), strsplit(.$loc1, ";"), strsplit(.$loc2, ";")) ) 
    
    #adding a progress bar
    dat <- data.frame(x = numeric(0), y = numeric(0))
    
    withProgress(message = 'Making plot', value = 0, {
      # Number of times we'll go through the loop
      n <- 10
      
      for (i in 1:n) {
        # Each time through the loop, add another row of data. This is
        # a stand-in for a long-running computation.
        dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("Doing part", i))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
    })
    
    plot(dat$x, dat$y)
    
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
  #converting row1's name to ID
    names(df)[1] <- "ID"
    
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
                   text = "Your analysis has completed! Click on View Data and Analysis ", type = "success",
                   closeOnClickOutside = TRUE, showCloseButton = FALSE)
    
  })
  
  #Genelocation analysis
  output$CP_location <- renderPlot ({
    
    req(mydata$CP_summary)
    
    return(ggplot(mydata$CP_summary, aes(x=frequency, y=reorder(CP_loc, frequency)))
           + geom_bar(stat = "identity") + theme_minimal() + labs(x = "Frequency", y = ""))
  })

    
  
  
   #Fgsea enrichment analysis
  output$fgseainput <- renderPlot({
    
    req(mydata$array)
    
    fgseaRes <- fgsea(pathways =examplePathways,
          stats = mydata$array,
          minSize = 15,
          maxSize = 500)
    
    topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
    topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    
    return(plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, 
                  gseaParam=0.5))
    
    
    
    
  })
  
  #Gene Ontology Analysis
  output$goinput <- renderPlot({
    
    req(mydata$go_array)
    #add message "starting analysis", disable the function of closing
    geneList <- mydata$go_array
    gene <-  names(geneList)[abs(geneList) > 1]
    
    
    
    CLP$ego <- enrichGO(gene          = gene,
                    universe      = names(geneList),
                    OrgDb         = org.Mm.eg.db::org.Mm.eg.db,
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
    
    return(goplot(CLP$ego))
  })
    #message = analysis completed
    
  #Show the Gene Ontology result table
    output$gotable <- renderTable({
      
      req(CLP$ego)
    z3 <<- CLP$ego
    return((CLP$ego))
    
    
  })
  
  
  # Creates an overview of our mapped data
  output$contents <- renderTable({

    req(mydata$protdf)
  
    return(mydata$protdf[1:10,])
   
    
  })
  

  
  # Render a histogram of pvalues
  output$histogram_pvalue <- renderPlot({
    
    req(mydata$protdf)
    
    return(hist(mydata$protdf[,2]))
  })
  
  # A Density plot for P value
  output$density_pvalue <- renderPlot({
    
    req(mydata$protdf)
    
    return(ggplot(mydata$protdf,aes(y=pvalue, x=for_change)) + geom_point() +
             ggtitle("Foldchange VS P_value") + theme(plot.title = element_text(
               hjust = 0.5, size = 15, face = 'bold', color = 'blue')
               
             )) 
  })
  # Render a boxplot of fold-change
  output$boxplot_fc <- renderPlot({
    
    req(mydata$protdf)
    
    return(boxplot(mydata$protdf[,3]))    
  })
  
  # Render a cool plot of fold-change
  output$zplot_fc <- renderPlot({
    
    req(mydata$protdf)
    
    return(ggplot(mydata$protdf, aes(y=for_change, x=pvalue)) + geom_polygon()
                  )
  })
  
  # Render a bargraph for drug interaction data
  
  output$bargraph_drug <- renderPlot({
    
    req(mydata$protdf)
    
    return(ggplot(mydata$protdf, aes(x= !is.na(drug_name))) + geom_bar())
  })
  
  # Render a graph for protein localization data
  output$bargraph_location <- renderPlot({
    
  })
 
  # Downloader for the table
  output$dl_table <- cp_dl_table_csv(mydata$protdf, "annotated_data.csv")
  
  # Case_when function
  #x5 %>% mutate(significant = case_when( (pvalue < 0.05 & for_change > 1 ) ~ "significant", 
                                           #TRUE ~ "not significant" ))
}


# Create Shiny app ----
shinyApp(ui, server)

