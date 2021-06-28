library(shiny)
library(shinydashboard)
library(fgsea)
library(data.table) 
library(ggplot2)
library(dplyr)
library(shinyWidgets)
library(markdown)
# Define UI foror data upload app ----
ui <- dashboardPage(
  
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
              menuSubItem("Forchange", tabName = "forchange"),
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
    #fluidRow()
    
    
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
              h2("Drug interaction Analysis")),
      tabItem(tabName = "genel", 
              h2("Gene location Analysis")),
      tabItem(tabName = "genep",
              h2("Gene pathway Analysis")),
      tabItem(tabName = "fgseaa", h2("Fgsea Analysis"), plotOutput("fgseainput")),
      tabItem(tabName = "forchange", h2("For change Analysis"),  plotOutput("boxplot_fc"),
              downloadButton("dl_table", "Download your table here!"))
    ),
    # Output: Data file ----
    
    
    
    # Output: Boxplot ----
   
   
    
    # Output: fgsea ----
    
    
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
  
  # The following code runs when a file is uploaded
  # We want to load data, and annotate it with our databases here
  # We first load the data based on its extension (.csv, .txt, .xlsx)
  # Then we load our databases and join the user data with out provided databases
  # Finally we save the df as a reactive object
  observeEvent(input$button, {
    
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
    
    x4 <<- df
  ## this step for
    array <- df %>% select(3) %>% unlist()
    names(array) <- df %>% select(EntrezGeneID) %>% unlist()
    array <- sort(array)
    
    x5 <<- array
    
   
#attaching df to "mydata" reactive value
    mydata$protdf <- df 
    mydata$fgsea_array <- array 
    
    print(df)
    
    sendSweetAlert(session = session, title = "Notification", 
                   text = "Your analysis has completed! Click on View Data and Analysis ", type = "success",
                   closeOnClickOutside = TRUE, showCloseButton = FALSE)
    
  })
  
   #Fgsea enrichment analysis
  output$fgseainput <- renderPlot({
    
    req(mydata$fgsea_array)
    
    fgseaRes <- fgsea(pathways =examplePathways,
          stats = mydata$fgsea_array,
          minSize = 15,
          maxSize = 500)
    
    topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
    topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    
    return(plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, 
                  gseaParam=0.5))
    
    
    
    
  })
  

  
  
  # The table we make can use our reactive protdf dataframe
  output$contents <- renderTable({

    req(mydata$protdf)
  #returns the epitome of uploaded value
    return(mydata$protdf[1:10,])
    #Q=why do we need this line of command
    
  })
  

  
  # Render a histogram of pvalues
  output$histogram_pvalue <- renderPlot({
    
    req(mydata$protdf)
    
    return(hist(mydata$protdf[,2]))
  })
  
  # A Density plot for P value
  output$density_pvalue <- renderPlot({
    
    req(mydata$protdf)
    
    return(ggplot(mydata$protdf,aes(y=pvalue, x=for_change)) + geom_point() )
  })
  # Render a boxplot of fold-change
  output$boxplot_fc <- renderPlot({
    
    req(mydata$protdf)
    
    return(boxplot(mydata$protdf[,3]))    
  })
  
  # Downloader for the table
  output$dl_table <- cp_dl_table_csv(mydata$protdf, "annotated_data.csv")
  
  
}

# Create Shiny app ----
shinyApp(ui, server)

