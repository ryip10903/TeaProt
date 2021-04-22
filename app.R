library(shiny)
library(shinydashboard)
library(fgsea)
library(data.table) 
library(ggplot2)
library(dplyr)

# Define UI foror data upload app ----
ui <- dashboardPage(
  
  # App title ----
  dashboardHeader(title = "Uploading Files"),
  
  
  
  # Sidebar panel for inputs ----
  dashboardSidebar(
    sidebarMenu(
      menuItem("START", tabName = "start", icon = icon("start")),
    # Input: Select a file ----
    fileInput("file", "Choose CSV File",
              multiple = FALSE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv"))
    )
  ),
  
  # Main panel for displaying outputs ----
  dashboardBody(
    tabItems(
      tabItem(tabName = "start")
    ),
    # Output: Data file ----
    tableOutput("contents"),
    downloadButton("dl_table", "Download your table here!"),
    
    # Output: Boxplot ----
    plotOutput("boxplot_fc"), 
    plotOutput("histogram_pvalue"),
    
    # Output: fgsea ----
    plotOutput("fgseainput")
    
  )
  
)

# Define server logic to read selected file ----
server <- function(input, output) {
  
  # Create a reactiveValues object called mydata
  mydata <- reactiveValues()
  forout_reactive <- reactiveValues()
  
  # The following code runs when a file is uploaded
  # We want to load data, and annotate it with our databases here
  # We first load the data based on its extension (.csv, .txt, .xlsx)
  # Then we load our databases and join the user data with out provided databases
  # Finally we save the df as a reactive object
  observeEvent(input$file, {
    
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
    return(head(mydata$protdf)) 
    #Q=why do we need this line of command
    
  })
  

  
  # Render a histogram of pvalues
  output$histogram_pvalue <- renderPlot({
    
    req(mydata$protdf)
    
    return(hist(mydata$protdf[,2]))
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

