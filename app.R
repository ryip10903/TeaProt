library(shiny)
library(shinydashboard)
library(dplyr)
library(org.Mm.eg.db)
# Define UI for data upload app ----
ui <- dashboardPage(
  
  # App title ----
  dashboardHeader(title = "Uploading Files"),
  
  
  
  # Sidebar panel for inputs ----
  dashboardSidebar(
    
    # Input: Select a file ----
    fileInput("file", "Choose CSV File",
              multiple = FALSE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv"))
  ),
  
  # Main panel for displaying outputs ----
  dashboardBody(
    
    # Output: Data file ----
    tableOutput("contents"),
    downloadButton("dl_table", "Download your table here!"),
    
    # Output: Boxplot ----
    plotOutput("boxplot_fc"), 
    plotOutput("histogram_pvalue")
    
  )
  
)

# Define server logic to read selected file ----
server <- function(input, output) {
  
  # Create a reactiveValues object called mydata
  mydata <- reactiveValues()
  
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
    
    entrezids <- as.list(org.Mm.egALIAS2EG)
    names(df)[1] <- "ID"
    df <- cp_idconvert(df, cp_idtype(df$ID))
    #insert Jeff's functino (cp_idconverter)
      #changing uniprot id to gene names
    #provided_data <- read.csv(file='database/knocktf.csv')
    
    #df <- left_join(df, provided_data, by = c('genenames'='genenames'))
      #Q= why do we need to add $protdf onto mydata
     
    
    
    #provided_data <- read.csv(file='knocktf/differential expression of genes in all datasets.txt')
    #df <- left_join(df, provided_data, by = c('ID'='genenames'))
  
    mydata$protdf <- df 
    print(df)
  })
  
  
  
  
  # The table we make can use our reactive protdf dataframe
  output$contents <- renderTable({

    req(mydata$protdf)
    
    return(head(mydata$protdf)) 
    #Q=why do we need this line of command
    
  })
  # The table we make can use our reactive protdf dataframe
  output$contents <- renderTable({
    
    req(mydata$protdf)
    
    return(head(mydata$protdf)) 
    #Q=why do we need this line of command
    
  })
  
  # Render a histogram of pvalues
  output$histogram_pvalue <- renderPlot({
    
    req(mydata$protdf)
    
    x1 <<- mydata$protdf
    
    print(head(mydata$protdf))
    
    return(hist(mydata$protdf[,3]))
  })
  
  # Render a boxplot of fold-change
  output$boxplot_fc <- renderPlot({
    
    req(mydata$protdf)
    
    return(boxplot(mydata$protdf[,4]))    
  })
  
  # Downloader for the table
  output$dl_table <- cp_dl_table_csv(mydata$protdf, "annotated_data.csv")
  
  
}

# Create Shiny app ----
shinyApp(ui, server)

