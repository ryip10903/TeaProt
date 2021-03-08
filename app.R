library(shiny)
library(shinydashboard)
library(dplyr)

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
    plotOutput('contents'),
    # Output: Data file ----
    tableOutput("contents")
    
  )
  
)

# Define server logic to read selected file ----
server <- function(input, output) {
  
  # Create a reactiveValues object called mydata
  mydata <- reactiveValues()
  
  # The following code runs when a file is uploaded
  # We want to load data, and annotate it with out databases here
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
    
    
    
    #provided_data <- read.csv(file='database/knocktf.csv')
    
    #df <- left_join(df, provided_data, by = c('genenames'='genenames'))
      #Q= why do we need to add $protdf onto mydata
     
    
    
    provided_data <- read.csv(file='database/knocktf.csv')
    df <- left_join(df, provided_data, by = c('genenames'='genenames'))
    mydata$protdf <- df 
  })
  
  
  
  # The table we make can use our reactive protdf dataframe
  output$contents <- renderTable({
    boxplot(
    req(mydata$protdf),
    
    return(head(mydata$protdf)) )
    #Q=why do we need this line of command
    
  })
  
}

# Create Shiny app ----
shinyApp(ui, server)

