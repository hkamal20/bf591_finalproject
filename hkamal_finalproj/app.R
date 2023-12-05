#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
library(bslib)
library(ggplot2)
library(colourpicker) 
library(tidyverse)

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Haniya Kamal BF591 Final Project"),
  sidebarLayout(
    sidebarPanel(
      fileInput('csv_file', label = paste0('Upload your metadata matrix'), accept = '.csv'),
      fileInput('csv_file2', label = paste0('Upload your counts matrix'), accept = '.csv'),
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Samples",
          tabsetPanel(
            tabPanel("Summary", paste0("Summary table of Dataset:"), tableOutput("summary")), #Tab with a summary of the table that includes a summary of the type and values in each column
            tabPanel("Table", DTOutput("table")), #Tab with a data table displaying the sample information, with sortable columns
            tabPanel("Plots", #Tab with histograms, density plots, or violin plots of continuous variable3
              plotOutput("plot1"),
              plotOutput("plot2"),
              plotOutput("plot3")
            )
          )
        ),
        tabPanel("Counts",
          tabsetPanel(
            tabPanel("Counts Table", tableOutput("countstable")),
            tabPanel("Scatterplots",
              plotOutput("plot4"), # Median count vs variance
              plotOutput("plot5")  # Median count vs number of zeros
                   ),
            tabPanel("Heatmap",
              plotOutput("heatmap")  # Clustered heatmap of counts remaining after filtering
                   ),
            tabPanel("PCA plot",
              plotOutput("pcaplot")
            )
          )
        ),
        tabPanel("Differential Expression",
            tabsetPanel(
              tabPanel("Differential Expression Results", tableOutput("table3")), #sortable table displaying differential expression results
                #Tab with content similar to that described in [Assignment 7]
          )
        ),
        tabPanel("Gene Set Enrichment Analysis"),
      )
    )
  )
)

#
server <- function(input, output) {
  options(shiny.maxRequestSize = 30*1024^2)
  
  load_data <- reactive({ #CSV file can be uploaded for the metadata
    req(input$csv_file)
    data <- input$csv_file
    if (is.null(data)) {
      return(NULL)
    }
    d<-read.csv(data$datapath, header = TRUE)
    return(d)
  })
  
  load_counts <- reactive({ #CSV file can be uploaded for the counts matrix
    req(input$csv_file2)
    data2 <- input$csv_file2
    if (is.null(data2)) {
      return(NULL)
    }
    d2<-read.csv(data2$datapath, header = TRUE)
    return(d2)
  })
  
  datatable1 <- function(dataf) { #making the datatable for samples tab
    df <- dataf 
    return(df)
  }
  
  summarytable <- function(dataf) {
    notsum <- datatable1(dataf)
    df <- notsum[, !(names(dataf) %in% c("Sample.Name", "Sample.Geo.Accession"))]
    
    #function to calculate Mean (sd) or Distinct Values
    summarization <- function(x) {
      # Remove NAs and then calculate mean and sd for numeric columns
      x <- x[!is.na(x)]
      
      if (is.numeric(x)) {
        mean_val <- mean(x)
        sd_val <- sd(x)
        return(paste(round(mean_val, 2), " (+/- ", round(sd_val, 2), ")"))
      } else {
        # Count distinct values for non-numeric columns
        distinct_values <- length(unique(x))
        return(paste(distinct_values, " distinct values"))
      }
    }
    
    # Create a logical vector for character columns
    is_character <- sapply(df, is.character)
    
    # Conditionally set values based on data type
    type_values <- ifelse(is_character, "character", sapply(df, class))
    
    # Create the final data frame
    d_sum <- data.frame(
      Columns = names(df),
      Type = type_values,
      Mean_or_Distinct_Value = ifelse(is_character, "Neurologically Normal, Huntington's Disease", sapply(df, summarization))
    )
    
    return(d_sum)
  }

  output$summary <- renderTable({ #pulling the summary table from the function, part of samples tabset
    summary_df <- load_data()
    summarytable(summary_df)
    
  })
  
  output$table <- renderDT({ #pulling the table from the Samples tabset
    data2 <- load_data()             #table will be sortable
    datatable1(data2) #loading the sample table 
  })
  
  output$plot1 <- renderPlot({
    violin_data <- load_data()
    plot <- ggplot(violin_data) +
      geom_violin(aes(x=!!sym('Diagnosis'),y=!!sym('Age.of.Death'),fill=Diagnosis))
    return(plot)
    
  })
  
  output$countstable <- renderTable({ #show the counts data in the counts matrix tabset
    counts <- load_counts()
    
  })
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
