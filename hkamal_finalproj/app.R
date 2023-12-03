#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Haniya Kamal Bf591 Final Project"),
  sidebarLayout(
    sidebarPanel(
      fileInput('csv_file', label = paste0('Upload your file')),
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Samples",
          tabsetPanel(
            tabPanel("Summary", verbatimTextOutput("summary")),
            tabPanel("Table", tableOutput("table")),
            tabPanel("Plots",
              plotOutput("plot1"),
              plotOutput("plot2"),
              plotOutput("plot3")
            )
          )
        ),
        tabPanel("Counts",
          tabsetPanel(
            tabPanel("Counts Table", tableOutput("table2")),
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

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white',
             xlab = 'Waiting time to next eruption (in mins)',
             main = 'Histogram of waiting times')
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
