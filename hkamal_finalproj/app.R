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
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(gplots)


# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Haniya Kamal BF591 Final Project"),
  theme = bslib::bs_theme(bootswatch = "pulse", version = 5),
  sidebarLayout(
    sidebarPanel(
      fileInput('csv_file', label = paste0('Upload your metadata matrix.'), accept = '.csv'),
      fileInput('csv_file2', label = paste0('Upload your counts matrix.'), accept = '.csv'),
      fileInput('csv_file3', label = paste0('Upload your differential expression results.'), accept = '.csv' )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Samples",
          tabsetPanel(
            tabPanel("Summary", paste0("Summary table of Dataset:"), tableOutput("summary")), #Tab with a summary of the table that includes a summary of the type and values in each column
            tabPanel("Table", DTOutput("table")), #Tab with a data table displaying the sample information, with sortable columns
            tabPanel("Plots", #Tab with histograms, density plots, or violin plots of continuous variable3
              radioButtons("sample_plots", "Choose a plot that you would like to see.",
                           choices = c('Diagnosis vs. Age of Death Violin Plot', 
                                       'RNA Integrity Number (RIN) Histogram', 'PMI Histogram')),
              plotOutput("plot1")
            )
          )
        ),
        tabPanel("Counts",
          tabsetPanel(
            sliderInput("variance", "Slider to include genes with at least X percentile of variance", min = 0, max = 100, value = 50),
            sliderInput("num_genes", "Slider to include genes with at least X samples that are non-zero", min = 0, max = 69, value = 69),
            tabPanel("Counts Summary Statistics",
                     verbatimTextOutput("summary_stats")),
            tabPanel("Scatterplots",
                     plotOutput("combined_scatter")),
            tabPanel("Heatmap",
              plotOutput("heatmap")  # Clustered heatmap of counts remaining after filtering
                   ),
            tabPanel("PCA plot", numericInput('x_axis', 'Choose a Principle Component as the X axis:', value = 1, min = 1, max = 10),
                     numericInput('y_axis', 'Choose a Principle Component as the Y axis:', value = 2, min = 1, max = 10),
                     plotOutput("pcaplot")
            )
          )
        ),
        tabPanel("Differential Expression",
            tabsetPanel(
              tabPanel("Differential Expression Results Table", DTOutput("deseqtable")), #sortable table displaying differential expression results
              tabPanel("Differential Expression Plots", 
                radioButtons("deseq_plot_button", "Choose an analysis that you would like to see.",
                             choices = c('Histogram of raw pvalues', 'Histogram of Log2FoldChanges', 'LogFC Volcano Plot')),
                plotOutput("deseq_plots"))
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
  options(shiny.maxRequestSize = 30*1024^2) # so program can take a larger file size
  
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
    d2 <- read.csv(data2$datapath, header = TRUE, row.names = 1)
    #print (d2)
    return(d2)
  })
  
  load_deseq_data <- reactive({ #CSV file can be uploaded to show the deseq results
    req(input$csv_file3)
    deseq_data <- input$csv_file3
    if (is.null(deseq_data)) {
      return(NULL)
    }
    deseq <- read.csv(deseq_data$datapath, header = TRUE)
    return(deseq)
  })
   
  labeled_deseq <- reactive({             #creating a dataframe with a column to indicate if the gene is 
    deseq_results <- load_deseq_data()    #up or down regulated, or neither
    
    lab_df <- deseq_results %>%
      rownames_to_column(var = "Gene") %>%
      mutate(
        log2FC = log2FoldChange,
        padjust = (padj),
        volcplot_status = case_when(
          padjust < 0.01 & log2FC > 0 ~ "UP",
          padjust < 0.01 & log2FC < 0 ~ "DOWN",
          TRUE ~ "NS"
        )
      ) %>%
      select(Gene, log2FC, padjust, volcplot_status)
    #print(lab_df)
    return(lab_df)
    
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
  
  filteredData <- reactive({
    req(load_counts())
    variance_threshold <- quantile(apply(load_counts(), 1, var), probs = input$variance / 100)
    num_genes_threshold <- input$num_genes
    
    
    filtered_genes_var <- rownames(load_counts())[apply(load_counts(), 1, var) >= variance_threshold]
    filtered_genes_non_zero <- rownames(load_counts())[rowSums(load_counts() > 0) >= num_genes_threshold]
    #filtered_genes <- intersect(filtered_genes_var, filtered_genes_non_zero)
    
    
    #f <- filter(load_counts(),apply(load_counts(), 1, var) >= variance_threshold, rowSums(load_counts() > 0) >= num_genes_threshold)
    #  filtered <- load_counts()[filtered_genes, , drop = FALSE]
    
    f <- load_counts() %>% 
      mutate(Passes_Filter = ifelse(rowSums(. > 0) >= num_genes_threshold & apply(., 1, var) >= variance_threshold, "Yes", "No"))
    return(f)
  })
  
  pca_data <- function(dataf){
    #counts <- load_counts()
    #metadata <- load_data()
    
    numeric_data <- dataf[, sapply(dataf, is.numeric)]
    expr_mat <- as.matrix(numeric_data) #pca plot requires numeric data
    expr_mat_centered <- scale(expr_mat)
    pca <- prcomp(expr_mat_centered, center = FALSE, scale = TRUE)
    return(pca)
  }
  
  
  deseq_func <- function(data) { #this is the function that the output will pull from to show deseq data
    deseq_df <- data
    return(deseq_df)
  }
  
  
  output$summary <- renderTable({ #pulling the summary table from the function, part of samples tabset
    summary_df <- load_data()
    summarytable(summary_df)
    
  })
  
  output$table <- renderDT({ #pulling the table from the Samples tabset
    data2 <- load_data()             #table will be sortable
    datatable1(data2) #loading the sample table 
  })
  
  output$plot1 <- renderPlot({ #plots for the sample tab
    plot_type <- input$sample_plots #setting the radio button 
    sample_tab_data <- load_data()
    
    if (plot_type == 'Diagnosis vs. Age of Death Violin Plot') {
      plot <- ggplot(sample_tab_data) +
        geom_violin(aes(x=!!sym('Diagnosis'), y=!!sym('Age.of.Death'), fill=Diagnosis))
      return(plot)
    } else if (plot_type == 'RNA Integrity Number (RIN) Histogram') {
      plot2 <- ggplot(sample_tab_data) +
        geom_histogram(aes(x=!!sym('Rin'), fill=Diagnosis))
      return(plot2)
    } else if (plot_type == 'PMI Histogram') {
      plot3 <- ggplot(sample_tab_data) +
        geom_histogram(aes(x=!!sym('PMI'), fill=Diagnosis))
      return(plot3)
    }
  })
  
  
  output$countstable <- renderTable({ #show the counts data in the counts matrix tabset
    counts <- load_counts()
    
    if (!is.null(counts)) {
      non_zero_counts <- apply(counts[, -1], 1, function(x) sum(x > 0))  # Count non-zero values for each gene
      #filter genes based on the slider input
      selected_genes <- counts[non_zero_counts >= input$num_genes, ]
      
      # Use the selected genes to filter the original counts data
      counts <- counts[counts$gene %in% selected_genes$gene, ]
    }
    
    return(counts)
  })
  
  
  output$summary_stats <- renderText({
    counts_data <- load_counts()
    filtered_data <- filteredData()
    
    #print(filtered_data)
    paste0("Number of samples: ", ncol(counts_data),
           "\nTotal number of genes: ", nrow(counts_data),
           "\nNumber of genes passing current filter: ", sum(filtered_data$Passes_Filter == "Yes"),
           "\nPercent of genes passing current filter: ", sum(filtered_data$Passes_Filter == "Yes") / nrow(load_counts()) * 100, "%",
           "\nNumber of genes not passing current filter: ", sum(filtered_data$Passes_Filter == "No"),
           "\nPercent of genes not passing current filter: ", sum(filtered_data$Passes_Filter == "No") / nrow(load_counts()) * 100, "%")
  })
  
  output$combined_scatter <- renderPlot({
    filtered_data <- filteredData()
    numeric_cols <- Filter(is.numeric, filtered_data) #excluding the gene and passes_filter column
    
    # test <- apply(numeric_cols, 1, function(x) median(x, na.rm = TRUE))
    stats_data <- numeric_cols %>% 
      mutate(Medians = apply(numeric_cols, 1, function(x) median(x, na.rm = TRUE)),
             variances = apply(numeric_cols, 1, var, na.rm = TRUE),
             num_zeros = apply(numeric_cols, 1, function(x) sum(x==0, na.rm = TRUE)))
    
    stats_data$Passes_Filter <- filtered_data$Passes_Filter
    
    
    # median count vs variance
    plotA <- ggplot(stats_data, aes(x = log(Medians),
                                    y = log(variances),
                                    color = Passes_Filter)) +
      geom_point() +
      labs(title = "Median Count vs Variance",
           x = "Log(Median Count)",
           y = "Log(Variance)") +
      scale_color_manual(values = c("Yes" = "darkblue", "No" = "red"))
    
    # median count vs number of zeros
    plotB <- ggplot(stats_data, aes(x = Medians, y= num_zeros, color = Passes_Filter)) +
      geom_point() +
      labs(title = "Median Count vs Number of Zeros",
           x = "Median Count",
           y = "Number of Zeros") +
      scale_color_manual(values = c("Yes" = "darkblue", "No" = "red"))
    
    two_plots <- cowplot::plot_grid(plotA, plotB, ncol = 1, align = 'v')                
    
    return(two_plots)
  })
  
  output$heatmap <- renderPlot({
    filtered_data <- filteredData()
    numeric_cols <- Filter(is.numeric, filtered_data)
    
    matrix <- as.matrix(numeric_cols)
    color_palette <- colorRampPalette(c("yellow2", "goldenrod", "darkred"))(50)
    
    # normalizing rows
    matrix <- scale(matrix, center = TRUE, scale = TRUE)
    
    map <- reactive({
      heatmap.2(matrix, col = color_palette, dendrogram = 'none', trace = 'none', breaks = seq(-10, 10, length = 51))
    })()
    
    return(map)
  })
  
  output$pcaplot <- renderPlot({
    pca <- pca_data(load_counts()) 
    
    x_pc <- input$x_axis 
    y_pc <- input$y_axis 
    pc_data <- as.data.frame(pca$x[, c(x_pc, y_pc)])  
    
    percent_variance <- (pca$sdev^2) / sum(pca$sdev^2) * 100  
    percent_variance_x <- round(percent_variance[x_pc], 2)  
    percent_variance_y <- round(percent_variance[y_pc], 2)  
    
    # Creating plot
    ggplot(pc_data, aes(x = pc_data[, 1], y = pc_data[, 2])) +
      geom_point() +
      labs(x = paste0("PC", x_pc, " (", percent_variance_x , "%)"),
           y = paste0("PC", y_pc, " (", percent_variance_y , "%)"),
           title = "PCA Plot")
  })
  
  
  output$deseqtable <- renderDT({ #show the sortable Deseq table
    deseq_df_results <- load_deseq_data()
    deseq_func(deseq_df_results)
  })
  
 
  output$deseq_plots <- renderPlot({ 
    deseq_plot_type <- input$deseq_plot_button
    deseq_tab_data <- load_deseq_data()
    #input$button
    #isolate
    if (deseq_plot_type == 'Histogram of raw pvalues') { #a histogram of the raw p-values from the DESeq2 results
      his_plot1 <- ggplot(deseq_tab_data) +
        geom_histogram(mapping=aes(x=!!sym('pvalue')), bins = 50,  fill = 'lightblue', color = 'black') +
        ggtitle('Histogram of raw pvalues obtained from DE analysis')
      return(his_plot1)
    } else if (deseq_plot_type == 'Histogram of Log2FoldChanges') {   #log2foldchange from DESeq2 results in a histogram
      his_plot2 <- ggplot(deseq_tab_data) +
        geom_histogram(mapping = aes(x=!!sym('log2FoldChange')), bins = 100, fill = 'lightblue', color = 'black') +
        ggtitle('Histogram of Log2FoldChanges for DE Genes')
      return(his_plot2)
    } else if (deseq_plot_type == 'LogFC Volcano Plot') {  #volcano plot that displays log2foldchange vs -log10(padj) and labeled by status
      volc <- ggplot(data = labeled_deseq()) + #passing the data dataframe from a function i created above 
        geom_point(aes(x = log2FC, y = -log10(padjust), color = volcplot_status)) +
        ggtitle('Volcano plot of DESeq2 differential expression results')
      return(volc)
    }
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
