###>>>>> WELCOME TO THE R-SHINY APP! <<<<<###

# Install necessary libraries
library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(reshape2)
library(stats)
library(ggbeeswarm)
library(gridExtra)
library(viridis)




# UI
ui <- fluidPage( #creating our fluid layout
  titlePanel("RNA-seq Data Analysis"), #title of the app
  
  # Element 1: Metadata analysis
  tabsetPanel(
    tabPanel("Metadata Analysis",
             sidebarLayout(
               sidebarPanel( # the following elements will be in the sidebar panel 
                 h3("Summary and Visualization"), #creating a level 3 heading
                 p("This tab gives you a summary of your data and allows you to visualize density plots for different characteristics of your samples"),
                 # a short description about the characteristics of the tab  
                 fileInput("metadata_file_upload", "Upload Metadata CSV File", accept = ".csv") # let's the user upload metadata
               ),
               
               #components of this tab
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary", tableOutput("metadata_summary_table")), #first sub-tab, summarizes the metadata
                   tabPanel("Data Table", dataTableOutput("metadata_data_table")), #second sub-tab displays all the complete metadata
                   tabPanel(
                     "Plots", #third sub-tab provides density-plot depending on the below input values
                     selectInput("metadata_x_col", "Select Column to Plot:", choices = NULL), #these values will be plotted
                     selectInput("metadata_group_col", "Select Grouping Column:", choices = NULL),# the values being plotted will be grouped according to this 
                     plotOutput("metadata_plot")
                   )
                 )
               )
             )
    ),
    
    # Element 2: Counts Matrix Analysis
    tabPanel("Counts Matrix Analysis", 
             sidebarLayout(
               sidebarPanel(
                 #the sidebar contains the following selections:
                 
                 #let the user upload their counts file 
                 fileInput("counts_file_upload", "Upload Normalized Counts Matrix (CSV)", accept = ".csv"), 
                 
                 #heading for clarity
                 h3("Exploring Counts Matrix"),
                 
                 #tab description
                 p("This tab allows filtering genes based on variance percentile. You can also visualise your filtered genes in the form of diagnostic plots
                   heatmap and principal components"),
                 
                 #this slider with allow the user to set a variance percentile threshold, the default is 50
                 sliderInput("variance_threshold", "Variance Percentile Threshold:", 
                             min = 0, max = 100, value = 50, step = 1),
                 
                 #this slider should allow the user to set a minimum of number non-zero samples for a gene to be filtered
                 #Default values for this slider are adjusted dynamically, meaning if the sample has 72 samples, the default will set to 36
                 sliderInput("nonzero_threshold", "Minimum Non-Zero Samples:", 
                             min = 0, max = 100, value = 1, step = 1), 
                 
                 #PCA radiobuttons are adjusted in the sidebar itself to avoid cluttering of the plot in the PCA sub-tab
                 radioButtons("pca_type", "PCA Visualization Type:",
                              choices = c("Scatter Plot (Custom PCs)" = "scatter", #user can display either of the PCA plots
                                          "Beeswarm of Top PCs" = "beeswarm")),
                 
                 conditionalPanel(
                   # this will allow the user to change the principal component being displayed on either of the axes
                   condition = "input.pca_type == 'scatter'",
                   fluidRow(
                     column(6, uiOutput("pc_x_selector")),
                     column(6, uiOutput("pc_y_selector"))
                   )
                 ),
                 
                 conditionalPanel(
                   # here the user can choose between the number of PCs they want to make a beeswarm plot for
                   condition = "input.pca_type == 'beeswarm'",
                   numericInput("top_n_pcs", "Number of Top Principal Components:", 
                                value = 5, min = 2, max = 10, step = 1)
                 )
               ),
               
               #Components of this tab:
               mainPanel(
                 tabsetPanel(
                   
                   #first sub-tab displays a summary of the genes being filtered out according to the set threshold of variance and non-zero samples
                   tabPanel("Filter Summary", 
                            verticalLayout(
                              tableOutput("filter_summary"),
                              textOutput("filter_description")
                            )),
                   
                   #second tab displays the diagnostic plots that can give a correlation between the variance/number of zero samples and the median
                   tabPanel("Diagnostic Plots", 
                            plotOutput("diagnostic_plots", height = "600px")),
                   
                   #the third tab provides a clustered heatmap
                   tabPanel("Clustered Heatmap", 
                            plotOutput("heatmap", height = "800px")),
                   
                   #the fourth tab provides the PCA scatterplot and beeswarm plot
                   tabPanel("PCA Visualization", 
                            plotOutput("pca_plot", height = "600px"))
                 )
               )
             )
    ),
    
    #Element 3: Differential Gene Expression
    tabPanel("Differential Expression Analysis",
             sidebarLayout(
               sidebarPanel(
                 # components of side panel:
                 
                 # allows user to upload the DESeq results file
                 fileInput("deseqFile", "Upload DESeq Results CSV", accept = ".csv"),
                 
                 #heading of the tab
                 h3("Differential Gene Expression"),
                 
                 #description
                 p("This tab lets you visualize your differential gene expression data"),
                 
                 #to change the stringency for evaluation of differentially expressed genes, two sliders have been included:
                 #slider 1 allows altering the log2foldchange threshold, threshold is 1.5
                 sliderInput("log2fcThreshold", "Log2 Fold Change Threshold:",
                             min = 0, max = 5, value = 1.5, step = 0.1),
                 
                 #slider 2 allows altering the threshold for the adjusted p-value, threshold is 0.05 
                 sliderInput("padjThreshold", "Adjusted p-value Threshold:",
                             min = 0, max = 0.1, value = 0.05, step = 0.01)
               ),
               mainPanel(
                 tabsetPanel(
                   
                   #first tab displays the results in the uploaded file in a table format 
                   tabPanel("Table View", 
                            DTOutput("resultsTable")
                   ),
                   
                   #second tab holds all the four plots for the representation of the DEGs
                   tabPanel("Plots", 
                            fluidRow(
                              column(6, plotOutput("pvalHist")),
                              column(6, plotOutput("log2fcHist"))
                            ),
                            plotOutput("volcanoPlot"),
                            plotOutput("topGenesPlot")
                   )
                 )
               )
             )
    ),
    
    # Element 4: Gene set Enrichment Analysis (with FGSEA)
    
    tabPanel("FGSEA Analysis Viewer",
             sidebarLayout(
               sidebarPanel(
                 #sidebar components:
                 #user can upload their FGSEA results here
                 fileInput("upload_fgsea", 
                           "Upload FGSEA Results CSV File",
                           accept = c(".csv", ".tsv", ".txt")),
                 #heading
                 h3("Gene Enrichment Analysis"),
                 #description of the tab
                 p("This tab lets you visualize your top enriched pathways, differentiate between Positive and Negative NES pathways and 
               download and visualize this filtered data"),
                 uiOutput("pathway_controls")
               ),
               
               mainPanel(
                 tabsetPanel(
                   # First sub-tab: Top Pathways Barplot
                   tabPanel("Top Pathways",
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("top_pathways", 
                                            "Number of Top Pathways:", 
                                            min = 5, max = 50, 
                                            value = 20)
                              ),
                              mainPanel(
                                plotOutput("fgsea_barplot")
                              )
                            )
                   ),
                   # Second sub-tab: Detailed Results Table
                   tabPanel("Results Table",
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("p_adj_filter", 
                                            "Adjusted P-value Threshold:", 
                                            min = 0, max = 1, 
                                            value = 0.25, step = 0.05),
                                radioButtons("nes_filter", 
                                             "NES Filter:",
                                             choices = c("All", "Positive", "Negative"),
                                             selected = "All"),
                                downloadButton("download_table", "Download Results")
                              ),
                              mainPanel(
                                DTOutput("results_table")
                              )
                            )
                   ),
                   # Third sub-tab: NES vs P-value Scatter Plot
                   tabPanel("NES Scatter Plot",
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("scatter_p_adj", 
                                            "Adjusted P-value Threshold:", 
                                            min = 0, max = 1, 
                                            value = 0.25, step = 0.05)
                              ),
                              mainPanel(
                                plotOutput("nes_scatter")
                              )
                            )
                   )
                 )
               )
             )
    )
  )

)




# SERVER

server <- function(input, output, session) {
  

  
   
 ##### METADATA ANALYSIS SERVER #####
  
  # Increased file upload size to 15MB since counts file tend to be computationally heavy
  options(shiny.maxRequestSize = 15*1024^2)
  
  # Reactive expression to load the metadata file
  metadata <- reactive({
    req(input$metadata_file_upload)
    
    # then we read CSV and convert empty strings to NA to make these values are excluded before plotting
    read.csv(input$metadata_file_upload$datapath, stringsAsFactors = FALSE, 
             na.strings = c("", "NA", "N/A", " "))
  })
  
  #  then we extract numeric values from columns while again handling NA, maximum robustness of the code is appreciated so...
  extract_numeric <- function(column) {
    
    column <- column[!is.na(column)] # removing NA
    
    # extracting numeric values from the metadata file to be able to use them for computation or plotting
    # this way the column containing sample type or condition 
    numeric_values <- as.numeric(gsub("[^0-9.]+", "", column))
    
    # Remove any resulting NA values
    numeric_values <- numeric_values[!is.na(numeric_values)]
    
    return(numeric_values)
  }
  # Update metadata column choices dynamically
  observe({
    req(metadata())
    
    # Get all column names
    all_columns <- names(metadata())
    
    # Removing the first column (Sample Condition) from plot column choices, this was to ensure that it doesn't cause any issues in the code during plotting
    plot_columns <- all_columns[-1]
    
    # Updating plot column dropdown without the first column
    # this section specifically provides the columns that can be plotted (numerical values)
    updateSelectInput(session, "metadata_x_col", choices = plot_columns)
    
    # Updating group column dropdown with all columns
    # This section can be used to group the plotted values
    updateSelectInput(session, "metadata_group_col", choices = all_columns)
  })
  
  # Generating summary table for metadata
  output$metadata_summary_table <- renderTable({
    req(metadata())
    # the lapply function applies the below functions to the data, column-by-column
    summary <- lapply(metadata(), function(column) {
     
      column <- column[!is.na(column)] # more NA handling
      
      numeric_values <- extract_numeric(column)
      if (length(numeric_values) > 0) {
        mean_val <- mean(numeric_values, na.rm = TRUE) #calculating the mean
        sd_val <- sd(numeric_values, na.rm = TRUE) # checking for all distinct values 
        c("Type" = "integer", "Mean/Distinct value" = sprintf("%.2f (+/- %.2f)", mean_val, sd_val)) #this part returns a summary that contains the type of data in a column
      } else { #in case of non-numerical values, we identify all distinct or unique values in the table 
        distinct_vals <- unique(column)
        c("Type" = "character", "Mean/Distinct value" = paste(distinct_vals, collapse = ", "))
      }
    })
    summary_df <- do.call(rbind, summary) %>% as.data.frame() #Combines the summaries for all columns (stored as lists) into a single data frame
    summary_df <- cbind(Column = names(metadata()), summary_df) #Adds a new column with the names of the original columns
    summary_df
  })
  
  # Rendering data table for metadata 
  output$metadata_data_table <- renderDataTable({
    req(metadata())
    datatable(metadata())
  })
  
  # Then we generate density plots for different columns of the metadata
  output$metadata_plot <- renderPlot({
    req(metadata(), input$metadata_x_col, input$metadata_group_col)
    
    # Preparing plot data
    plot_data <- metadata() %>% 
      # Removing rows where either x_col (user selected column for plotting) or group_col(user selected column for grouping) is NA
      filter(!is.na(.[[input$metadata_x_col]]) & !is.na(.[[input$metadata_group_col]]))
    
    # Extracting numeric values for the selected x column, these will be plotted
    numeric_values <- extract_numeric(plot_data[[input$metadata_x_col]])
    
    # Preparing data for plotting
    plot_data <- data.frame(
      x = numeric_values,
      group = plot_data[1:length(numeric_values), input$metadata_group_col] #this extracts corresponding values from the grouping column
    )
    
    # Creating the density plot
    ggplot(plot_data, aes(x = x, fill = group)) +
      geom_density(alpha = 0.5) +
      theme_minimal() +
      labs(
        title = paste("Density Plot of", input$metadata_x_col, "by", input$metadata_group_col), 
        x = input$metadata_x_col,
        fill = input$metadata_group_col
      ) +
      scale_fill_manual(values = c("blue", "red"))  
  })
  
 
  
  
   ##### COUNTS MATRIX ANALYSIS SERVER #####
  
  # Reactive expression to load counts data
  counts_data <- reactive({
    req(input$counts_file_upload)
    read.csv(input$counts_file_upload$datapath, row.names = 1)
  })
  
  # Dynamic slider adjustment so that the slider scale can change according to the total number of samples present
  observe({
    req(counts_data())
    updateSliderInput(session, "nonzero_threshold",
                      max = ncol(counts_data()),
                      value = max(1, round(ncol(counts_data()) * 0.5))) #the default threshold is set to half of the total number of samples
  })
  # Then we have to generate our filtered data (Reactive expression for filtered data)
  filtered_data <- reactive({
    req(counts_data())
    counts <- counts_data()
    
    # Then we calculate variances and non-zero counts
    variances <- apply(counts, 1, var) # computes variance for each row (gene) across samples
    variance_cutoff <- quantile(variances, probs = input$variance_threshold / 100) #determines the threshold variance, based on the user-specified percentile
    
    nonzero_counts <- rowSums(counts > 0) # how many samples each gene is present in (non-zero values)
    nonzero_cutoff <- input$nonzero_threshold #user-selected minimum count threshold from the slider
    
    # Filter genes
    # we store genes that meet both variance and non-zero thresholds
    passing_genes <- (variances >= variance_cutoff) & (nonzero_counts >= nonzero_cutoff)
    filtered_counts <- counts[passing_genes, ] # then subset the counts dataframe to include only filtered genes
    
    list(
      filtered_counts = filtered_counts,
      total_genes = nrow(counts),
      filtered_genes = nrow(filtered_counts),
      passing_genes = rownames(filtered_counts)
    )
  })
  
  # Filtering summary table
  output$filter_summary <- renderTable({
    req(counts_data())
    filtered_info <- filtered_data() #extracts filtered data summary
    
    # Calculating additional metrics for displaying in the summary
    total_genes <- filtered_info$total_genes
    filtered_genes <- filtered_info$filtered_genes
    not_filtered_genes <- total_genes - filtered_genes #direct difference gives us the number of non-filtered genes
    not_filtered_percentage <- sprintf("%.2f%%", (not_filtered_genes / total_genes) * 100) #computes the percentage of genes that did not pass filtering
    
    # Creating the summary table
    data.frame(
      Metric = c(
        "Total Samples", 
        "Total Genes", 
        "Filtered Genes", 
        "Filtered Percentage",
        "Genes Not Passing Filter",
        "Not Passing Percentage"
      ),
      Value = c(
        ncol(counts_data()), 
        total_genes, 
        filtered_genes, 
        sprintf("%.2f%%", (filtered_genes / total_genes) * 100),
        not_filtered_genes,
        not_filtered_percentage
      )
    )
  })
  
  # printing the filter description for the summary table
  output$filter_description <- renderText({
    req(filtered_data())
    filtered_info <- filtered_data()
    
    total_genes <- filtered_info$total_genes
    filtered_genes <- filtered_info$filtered_genes
    not_filtered_genes <- total_genes - filtered_genes
    
    paste0(
      "Filtering genes with variance above the ", input$variance_threshold, 
      "th percentile and present in at least ", input$nonzero_threshold, 
      " samples. ",
      filtered_genes, " genes (", sprintf("%.2f%%", (filtered_genes / total_genes) * 100), 
      ") passed the filter, while ", not_filtered_genes, 
      " genes (", sprintf("%.2f%%", (not_filtered_genes / total_genes) * 100), 
      ") did not."
    )
  })
  
  # Then we move to the diagnostic plots
  output$diagnostic_plots <- renderPlot({
    req(counts_data())
    counts <- counts_data()
    filtered_info <- filtered_data()
    
    # Computing metrics
    variances <- apply(counts, 1, var)
    nonzero_counts <- rowSums(counts > 0)
    medians <- apply(counts, 1, median) #here, we calculate the median expression level for each gene
    
    # Preparing plot data
    plot_data <- data.frame(
      Median = medians,
      Variance = variances,
      Zero = ncol(counts) - nonzero_counts,
      Filtered = rownames(counts) %in% filtered_info$passing_genes #thi is basically a boolean that indicates whether a gene has passed the filter or not
    )
    
    # Plot 1: Median vs Variance plot (log scale)
    # we color code the points according to their filtered status to depict how many genes pass at a given threshold of variance
    p1 <- ggplot(plot_data, aes(x = Median, y = Variance, color = Filtered)) +
      geom_point(alpha = 0.7, size = 2) +
      scale_color_manual(values = c("TRUE" = "#1E90FF", "FALSE" = "lightgrey")) +
      scale_y_log10(labels = scales::scientific) +
      scale_x_log10(labels = scales::scientific) +
      theme_minimal() +
      labs(
        title = "Median Count vs Variance",
        x = "Median Count",
        y = "Variance (log10)",
        color = "Passes Filter"
      )
    
    # Plot 2: Median vs Number of Zero Samples plot
    p2 <- ggplot(plot_data, aes(x = log10(Median + 1), y = Zero, color = Filtered)) +
      geom_point(alpha = 0.7, size = 2) +
      scale_color_manual(values = c("TRUE" = "#1E90FF", "FALSE" = "lightgrey")) +
      theme_minimal() +
      labs(
        title = "Median Count vs Number of Zero Samples",
        x = "Median Count (log10)",
        y = "Number of Zero Samples",
        color = "Passes Filter"
      )
    
    # Arranging plots side by side for aesthetics
    gridExtra::grid.arrange(p1, p2, ncol = 2)
  })
  
  # Here I am just using a tiny notification because my heatmap takes a substantial amount of time to load!
  observeEvent(filtered_data(), {
    showNotification(
      "Please be patient, Miss Heatmap takes her own time to doll up!",
      duration = NULL,  # Notification stays till prompted to disappear
      id = "heatmap_notify"
    )
  })
  
  # Clustered Heatmap
  output$heatmap <- renderPlot({
    req(filtered_data())
    counts <- filtered_data()$filtered_counts
    
    # Here we apply log2 transformation with a pseudo-count of 1 to avoid log(0)
    counts_log2 <- log2(counts + 1)
    
    # Then we normalize counts by row (gene-wise z-score)
    # this scales each row (gene) to have mean 0 and standard deviation 1
    counts_scaled <- t(scale(t(counts_log2)))
    
    # Handle any NaN values
    counts_scaled[is.na(counts_scaled)] <- 0
    
    # Render the heatmap
    pheatmap(counts_scaled, 
             cluster_rows = TRUE, 
             cluster_cols = FALSE, 
             color = viridis(50),
             main = "Clustered Heatmap of Filtered Genes (Log2 Normalized)",
             fontsize = 8,
             show_rownames = FALSE)
  })
  
  # PC X-axis selector UI
  # this dynamically generates a dropdown menu for selecting the X-axis PCA component
  # Limits options to the first 10 principal components or the total number of samples
  output$pc_x_selector <- renderUI({
    req(counts_data())
    n_pcs <- min(10, ncol(data()))
    selectInput("selected_pc_x", "X-axis Principal Component:", 
                choices = paste0("PC", 1:n_pcs), 
                selected = "PC1")
  })
  
  # PC Y-axis selector UI
  output$pc_y_selector <- renderUI({
    req(counts_data())
    n_pcs <- min(10, ncol(data()))
    selectInput("selected_pc_y", "Y-axis Principal Component:", 
                choices = paste0("PC", 1:n_pcs), 
                selected = "PC2")
  })
  
  # PCA plot with improved visualization options
  output$pca_plot <- renderPlot({
    req(filtered_data())
    counts <- log2(filtered_data()$filtered_counts + 1) #log transforming the filtered data
    
    # Perform PCA
    pca <- prcomp(t(counts), center = TRUE, scale. = TRUE) # computes principal computes
    
    # Calculate variance explained
    variance_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100 # proportion of variance explained by each PC
    
    # Beeswarm plot for top PCs
    if (input$pca_type == "beeswarm") {
      top_pcs <- as.data.frame(pca$x[, 1:input$top_n_pcs]) # selects the top N principal components based on the user input
      top_pcs$Sample <- colnames(counts)
      melted_pcs <- melt(top_pcs, id.vars = "Sample") # converts the top_pcs data from a wide format to a long format 
      
      ggplot(melted_pcs, aes(x = variable, y = value, color = variable)) +
        geom_beeswarm(alpha = 0.7) +
        theme_minimal() +
        scale_color_viridis_d() +
        labs(
          title = "Beeswarm Plot of Top Principal Components",
          x = "Principal Component",
          y = "Projection Value",
          color = "PC"
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    } 
    
    # Scatter plot for selected PCs
    else {
      # these lines extract the numeric indices of the PCs selected by the user for plotting
      pc_x_num <- as.numeric(sub("PC", "", input$selected_pc_x))
      pc_y_num <- as.numeric(sub("PC", "", input$selected_pc_y))
      
      pca_data <- data.frame(
        PCX = pca$x[, pc_x_num],
        PCY = pca$x[, pc_y_num],
        Sample = colnames(counts)
      )
      
      ggplot(pca_data, aes(x = PCX, y = PCY, label = Sample)) +
        geom_point(size = 3, alpha = 0.7, color = "#1E90FF") +
        geom_text(vjust = -1.5, size = 3) +
        theme_minimal() +
        labs(
          title = "PCA Scatter Plot",
          x = sprintf("%s (%.2f%% Variance)", # this dynamically inserts values into the label strings
                      input$selected_pc_x, 
                      variance_explained[pc_x_num]),
          y = sprintf("%s (%.2f%% Variance)", 
                      input$selected_pc_y, 
                      variance_explained[pc_y_num])
        )
    }
  })
  
  
  
  
  
  ####DESEQ ANALYSIS SERVER####
  
  deseq_data <- reactive({
    req(input$deseqFile)
    deseq <- read.csv(input$deseqFile$datapath, stringsAsFactors = FALSE)
    deseq <- deseq[complete.cases(deseq), ] # another function that removes any rows containing NA values to ensure clean data
    return(deseq)
  })
  
  output$resultsTable <- renderDT({
    req(deseq_data())
    datatable(deseq_data(), options = list(pageLength = 100, search = list(search = "")), rownames = FALSE) #displays 10 rows in a page for better readability
  })
  
  # First we make a p-value vs frequency histogram 
  output$pvalHist <- renderPlot({
    req(deseq_data())
    ggplot(deseq_data(), aes(x = pvalue)) +
      geom_histogram(binwidth = 0.01, fill = "maroon", color = "black") +
      theme_minimal() +
      labs(title = "Histogram of P-values", x = "P-value", y = "Frequency")
  })
  
  # Then we make a lof2foldchange vs frequency histogram
  output$log2fcHist <- renderPlot({
    req(deseq_data())
    filtered <- deseq_data()[deseq_data()$padj <= input$padjThreshold, ]
    ggplot(filtered, aes(x = log2FoldChange)) +
      geom_histogram(binwidth = 0.1, fill = "darkblue", color = "black") +
      theme_minimal() +
      labs(title = "Histogram of Log2 Fold Change", x = "Log2 Fold Change", y = "Frequency")
  })
  
  # Then we make a volcano plot to differentiate the genes into upregulated, downregulated and non-significant categories
  output$volcanoPlot <- renderPlot({
    req(deseq_data())
    volcano <- deseq_data()
    volcano$category <- ifelse(volcano$padj <= input$padjThreshold & volcano$log2FoldChange >= input$log2fcThreshold, "Upregulated",
                               ifelse(volcano$padj <= input$padjThreshold & volcano$log2FoldChange <= -input$log2fcThreshold, "Downregulated", "Not Significant"))
    ggplot(volcano, aes(x = log2FoldChange, y = -log10(padj), color = category)) +
      geom_point() +
      scale_color_manual(values = c("Upregulated" = "darkred", "Downregulated" = "darkgreen", "Not Significant" = "darkgray")) +
      theme_minimal() +
      labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10 Adjusted P-value")
  })
  
  # here we select the top 10 genes with lowest p-adjusted values, these are the top 10 differentially expressed genes
  # for each gene, we plot the normalized count value of that gene in control sample and diseased sample
  output$topGenesPlot <- renderPlot({
    req(deseq_data())
    topGenesData <- deseq_data()[order(deseq_data()$padj), ][1:10, ]
    topGenes <- topGenesData$Gene.ID
    plotData <- data.frame(
      Gene = rep(topGenes, each = 2),
      SampleType = rep(c("Control", "Huntington's Disease"), times = 10),
      NormCount = c(topGenesData$Control.mean, topGenesData$HD.mean)
    )
    
    ggplot(plotData, aes(x = Gene, y = log10(NormCount + 1), color = SampleType)) +
      geom_point(size = 3) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(title = "Plot of Log10(normalized counts) for top ten DE genes", 
           x = "Gene", y = "log10(Normalized Counts)", color = "Sample Type")
  })
  

  
  
  
  
  ###FGSEA ANALYSIS###
  
  fgsea_data <- reactive({
    req(input$upload_fgsea)
    
    # Read the uploaded file
    infile <- input$upload_fgsea
    
    # Read the CSV file
    data <- read.csv(infile$datapath, stringsAsFactors = FALSE)
    
    return(data)
  })
  
  # Reactive data filtering
  filtered_data_fgsea <- reactive({
    req(fgsea_data())
    data <- fgsea_data()
    
    # P-value filtering
    data <- data[data$padj <= input$p_adj_filter, ] #filters rows where the padj is less than or equal to the user-specified threshold 
    
    # NES filtering - we further filter the data
    if (input$nes_filter == "Positive") {
      data <- data[data$NES > 0, ]
    } else if (input$nes_filter == "Negative") {
      data <- data[data$NES < 0, ]
    }
    
    data
  })
  
  # Tab 1: Top Pathways Barplot
  output$fgsea_barplot <- renderPlot({
    req(filtered_data_fgsea())
    top_paths <- head(filtered_data_fgsea() %>% 
                        arrange(padj), 
                      input$top_pathways)
    
    ggplot(top_paths, aes(x = reorder(pathway, NES), y = NES, fill = NES)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      theme_minimal() +
      labs(title = "Top Pathways by Normalized Enrichment Score",
           x = "Pathway", y = "Normalized Enrichment Score")
  })
  
  # Tab 2: Detailed Results Table
  output$results_table <- renderDT({
    req(filtered_data_fgsea())
    filtered_data_fgsea() %>%
      select(pathway, NES, pval, padj, size, leadingEdge)
  }, options = list(pageLength = 10))
  
  # Download handler for results table
  output$download_table <- downloadHandler(
    filename = "fgsea_filtered_results.csv",
    content = function(file) {
      write.csv(filtered_data_fgsea(), file, row.names = FALSE)
    }
  )
  
  # Tab 3: NES vs P-value Scatter Plot
  output$nes_scatter <- renderPlot({
    req(filtered_data_fgsea())
    data <- fgsea_data()  # Use full dataset for consistent plotting
    
    ggplot(data, aes(x = NES, y = -log10(padj))) +
      geom_point(aes(color = padj <= input$scatter_p_adj), alpha = 0.7) +
      scale_color_manual(values = c("grey", "red"), 
                         name = "Significant",
                         labels = c("Not Significant", "Significant")) +
      scale_x_continuous(limits = c(-1, 3)) +  # Fixed x-axis from -1 to 3
      theme_minimal() +
      labs(title = "NES vs Adjusted P-value",
           x = "Normalized Enrichment Score", 
           y = "-log10(Adjusted P-value)")
  }, res = 96)  # Added resolution for clearer plotting

  }
  

  # Run the application
  shinyApp(ui = ui, server = server)