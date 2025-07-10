# Load required packages
#remotes::install_version("rstatix", version = "0.7.0")

library(shiny)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(DT)
library(dplyr)
library(readxl)
library(tidyverse)

#library(rstatix)

#metadataColumns<-c("SampleID","Groups","ProteinConc","ReconstitutionVolume","AdjustingValue","DilutionValue","Group2","Group3","comments")
metadataColumns<-c("SampleID","ProteinConc","ReconstitutionVolume","AdjustingValue","DilutionValue","comments")


compute_pvalues <- function(data) {
  # Create an empty data frame to store results
  p_values_full <- data.frame(variable = character(),
                         group1 = character(),
                         group2 = character(),
                         p = numeric(),
                         p.adj = numeric(),
                         stringsAsFactors = FALSE)
  
  # Get unique variables
  variables <- unique(data$variable)
  
  # Loop through each variable
  for (var in variables) {
    # Subset data for the current variable
    subset_data <- data[data$variable == var, ]
    if(data$groupsIncluded[1]=="T"){
      # Create a combined identifier for Groups and CombinedMeta
      subset_data$CombinedGroup <- paste(subset_data$Groups, subset_data$CombinedMeta, sep = "_")
    }else{
      subset_data$CombinedGroup <- subset_data$CombinedMeta
    }
    
    # # Create a combined identifier for Groups and CombinedMeta
    # subset_data$CombinedGroup <- paste(subset_data$Groups, subset_data$CombinedMeta, sep = "_")
    
    # Get counts of elements in each CombinedGroup
    group_counts <- table(subset_data$CombinedGroup)
    
    # Filter for groups with at least 3 elements
    valid_groups <- names(group_counts[group_counts >= 3])
    subset_data <- subset_data[subset_data$CombinedGroup %in% valid_groups, ]
    
    # Skip if less than two valid groups remain
    if (length(valid_groups) < 2) next
    
    # Perform pairwise t-tests for all combinations of valid groups
    group_combinations <- combn(valid_groups, 2, simplify = FALSE)
    p_values <- data.frame(variable = character(),
                                group1 = character(),
                                group2 = character(),
                                p = numeric(),
                                p.adj = numeric(),
                                stringsAsFactors = FALSE)
    for (comb in group_combinations) {
      group1 <- comb[1]
      group2 <- comb[2]
      
      # Extract values for the two groups
      values1 <- subset_data$value[subset_data$CombinedGroup == group1]
      values2 <- subset_data$value[subset_data$CombinedGroup == group2]
      
      
      # Perform t-test with error handling
      test_result <- tryCatch(
        {t.test(values1, values2)},error = function(e) {NULL}
      )
      
      # Skip this iteration if test_result is NULL
      if (is.null(test_result)) {next}
     
      
      # Store p-value result
      p_values <- rbind(p_values, data.frame(
        variable = var,
        group1 = group1,
        group2 = group2,
        p = test_result$p.value,
        stringsAsFactors = FALSE
      ))
    }
    
    # Adjust p-values within the current variable using the Benjamini-Hochberg method
    p_values$p.adj[p_values$variable == var] <- p.adjust(
      p_values$p[p_values$variable == var],
      method = "BH"
    )
    p_values_full <- rbind(p_values_full,p_values)
  }
  
  # Sort results by variable and adjusted p-value
  p_values_full <- p_values_full[order(p_values_full$variable, p_values_full$p.adj), ]
  
  return(p_values_full)
}

# Function to prepare the comparison table
prepare_comparison_table <- function(data) {
  # Filter for significant adjusted p-values
  significant_data <- data[data$p.adj < 0.05, ]
  #significant_data <- data
  
  # Sort by 'variable', 'group1', and 'group2'
  sorted_data <- significant_data[order(
    significant_data$variable,
    significant_data$group1,
    significant_data$group2
  ), ]
  
  return(sorted_data)
}


# Define UI for the application
ui <- fluidPage(
  titlePanel("Standard Curve and Sample Measurements"),
  
  # Tab layout for different sections including the 'How To Use' tab
  tabsetPanel(
    tabPanel(
      "Figure Generation",
      sidebarLayout(
        sidebarPanel(
          fileInput('file1', 'Choose CSV File', accept = c('text/csv', 'text/comma-separated-values,text/plain')),
          
          # UI for selecting a sheet if it's an Excel file
          uiOutput("sheetSelector"),
          
          # Dynamic inputs for each sample column will appear here
          uiOutput("dynamicInputs"),
          
          # A button to confirm the inputs and trigger the next steps
          actionButton("confirmValues", "Confirm Values"),
          uiOutput("checkboxesUI"), # Placeholder for dynamic checkboxes
          # Selection options for group and sample filtering
          h3("Customize Plot"),
          uiOutput("groupSelector"),
          actionButton("applyFilter", "Apply Group Filter"),
          
          # Update data table by removing outlier samples manually
          h3("Remove Selected Sample by ID"),
          uiOutput("sampleRemover"),
          actionButton("applysampleRemover", "Remove Selected Sample"),
          
          # Inputs for PDF filename and dimensions
          h3("Download Options"),
          #textInput("filename", "Enter file name (without extension):", value = "boxplot"),
          #numericInput("pdfWidth", "PDF Width (inches):", value = 8, min = 1, step = 0.5),
          #numericInput("pdfHeight", "PDF Height (inches):", value = 8, min = 1, step = 0.5),
          uiOutput("downloadPlotUI")
        ),
        
        mainPanel(
          h3("Standard Curves"),
          uiOutput("curvePlots"),
          h3("Sample Concentrations"),
          plotOutput("boxPlot"),
          checkboxInput('showDots', 'Show Dots on Boxplot', FALSE),
          checkboxInput('showNames', 'Show SampleIDs', FALSE),
          p("only works for some plot types (Breeder plots)"),
          numericInput("y_max", "Y-axis max value:", value = NA, min = 0, step = 1),
          numericInput("x_min", "X-axis min:", value = 0, min = 0, step = 1),
          numericInput("x_max", "X-axis max:", value = NA, min = 0, step = 1),
          
          sliderInput('dotSize', 'Dot Size', min = 0.25, max = 2, value = 1),
          # Display default min and max values
          textOutput("defaultYRange"),
          # Y-axis customization UI
          uiOutput("yAxisControls"), # Placeholder for dynamic controls
          h3("Data Preview"),
          DT::dataTableOutput("dataPreview"),
          h3("Group Comparisons"),
          p("Shows only comparison values with significant adj P values and at least 3 datapoints"),
          DT::dataTableOutput("pvalueTable")
        )
      )
    ),
    
    # New Tab for 'How To Use' with PNG images and descriptions
    tabPanel(
      "How To Use",
      h3("How to Use This Application"),
      h4("The input file type required is a .CSV file. Below are five examples of different acceptable input formats"),
      fluidRow(
        column(12,
               p("Input 1: Basic input which will plot raw value figures:"),
               img(src = "inputFormat1.png", height = "700px", width = "700px")
        )
      ),
      fluidRow(
        column(12,
               p("Input 2: This input type can be used to view comparisons between different characteristics of a dataset."),
               p("When using this input type, the columns can be in any order, but the blue ones are required."),
               p("If you want to add meta data to select and group data by, add META to the front of that column."),
               p("Remaining columns will be viewed as samples to be plotted."),
               img(src = "inputFormat2_updated.png", height = "550px", width = "1300px")
        )
      )
      ,
      fluidRow(
        column(12,
               p("Input 3: This input type can be used for HPLC or other assays and can calculate concentrations based on standards."),
               p("When using this input type, the application will ask for different values needed."),
               img(src = "inputFormat2.png", height = "900px", width = "900px")
        )
      ),
      
      fluidRow(
        column(12,
               p("Input 4: This input type is for HPLC and will require no further inputs."),
               img(src = "inputFormat3.png", height = "800px", width = "800px")
        )
      ),
      fluidRow(
        column(12,
               p("Input 5: For breeder data: use MMDDYYYY format. First row first column will be used as the name."),
               p("Second row: always include Breeder ID, Breeder DOB, Pups DOB, #Pups born"),
               img(src = "input5.png", height = "300px", width = "800px")
        )
      )
      
      
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  sheetNames <- reactiveVal(NULL)

  
  metadataColumns <- reactiveValues(data = c("SampleID", "Groups","ProteinConc", "ReconstitutionVolume", "AdjustingValue", "DilutionValue", "comments"))
  originalMetaCols <- reactiveValues(meta_cols = NULL)
  yAxisLimits <- reactiveValues(ymin = NULL, ymax = NULL)
  # Reactive expression to load data
  
  # Observe when a file is uploaded
  observeEvent(input$file1, {
    req(input$file1)
    file_ext <- tools::file_ext(input$file1$name)
    
    if (file_ext %in% c("xls", "xlsx")) {
      sheets <- readxl::excel_sheets(input$file1$datapath)
      sheetNames(sheets)  # Store the sheet names
    } else {
      sheetNames(NULL)  # Reset if not an Excel file
    }
  })
  
  # Render UI for selecting a sheet if applicable
  output$sheetSelector <- renderUI({
    sheets <- sheetNames()
    if (!is.null(sheets) && length(sheets) > 1) {
      selectInput("selectedSheet", "Select a Sheet:", choices = sheets, selected = sheets[1])
    }
  })
  
  
  # Reactive function to read the data
  dataInput <- reactive({
    req(input$file1)
    file_ext <- tools::file_ext(input$file1$name)
    print(file_ext)
    if (file_ext %in% c("xls", "xlsx")) {
      sheets <- excel_sheets(input$file1$datapath)  # Get sheet names
      selected_sheet <- ifelse(length(sheets) == 1, sheets[1], input$selectedSheet)  # Use first sheet if only one
      
      req(selected_sheet)  # Ensure a sheet is selected
      print(selected_sheet)
      
      data <- read_excel(input$file1$datapath, sheet = selected_sheet)
      print(head(data))
    } else {
      data <- read.csv(input$file1$datapath)
      print("is csv?")
      print(head(data))
    }
    print(paste("TEST268"))
    # Remove "comments" column if it exists
    if ("comments" %in% names(data)) {
      data <- data[ , !names(data) %in% "comments"]
    }
    if(length(grep("weight",colnames(data)[1]))>0){
      #data<-read_excel("~/Downloads/weightinput_example1.xlsx")
      dateCols<-grep("^\\d",data[1,])
     # data1<-data2
      colnames(data) <- as.character(data[1, ])  # Set first row as column names
      data <- data[-1, ]  # Remove first row from data
      
      # Convert weight columns to numeric
      data <- data %>%
        pivot_longer(cols = dateCols,  # Select columns starting with 45 (dates)
                     names_to = "weigh_date", 
                     values_to = "weight") %>%
        mutate(
          weigh_date = as.numeric(weigh_date),  # Convert weigh_date to numeric
          weight = as.numeric(weight)  # Convert weight to numeric
        )
      data$weigh_date<-as.numeric(data$weigh_date)
      data$`Pups DOB`<-as.numeric(data$`Pups DOB`)
      data$Groups<-floor((data$weigh_date-data$`Pups DOB`)/7)
      data<-data[,-grep("Pups DOB|weigh_date",colnames(data))]
      meta_cols <- grep("^META", names(data), value = TRUE)
      
      if(length(meta_cols)>0){
        originalMetaCols$meta_cols <- meta_cols 
        # Rename columns in the data
        clean_meta_cols <- gsub("^META", "", meta_cols)
        names(data)[names(data) %in% meta_cols] <- clean_meta_cols
        
        # Update metadataColumns$data with unique values
        metadataColumns$data <- unique(c(metadataColumns$data, clean_meta_cols))
      }
      
      
      print(head(data, 50))
      
      return(data)
    }else{if (any(grepl("breeder", tolower(data[1, ])))) {
    
      print(head(data, 50))
      
      return(data)
    }else {
    # Identify columns that start with "META"
    meta_cols <- grep("^META", names(data), value = TRUE)
    
    if(length(meta_cols)>0){
      originalMetaCols$meta_cols <- meta_cols 
    # Rename columns in the data
    clean_meta_cols <- gsub("^META", "", meta_cols)
    names(data)[names(data) %in% meta_cols] <- clean_meta_cols
    
    # Update metadataColumns$data with unique values
    metadataColumns$data <- unique(c(metadataColumns$data, clean_meta_cols))
    }
    
    
    print(head(data, 50))
    
    return(data)
    }
    }
  })
  
  

  metaInfo <- reactive({
    req(originalMetaCols$meta_cols)  # Use the original names
    clean_meta_cols <- gsub("^META", "", originalMetaCols$meta_cols)
    print(paste("#########################################",clean_meta_cols))
    return(c("Groups",clean_meta_cols))
  })
  
  # Dynamically generate checkboxes based on metaInfo
  output$checkboxesUI <- renderUI({
    req(metaInfo())
    checkboxGroupInput("selectedMeta", "Select Metadata Columns", choices = metaInfo())
  })
  
  # Create dynamic inputs for each sample column
  output$dynamicInputs <- renderUI({
    data <- dataInput()  # Load the uploaded data
    if (any(grepl("breeder", tolower(data[1, ])))) {
      
    }else{if((any(colnames(data) %in% c("ReconstitutionVolume","AdjustingVolume"))&any(data$SampleID %in% "Standard"))){
      
    }else{
    # Get the sample columns excluding 'SampleID' and 'Groups'
      sample_cols <- names(data)[!(names(data) %in% metadataColumns$data)]    # Generate numeric inputs for each sample column (Dilution, Reconstitution, Injection)
     print(paste(sample_cols))
       if((length(which(metadataColumns$data %in% names(data)))>2&!any(names(data) %in% "ProteinConc"))){
        input_list <- lapply(sample_cols, function(sample_col) {
          tagList(
            h4(paste("Inputs for", sample_col)),
            
            numericInput(inputId = paste0("dilution_", sample_col), 
                         label = "Dilution Factor", 
                         value = 1, 
                         min = 1)
          )
        })
        
        # Combine all numeric inputs into a single UI element
        do.call(tagList, input_list)
      }else{
        input_list <- lapply(sample_cols, function(sample_col) {
          fluidRow(
            column(2, h5(sample_col)),  # Display the sample name in the first column
            
            column(2, numericInput(inputId = paste0("dilution_", sample_col), 
                                   label = "Dilution Factor", 
                                   value = 1, min = 1)),
            
            column(2, numericInput(inputId = paste0("reconstitution_", sample_col), 
                                   label = "Reconst. Volume", 
                                   value = 300, min = 0)),
            
            column(2, numericInput(inputId = paste0("injection_", sample_col), 
                                   label = "Inject. Volume", 
                                   value = 25, min = 0)),
            
            column(2, numericInput(inputId = paste0("adjusting_", sample_col), 
                                   label = "Adjust. Value", 
                                   value = 5, min = 1))
          )
        })
        
        # Combine all rows into a single UI element
        do.call(tagList, input_list)
    }
  }}
  })
  
  # Reactive expression to store the fitted models
  lmFits <- reactive({
    data <- dataInput()
    if (any(grepl("breeder", tolower(data[1, ])))) {fits<-NULL}else{if(!any(data$SampleID %in% "Standard")){
      fits<-NULL
    }else{
    # Get the names of the sample columns excluding 'SampleID' and 'Groups'
    #sample_cols <- names(data)[!(names(data) %in% c("SampleID", "ProteinConc", "Groups"))]
    sample_cols <- names(data)[!(names(data) %in% metadataColumns$data)]
    print(paste("LINE273",sample_cols))
    # Create a list of fitted models for each sample column
    fits <- lapply(sample_cols, function(sample_col) {
      standard_data <- data[data$SampleID == "Standard", ]
      standard_data$Groups <- as.numeric(standard_data$Groups)
      
      # Perform log transformation for the sample column
      standard_data$log_Groups <- log10(standard_data$Groups)
      standard_data$log_sample <- log10(standard_data[[sample_col]])
      
      # Fit linear model
      lm_fit <- lm(log_sample ~ log_Groups, data = standard_data)
      
      return(lm_fit)
    })
    
    # Set the names of the fitted models to match the sample columns
    names(fits) <- sample_cols
  }}
    return(fits)
  })
  
  # Generate standard curve plots
  output$curvePlots <- renderUI({
    fits <- lmFits()
    if(is.null(fits)){}else{
    sample_cols <- names(fits)
    
    # Create a list to store rows of plots
    plot_output_list <- list()
    
    # Loop through the sample columns and add plotOutput elements
    for (i in seq_along(sample_cols)) {
      if (i %% 3 == 1) {
        # Start a new fluidRow every 3 plots
        plot_output_list[[length(plot_output_list) + 1]] <- fluidRow()
      }
      
      # Add plotOutput to the current row
      plot_output_list[[length(plot_output_list)]] <- tagAppendChild(
        plot_output_list[[length(plot_output_list)]], 
        column(4, plotOutput(paste0("plot_", sample_cols[i]), height = 300, width = 300))
      )
    }
    
    # Return the list of plot rows as UI elements
    do.call(tagList, plot_output_list)
  }
  })
  
  # Create individual renderPlot calls for each dynamically created plot
  observe({
    fits <- lmFits()
    if(is.null(fits)){}else{
    sample_cols <- names(fits)
    
    for (sample_col in sample_cols) {
      local({
        col <- sample_col  # Local variable to capture the sample_col value
        
        output[[paste0("plot_", col)]] <- renderPlot({
          data <- dataInput()
          standard_data <- data[data$SampleID == "Standard", ]
          standard_data$Groups <- as.numeric(standard_data$Groups)
          standard_data$log_Groups <- log10(standard_data$Groups)
          standard_data$log_sample <- log10(standard_data[[col]])
          
          lm_fit <- fits[[col]]
          formula <- sprintf("y = %.2fx + %.2f", coef(lm_fit)[2], coef(lm_fit)[1])
          r_squared <- summary(lm_fit)$r.squared
          r_squared_text <- sprintf("R² = %.4f", r_squared)
          
          ggplot(standard_data, aes(x = log_Groups, y = log_sample)) +
            geom_point() +
            geom_smooth(method = "lm", se = FALSE, color = "blue") +
            labs(
              title = paste(col, "(Log-Log Scale)"),
              x = "Log10 Concentration (log10 umol)",
              y = paste("Log10", col, "AUC")
            ) +
            theme_minimal() +
            annotate("text", x = min(standard_data$log_Groups, na.rm = TRUE), 
                     y = max(standard_data$log_sample, na.rm = TRUE), 
                     label = formula, hjust = 0, vjust = 1.5, size = 5, color = "black") +
            annotate("text", x = min(standard_data$log_Groups, na.rm = TRUE), 
                     y = max(standard_data$log_sample, na.rm = TRUE) - 0.2, 
                     label = r_squared_text, hjust = 0, vjust = 1.5, size = 5, color = "black")
        })
      })
    }
    }
  })
  
  # Reactive expression to store filtered groups
  selectedGroups <- reactiveVal()
  selectedSamples <- reactiveVal()
  # Render group selector UI
  output$groupSelector <- renderUI({
    data <- dataInput()
    sample_cols <- names(data)[!(names(data) %in% metadataColumns$data)]
    ignores<-unique(data$Groups[data$SampleID=="Standard"])
    toSearch<-unique(data$Groups)[!(unique(data$Groups) %in% ignores)]
    toSearch<-c(sample_cols,toSearch)
    selectInput("groupsFilter", "Select Groups to Display:", 
                choices = toSearch, multiple = TRUE)
    
  })
  
  # Observe group filtering
  observeEvent(input$applyFilter, {
    selectedGroups(input$groupsFilter)
  })
  
  # samples able to select and remove
  
  output$sampleRemover <- renderUI({
    data <- dataInput()
    toSearch<-unique(data$SampleID)[!(unique(data$SampleID)=="Standard")]
    selectInput("samplesFilter", "Select SampleIDs to remove:", 
                choices = toSearch, multiple = TRUE)
  })
  
  # Observe group filtering
  observeEvent(input$applysampleRemover, {
    selectedSamples(input$samplesFilter)
  })
  
  # Reactive value for transformed data
  transformedData <- reactiveVal()
  filteredDataForPlot<- reactiveVal()
  # Observe input confirmation
  observeEvent(input$confirmValues, {
    data <- dataInput()
    fits <- lmFits()
    #sample_cols <- names(data)[!(names(data) %in% c("SampleID", "ProteinConc", "Groups"))]
    sample_cols <- names(data)[!(names(data) %in% metadataColumns$data)]
    
    transformed_data <- data
    if (any(grepl("breeder", tolower(data[1, ])))) {}else{
    transformed_data <-transformed_data[transformed_data$SampleID != "Standard",]
    print(head(transformed_data,20))
    
    if(is.null(fits)){
      transformed_data <- data
      transformed_data <-transformed_data[transformed_data$SampleID != "Standard",]
      for (sample_col in sample_cols) {
        dilution <- as.numeric(input[[paste0("dilution_", sample_col)]])
        transformed_data[[sample_col]] <- round(transformed_data[[sample_col]]*dilution,4)
        }
    
      
    }
    if(!is.null(fits)){
      for (sample_col in sample_cols) {
      lm_fit <- fits[[sample_col]]
      if(any(names(transformed_data) %in% c("ReconstitutionVolume","AdjustingValue"))){
        #dilution <- as.numeric(input[[paste0("dilution_", sample_col)]])
        dilution<-1
        reconstitution <- as.numeric(transformed_data$ReconstitutionVolume)
        #injection <- as.numeric(input[[paste0("injection_", sample_col)]])
        injection<-25
        AdjustingValue<-transformed_data$AdjustingValue
      }
      if(!any(names(transformed_data) %in% c("ReconstitutionVolume","AdjustingValue"))){
        dilution <- as.numeric(input[[paste0("dilution_", sample_col)]]) 
        reconstitution <- as.numeric(input[[paste0("reconstitution_", sample_col)]]) 
        injection <- as.numeric(input[[paste0("injection_", sample_col)]])
        AdjustingValue<-as.numeric(input[[paste0("adjusting_", sample_col)]])
      }
      #print(paste(dilution, reconstitution, injection, AdjustingValue))
      #print(lm_fit)
      #print(paste("LINE 300"))
      # Convert log10 intensity to concentration using the inverse of the linear model formula
      tempdf <- transformed_data
      tempdf$AmountExtracted <- round((10^((log10(tempdf[[sample_col]]) - coef(lm_fit)[1]) / coef(lm_fit)[2])),4)
      
      # Use user inputs in the concentration calculation
      tempdf$in300uLmM <- round(dilution*(tempdf$AmountExtracted * reconstitution) / injection,4)
      tempdf$adjustedToSample <- round(tempdf$in300uLmM*AdjustingValue,4)
      if(any(names(transformed_data) %in% "ProteinConc")){
        tempdf$umolg <- round((1000 * tempdf$adjustedToSample / tempdf$ProteinConc),4)
        
      }else{
        tempdf$umolg <- round(tempdf$adjustedToSample,4)
        
      }
      
      # Append the calculated umolg back to the transformed_data
      transformed_data[[paste0(sample_col, "_AmountExtracted")]] <- tempdf$adjustedToSample
      transformed_data[[paste0(sample_col, "_umolg")]] <- tempdf$umolg
    }}}
    
    print(head(transformed_data))
    transformedData(transformed_data)
    
    
  })
  
  # Render DataTable only after data has been uploaded
  output$dataPreview <- DT::renderDataTable({
    req(transformedData())  # Ensure data exists before rendering
    transformedData()
  })
  
  # Render box plot only after data has been uploaded
  # Define a reactive plot object
  reactivePlot <- reactive({
    fits <- lmFits()
    req(transformedData())  # Ensure data exists before plotting
    data <- transformedData()
    if (any(grepl("breeder", tolower(data[1, ])))) {
       print(data)
      print(paste("this is breeder data"))
          # Example code to create a bar plot (adjust as needed)
          breeder_data <- data[2:nrow(data), ]
          colnames(breeder_data)<-as.character(data[1,]) # Exclude the header row
          breeder_data$AgeOfMother<-round((as.numeric(breeder_data$`Pups DOB`)-as.numeric(breeder_data$`Breeder DOB`))/7)
          breeder_data<-breeder_data[breeder_data$AgeOfMother>0,]
          breeder_data<-breeder_data[,c("AgeOfMother","#Pups born","Breeder ID")] # Assuming first column contains breeder data
          breeder_data$`#Pups born` <- as.numeric(breeder_data$`#Pups born`)
          
          summary_data <- breeder_data %>%
            group_by(AgeOfMother) %>%
            summarise(
              mean_pups = mean(`#Pups born`, na.rm = TRUE),
              se_pups = sd(`#Pups born`, na.rm = TRUE) / sqrt(n()),
              .groups = "drop"
            )
          print(paste("this is summary data"))
          print(head(summary_data,100))

          # Define automatic limits
          y_max_auto <- max(breeder_data$'#Pups born', na.rm = TRUE)
          x_min_auto <- min(breeder_data$AgeOfMother, na.rm = TRUE)
          x_max_auto <- max(breeder_data$AgeOfMother, na.rm = TRUE)
          
          # Handle potential NULL or NA values for user input
          y_max <- if (!is.null(input$y_max) && !is.na(input$y_max) && input$y_max > 0) input$y_max else y_max_auto
          x_min <- if (!is.null(input$x_min) && !is.na(input$x_min) && input$x_min >= 0) input$x_min else x_min_auto
          x_max <- if (!is.null(input$x_max) && !is.na(input$x_max) && input$x_max > x_min) input$x_max else x_max_auto
          
          
          # Generate plot
          # p <- ggplot(breeder_data, aes(x = AgeOfMother, y = .data['#Pups born'])) +
          #   geom_bar(stat = "identity", fill = "skyblue", color = "black", width = 0.7) +  # Bar plot
          #   geom_errorbar(aes(ymin = mean_pups - se_pups, ymax = mean_pups + se_pups), width = 0.2) +  # Error bars
          #   labs(
          #     title = paste0("Pups Born by Age of Mother ", colnames(data)[1]),
          #     x = "Age of Mother (weeks)",
          #     y = "# Pups Born"
          #   ) +
          #   theme_minimal() +
          #   coord_cartesian(ylim = c(0, y_max), xlim = c(x_min, x_max))  # Adjust axis limits dynamically
          # 
          # # Add jitter if requested
          # if (input$showDots) {
          #   p <- p + geom_jitter(size = input$dotSize, width = 0.2)
          # }
          ############################
          
          p <- ggplot(summary_data, aes(x = AgeOfMother, y = mean_pups)) +
            geom_bar(stat = "identity", fill = "skyblue", color = "black") +
            geom_errorbar(aes(ymin = mean_pups - se_pups, ymax = mean_pups + se_pups), width = 0.2) +
            labs(title = paste0("Pups Born by Age of Mother ", colnames(data)[1]), x = "Age of Mother (weeks)",
                 y = "# Pups Born") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+ coord_cartesian(ylim = c(0, y_max), xlim = c(x_min, x_max)) 
          
          # if (!is.null(input$ymin) && !is.null(input$ymax)) {
          #   p <- p + coord_cartesian(ylim = c(input$ymin, input$ymax))
          # }
          
          if (input$showDots) p <- p + geom_jitter(data = breeder_data, aes(x = AgeOfMother, y = `#Pups born`), 
                                                   size = input$dotSize, width = 0.2, color = "red")
          if (input$showNames) p <- p + geom_text(data = breeder_data, aes(x = AgeOfMother, y = `#Pups born`, label = `Breeder ID`), 
                                                  hjust = -0.2, size = 5,color = "red")
          
                  ############################
          print(paste("barplot should be showing"))
          return(p)  # Return the plot as a reactive object
        
    }else{
    if (is.null(fits)){
      data <- data
      print(paste("LINE488"))
      print(head(data))
    } else {
      data <- cbind(data[, grep("Extracted|SampleID|Groups|ProteinConc|Reconstitution|Adjusting|Dilution", colnames(data))],
                    data[, which(colnames(data) %in% gsub("META", "", originalMetaCols$meta_cols))])
    }
    
    sample_cols <- gsub("_AmountExtracted", "", names(data)[!(names(data) %in% metadataColumns$data)])
    sampleRemove <- selectedSamples()
    if (!is.null(sampleRemove)){
      data <- data[!(data$SampleID %in% sampleRemove), ]
    }
    
    # Filter data based on selected groups
    grouping <- selectedGroups()
    print(paste("LINE502"))
    print(grouping)
       if (is.null(selectedGroups())){
      data_filtered <- data
    } else {
      
      if (length(grouping) == 1 & (any(colnames(data) %in% paste0(grouping, "_AmountExtracted")) | any(colnames(data) %in% grouping))){
        if (is.null(fits)){
          #samples <- c(1:dim(data)[2])
          if(any(colnames(data) %in% grouping)){
            samples <- c(which(colnames(data) %in% metadataColumns$data), which(colnames(data) %in% grouping))
          }else{
            samples <- c(which(colnames(data) %in% metadataColumns$data), which(colnames(data) %in% paste0(grouping, "_AmountExtracted")))
          }
        } else {
          if(any(colnames(data) %in% grouping)){
            samples <- c(which(colnames(data) %in% metadataColumns$data), which(colnames(data) %in% grouping))
          }else{
            samples <- c(which(colnames(data) %in% metadataColumns$data), which(colnames(data) %in% paste0(grouping, "_AmountExtracted")))
          }
        }
        data_filtered <- data[, samples]
      } else {
        
        groups <- unique(data$Groups[data$Groups %in% grouping])
        print(paste("####################################"))
        print(paste("####################################"))
        print(paste("####################################"))
        print(paste("####################################"))
        print(paste("LINE 554grouping"))
        print(grouping)
        print(paste("####################################"))
        print(paste("####################################"))
        print(paste("####################################"))
        print(paste("####################################"))
        print(paste("LINE 554datagroup"))
        print(unique(data$Groups))
        print(paste("####################################"))
        print(paste("####################################"))
        print(paste("####################################"))
        print(paste("####################################"))
        if(any(colnames(data) %in% grouping)){
          print(paste("LINE545 TEST ^^^^^^^^^^^^^^^^^^^^^^"))
          if(!any(data$Groups %in% grouping)){samples <- c(which(colnames(data) %in% metadataColumns$data), which(colnames(data) %in% grouping))}else{
            print(paste("LINE545.B TEST ^^^^^^^^^^^^^^^^^^^^^^"))
            print(colnames(data))
            print(grouping)
            cols<-which(data$Groups %in% grouping)
            if(any(grouping %in% sample_cols)){samples <- c(which(colnames(data) %in% metadataColumns$data), which(colnames(data) %in% paste0(grouping, "_AmountExtracted")))}else{
              samples <- c(which(colnames(data) %in% metadataColumns$data),which(colnames(data) %in% sample_cols))
            }
            
            
            
            
          }
                  }
        
        #samples <- c(which(colnames(data) %in% metadataColumns$data), grep(matched_columns,colnames(data)))
        if(any(grouping %in% sample_cols)){
          if(length(which(colnames(data) %in% paste0(grouping, "_AmountExtracted")))==0){samples <- c(which(colnames(data) %in% metadataColumns$data), which(colnames(data) %in% grouping))}else{
            samples <- c(which(colnames(data) %in% metadataColumns$data), which(colnames(data) %in% paste0(grouping, "_AmountExtracted")))
          }
          
          
        
        }else{
          samples <- c(which(colnames(data) %in% metadataColumns$data),which(colnames(data) %in% sample_cols))
        }
        
        print(paste("####################################"))
        print(paste("LINE 555"))
        print(sample_cols)
        print(paste("Samples"))
        print(samples)
        print(paste("####################################"))
       
        #data_filtered <- data[data$Groups %in% groups, samples]
        if(length(groups)==0){data_filtered <- data[, samples]}else{
          
          data_filtered <- data[which(data$Groups %in% groups), samples]
        }
        #data_filtered <- data[, samples]
        
        print(paste("####################################"))
        print(paste("####################################"))
        print(paste("####################################"))
        print(paste("####################################"))
        print(paste("LINE 556"))
        print(head(data_filtered))
        print(paste("####################################"))
        print(paste("####################################"))
        print(paste("####################################"))
        print(paste("####################################"))
      }
    }
    
    # Combine selected metadata columns into a new column
    selectedMeta <- input$selectedMeta
    if (!is.null(selectedMeta) && length(selectedMeta) > 0) {
      data_filtered$CombinedMeta <- apply(data_filtered[, selectedMeta, drop = FALSE], 1, paste, collapse = "_")
    } else {
      data_filtered$CombinedMeta <- data_filtered$Groups
    }
    filteredDataForPlot(data_filtered)
    # Melt data to long format for plotting
    melted_data <- melt(data_filtered, id.vars = c(colnames(data)[colnames(data) %in% metadataColumns$data], "CombinedMeta"))
    
    # Set y-axis limits
    yAxisLimits$ymin <- round(min(melted_data$value, na.rm = TRUE), 3)
    yAxisLimits$ymax <- round(max(melted_data$value, na.rm = TRUE), 3)
    
    concTitle <- if (is.null(fits)) "concentration (raw value)" else "concentration in 200uL(mM)"
    if(length(grep("weight",colnames(data)))>0){concTitle <- "weight (g)"
    # Summarize data by Age (Groups) and CombinedMeta (for different lines)
    summary_data <- melted_data %>%
      group_by(Groups, CombinedMeta) %>%
      summarise(
        mean_value = mean(value, na.rm = TRUE),
        sd_value = sd(value, na.rm = TRUE),
        se_value = sd_value / sqrt(n()),  # Standard error
        .groups = "drop"
      )
    # Base plot: Line plot with error bars
    p <- ggplot(summary_data, aes(x = Groups, y = mean_value, color = CombinedMeta, group = CombinedMeta)) +
      geom_line(size = 1) +  # Line for each CombinedMeta group
      geom_point(size = 3) +  # Points for each group
      geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.3) +  # Error bars
      theme_minimal() +
      labs(title = "Sample Concentrations by Group", x = "Age (Weeks)", y = concTitle, color = "Sample") +
      scale_x_continuous(breaks = unique(summary_data$Groups))  # Ensure all ages appear on the x-axis
    
    # Add jittered dots for individual samples
    if (input$showDots) {
      p <- p + geom_jitter(data = melted_data, aes(x = Groups, y = value, color = CombinedMeta), 
                           size = input$dotSize, width = 0.2, alpha = 0.7)
    }
    
    # Add sample labels
    if (input$showNames) {
      p <- p + geom_text(data = melted_data, aes(x = Groups, y = value, label = SampleID, color = CombinedMeta), 
                         hjust = -0.2, size = 3)
    }
    
    # Apply y-axis limits if provided
    if (!is.null(input$ymin) && !is.null(input$ymax)) {
      p <- p + coord_cartesian(ylim = c(input$ymin, input$ymax))
    }
    
    # Return the plot
    return(p)
    }else{
    # Generate plot
    p <- ggplot(melted_data, aes(x = CombinedMeta, y = value, color = variable)) +
      geom_boxplot() +
      labs(title = "Sample Concentrations by Group", y = concTitle, color = "Sample") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    if (!is.null(input$ymin) && !is.null(input$ymax)) {
      p <- p + coord_cartesian(ylim = c(input$ymin, input$ymax))
    }
    
    if (input$showDots) p <- p + geom_jitter(size = input$dotSize, width = 0.2)
    if (input$showNames) p <- p + geom_text(aes(label = SampleID), hjust = -0.2)
    
    return(p)  # Return the plot as a reactive object
    }
  }
  })
  
  
  
  output$boxPlot <- renderPlot({reactivePlot()})
  
  # Update y-axis controls after plot is generated
  output$yAxisControls <- renderUI({
    req(yAxisLimits$ymin, yAxisLimits$ymax) # Ensure limits are computed
    tagList(
      div(
        class = "row",
        div(
          class = "col-sm-6",
          numericInput("ymin", "Y-axis Min:", value = yAxisLimits$ymin, step = 0.5)
        ),
        div(
          class = "col-sm-6",
          numericInput("ymax", "Y-axis Max:", value = yAxisLimits$ymax, step = 0.5)
        )
      ),
      actionButton("resetYAxis", "Reset Y-Axis")
    )
  })
  
  # Reset y-axis limits when reset button is clicked
  observeEvent(input$resetYAxis, {
    updateNumericInput(session, "ymin", value = yAxisLimits$ymin)
    updateNumericInput(session, "ymax", value = yAxisLimits$ymax)
  })
  
  # Display default y-axis range for reference
  output$defaultYRange <- renderText({
    req(yAxisLimits$ymin, yAxisLimits$ymax)
    paste("Default Y-axis Range: Min =", yAxisLimits$ymin, ", Max =", yAxisLimits$ymax)
  })
  
  output$downloadPlotUI <- renderUI({
    req(reactivePlot())  # Ensure that the plot is available before showing UI elements
    
    tagList(
      downloadButton("downloadPlot", "Download Plot as PDF"),
      textInput("filename", "Enter file name (without extension):", value = "boxplot"),
      numericInput("pdfWidth", "PDF Width:", value = 8, min = 1, step = 0.5),
      numericInput("pdfHeight", "PDF Height:", value = 8, min = 1, step = 0.5)
    )
  })
  
  # Download Handler
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0(input$filename, ".pdf")  # Use the filename from user input
    },
    content = function(file) {
      pdf(file, width = input$pdfWidth, height = input$pdfHeight,useDingbats = FALSE)
      print(reactivePlot())
      dev.off()
    }
    
  )
  
  # Render p-value table
  output$pvalueTable <- renderDT({
    fits <- lmFits()
    req(transformedData())  # Ensure data exists before plotting
    data <- transformedData()
    data2<-filteredDataForPlot()
    
    print(paste("This is the pvalueStuff"))
    print(head(data))
    print(paste("newinput"))
    print(head(data2))
    
if (any(grepl("breeder", tolower(data[1, ])))) {filtered_sorted_table<-NULL}else{
    data_filtered <- data2
    # Melt data to long format for plotting
    melted_data <- melt(data_filtered, id.vars = c(colnames(data)[colnames(data) %in% metadataColumns$data], "CombinedMeta")) 
    selectedMeta <- input$selectedMeta
    #print(paste("Line795"))
    #print(!is.null(selectedMeta))
    #print(length(selectedMeta) > 0)
    #print((selectedMeta %in% "Groups"))
    if (!is.null(selectedMeta) && length(selectedMeta) > 0  && any(selectedMeta %in% "Groups")) {
      melted_data$groupsIncluded <- "T"
    } else {
      melted_data$groupsIncluded <- "F"
    }
    print(paste(table(melted_data$CombinedMeta)))
    print(head(melted_data))
    # Compute p-values and prepare comparison table
    p_values <- compute_pvalues(melted_data)
    filtered_sorted_table <- prepare_comparison_table(p_values)  # Filter and sort
}   
    datatable(filtered_sorted_table, options = list(pageLength = 10, autoWidth = TRUE))

  })  
}

# Run the application
shinyApp(ui = ui, server = server)
