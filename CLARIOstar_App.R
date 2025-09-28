library(shiny)
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(plotly)
library(stringr)
library(ggplotify)
library(statmod)
library(nlme)

ui <- navbarPage(
  titlePanel("CLARIOstar Absorbance Visualization"),
  
  tabPanel("Plotting",
           sidebarLayout(
             sidebarPanel(
               fileInput("file", "Upload Excel File:", accept = ".xlsx"),
               fileInput("file2", "Upload Additonal Excel File (repeatable):", accept = ".xlsx"),
               actionButton("uncheck_all", "Uncheck All Samples"),
               uiOutput("sample_selector")
             ),
             
             mainPanel(
               plotlyOutput("absorbance_plot"),
               checkboxInput("toggle_average", "Show Average with Error Bars", TRUE),
               checkboxInput("toggle_mean_only", "Show Only Mean Values", FALSE),
               checkboxInput("toggle_groupcurves", "Show Growth Curve Group Means (Run Pairwise Comparison First)", FALSE),
               numericInput("max_timepoint", "Max Timepoint to Include:", value = NA, min = 0),
               actionButton("plot_button", "Generate Plot"),
               # Update data table by removing outlier samples manually
               fluidRow(
                 column(6,
                        h3("Remove Selected Sample by ID and column #"),
                        uiOutput("sampleRemover"),
                        actionButton("applysampleRemover", "Remove Selected Sample"),
                        downloadButton("download_plot", "Download Plot as PDF"),
                        numericInput("pdfWidth", "PDF Width:", value = 12, min = 1, step = 0.5),
                        numericInput("pdfHeight", "PDF Height:", value = 8, min = 1, step = 0.5),
                        uiOutput("x_range_inputs"),
                        uiOutput("y_range_inputs"),
                        sliderInput("point_size", "Point Size:", min = 0.25, max = 5, value = 1, step = 0.25)
                 ),
                 column(6,
                        h3("Pairwise Trend Comparison (p-values)"),
                        h5("May take 60+ seconds with many samples"),
                        
                        numericInput("nsim_input", "Number of Permutations for Growth Curve Comparison (nsim)", 
                                     value = 50, min = 10, max = 2000),
                        helpText("It will take very long to compare many samples when using a large n of permutations"),
                        actionButton("compare_button", "Run Pairwise Comparison"),
                        br(), br(),
                        dataTableOutput("pairwise_pvalues")
                 )
               )
             )
           )
  ),
  tabPanel(
    "How To Use",
    h3("How to Use This Application"),
    h4("The input file type required is an .xlsx file output from the MARS software."),
    
    fluidRow(
      column(12,
             p("Export the excel file using the following selections:"),
             img(src = "clariostarInput.png", height = "700px", width = "700px")
      )
    ),
    
    fluidRow(
      column(12,
             p("How to filter: You can hover over your lines to see which sample , row and COLUMN # they belong to then remove it:"),
             img(src = "clariostarHOWTO.png", height = "700px", width = "700px"),
             p("You can add additional files at any time through the second file input. Add one new one at a time."),
             p("Unselect all samples and select a few before running pairwise comparisons, use a smaller nsim for more samples."),
             p(
               "Using a large number will be more accurate but be much slower. Both growcurve from ",
               a("statmod", href = "https://cran.r-project.org/web/packages/statmod/statmod.pdf", target = "_blank"),
               " and anova from",
               a("stats", href = "https://r-universe.dev/manuals/stats.html#anova.lm", target = "_blank"),
               " p values are calculated."
             ),
             p("p values for statmod will fluctuate as they are based on simulated modelling and not a set value.")
      )
    )
    # 
    # fluidRow(
    #   column(12,
    #          p("Input 3: This input type is for HPLC and will require no further inputs."),
    #          img(src = "inputFormat3.png", height = "800px", width = "800px")
    #   )
    # )
  )
)

##########################
compareGrowthCurves <- function(group, y, levels = NULL, nsim = 100, fun = meanT, times = NULL, 
                                verbose = TRUE, adjust = "holm", n0 = 0.5) {
  group <- as.character(group)
  
  if (is.null(levels)) {
    tab <- table(group)
    tab <- tab[tab >= 2]
    lev <- names(tab)
  } else {
    lev <- as.character(levels)
  }
  
  nlev <- length(lev)
  if (nlev < 2) stop("Less than 2 groups to compare")
  
  if (is.null(dim(y))) stop("y must be matrix-like")
  y <- as.matrix(y)
  if (!is.null(times)) y <- y[, times, drop = FALSE]
  
  g1 <- g2 <- character()
  stat <- pvalue <- numeric()
  
  for (i in 1:(nlev - 1)) {
    for (j in (i + 1):nlev) {
      if (verbose) cat(lev[i], lev[j])
      
      sel <- group %in% c(lev[i], lev[j])
      
      result <- tryCatch({
        out <- compareTwoGrowthCurves(group[sel], y[sel, , drop = FALSE], 
                                      nsim = nsim, fun = fun, n0 = n0)
        if (verbose) cat(" ", round(out$stat, 2), "\n")
        list(stat = out$stat, pval = out$p.value)
      }, error = function(e) {
        if (verbose) cat(" -- skipped due to error:", conditionMessage(e), "\n")
        list(stat = NA, pval = NA)
      })
      
      g1 <- c(g1, lev[i])
      g2 <- c(g2, lev[j])
      stat <- c(stat, result$stat)
      pvalue <- c(pvalue, result$pval)
    }
  }
  
  tab <- data.frame(Group1 = g1, Group2 = g2, Stat = stat, P.Value = pvalue)
  tab$adj.P.Value <- p.adjust(pvalue, method = adjust)
  return(tab)
}
###########################

server <- function(input, output, session) {
  upload_counter <- reactiveVal(1)
  uploaded_files <- reactiveValues(extra_files = list(), counter = 2)
  
  observeEvent(input$uncheck_all, {
    updateCheckboxGroupInput(session, "samples", selected = character(0))
  })
  
  statmod_results <- eventReactive(input$compare_button, {
    req(filtered_data(), input$samples)
    library(nlme)
    
    dat <- filtered_data() |>
      filter(Sample_ID %in% input$samples) |>
      mutate(replicate = paste(Sample_ID, Row, Column, sep = "_"))
    
    if (!is.null(input$max_timepoint) && !is.na(input$max_timepoint)) {
      dat <- dat %>% filter(timepoint <= input$max_timepoint)
    }
    
    dat <- dat[dat$Sample_ID != "Empty", ]
    
    wide <- dat |>
      select(replicate, timepoint, Absorbance) |>
      pivot_wider(names_from = timepoint, values_from = Absorbance)
    
    y_mat <- as.matrix(wide |> select(-replicate))
    ord <- order(as.numeric(colnames(y_mat)))
    y_mat <- y_mat[, ord]
    colnames(y_mat) <- colnames(y_mat)
    
    grp <- factor(str_remove(wide$replicate, "_[A-H]_[0-9]+$"))
    keep_levels <- names(which(table(grp) >= 2))
    y_mat <- y_mat[grp %in% keep_levels, , drop = FALSE]
    grp <- droplevels(grp[grp %in% keep_levels])
    
    # ---- Run compareGrowthCurves ----
      pmat <- compareGrowthCurves(group = grp,
                            y = y_mat,
                            nsim = input$nsim_input,
                            adjust = "BH")
    
    # Rename for clarity
    ptab <- pmat |>
      dplyr::rename(Sample1 = Group1,
                    Sample2 = Group2,
                    p_value = P.Value,
                    adjP = adj.P.Value) |>
      dplyr::mutate(across(c(p_value, adjP), ~ signif(.x, 4)))
    
    # ---- Run ANOVA on each pair ----
    anova_results <- list()
    
    for (i in seq_len(nrow(ptab))) {
      samp1 <- ptab$Sample1[i]
      samp2 <- ptab$Sample2[i]
      
      subdat <- dat %>%
        filter(Sample_ID %in% c(samp1, samp2))
      
      subdat$Sample_ID <- factor(subdat$Sample_ID)
      subdat$replicate <- factor(paste(subdat$Sample_ID, subdat$Row, subdat$Column, sep = "_"))
      
      # Use lme model: y ~ time * group
      tryCatch({
        fit_full <- lme(Absorbance ~ timepoint * Sample_ID,
                        random = ~1 | replicate,
                        data = subdat)
        fit_reduced <- lme(Absorbance ~ timepoint + Sample_ID,
                           random = ~1 | replicate,
                           data = subdat)
        anova_p <- anova(fit_reduced, fit_full)$"p-value"[2]
      }, error = function(e) {
        anova_p <- NA
      })
      
      anova_results[[i]] <- anova_p
    }
    
    # Add ANOVA p-values and BH-adjust
    ptab$anova_p <- signif(unlist(anova_results), 4)
    ptab$anova_adjP <- signif(p.adjust(ptab$anova_p, method = "BH"), 4)
    
    return(ptab)
  })
  
  output$pairwise_pvalues <- renderDataTable({
    req(input$compare_button > 0)
    req(statmod_results())
    statmod_results()
  }, options = list(pageLength = 10))
  
  output$sampleRemover <- renderUI({
    req(processed_data())  # Ensure processed_data() is ready
    data <- processed_data()
    toSearch <- unique(data$uniqueSampleID)
    selectInput(
      inputId = "samplesFilter",
      label = "Select SampleIDs to remove:",
      choices = toSearch,
      multiple = TRUE
    )
  })
  
  selectedSamples <- reactiveVal()
  
  # Observe group filtering
  observeEvent(input$applysampleRemover, {
    selectedSamples(input$samplesFilter)
  })
  
  
  
  filtered_data <- reactive({
    req(processed_data(), input$samples)
    data <- processed_data() %>% dplyr::filter(Sample_ID %in% input$samples)
    
    if (!is.null(input$max_timepoint) && !is.na(input$max_timepoint)) {
      data <- data %>% filter(timepoint <= input$max_timepoint)
    }
    print(dim(data))
    sampleRemove <- selectedSamples()
    if (!is.null(sampleRemove)){
      print(sampleRemove)
      sampleIDs<-unlist(str_extract_all(sampleRemove,pattern="\\w.+(?=_wellLoc)"))
      rows<-unlist(str_extract_all(sampleRemove,pattern="(?<=_wellLoc)[A-Z]+"))
      columns<-unlist(str_extract_all(sampleRemove,pattern="(?<=_wellLoc\\w)\\d+"))
      
      
      for(i in 1:length(sampleIDs)){
        print(sampleIDs[i])
        print(columns[i])
        print(rows[i])
        print(head(data))
        data <- data[!((data$Sample_ID %in% sampleIDs[i])&(data$Column %in% columns[i])&(data$Row %in% rows[i])), ] 
        print(dim(data))
      }
      
      
    }
    return(data)# data <- data[!(data$SampleID %in% sampleRemove)&(), ]
  })
  
  axis_ranges <- reactive({
    req(filtered_data())
    data <- filtered_data()
    print(paste("axis"))
    print(dim(data))
    x_min <- min(data$timepoint, na.rm = TRUE)
    x_max <- max(data$timepoint, na.rm = TRUE)
    y_min <- min(data$Absorbance, na.rm = TRUE)
    y_max <- max(data$Absorbance, na.rm = TRUE)
    list(x_min = x_min, x_max = x_max, y_min = y_min, y_max = y_max)
  })
  
  output$x_range_inputs <- renderUI({
    req(axis_ranges())
    ranges <- axis_ranges()
    tagList(
      numericInput("x_min", "X-Axis Minimum:", value = ranges$x_min),
      numericInput("x_max", "X-Axis Maximum:", value = ranges$x_max)
    )
  })
  
  output$y_range_inputs <- renderUI({
    req(axis_ranges())
    ranges <- axis_ranges()
    tagList(
      numericInput("y_min", "Y-Axis Minimum:", value = ranges$y_min),
      numericInput("y_max", "Y-Axis Maximum:", value = ranges$y_max)
    )
  })
  
  uploaded_data <- reactive({
    req(input$file)
    file_path <- input$file$datapath
    data <- read_excel(file_path, sheet = 1)
    data <- data[, 1:13]  # Adjust columns as necessary
    return(data)
  })
  
  rename_duplicates <- function(ids) {
    new_ids <- ids
    id_counts <- list()
    
    for (i in seq_along(ids)) {
      id <- ids[i]
      if (id %in% names(id_counts)) {
        id_counts[[id]] <- id_counts[[id]] + 1
        new_ids[i] <- paste0(id, "_", id_counts[[id]])
      } else {
        id_counts[[id]] <- 1
      }
    }
    return(new_ids)
  }
  
  process_file <- function(file_path) {
    data <- readxl::read_excel(file_path, sheet = 1)
    data <- data[, 1:13]
    
    absorbance_start <- grep("Blank corrected based", data.frame(data)[, 2])
    sample_ids_start <- grep("Sample IDs", data.frame(data)[, 2])
    
    timepoint <- 0
    full_merged <- NULL
    
    for (i in seq_along(absorbance_start)) {
      absorbance_data <- data[(absorbance_start[i] + 2):(absorbance_start[i] + 9), ]
      sample_ids <- data[(sample_ids_start[i] + 2):(sample_ids_start[i] + 9), ]
      
      colnames(absorbance_data) <- c("Row", as.character(1:12))
      colnames(sample_ids) <- c("Row", as.character(1:12))
      
      absorbance_data <- as.data.frame(absorbance_data)
      absorbance_data[is.na(absorbance_data)] <- 0
      sample_ids <- sample_ids %>% mutate(across(everything(), as.character))
      sample_ids[is.na(sample_ids)] <- "Empty"
      
      absorbance_data <- absorbance_data %>% mutate(across(-Row, as.numeric))
      
      absorbance_long <- absorbance_data %>%
        pivot_longer(cols = -Row, names_to = "Column", values_to = "Absorbance")
      
      sample_ids_long <- sample_ids %>%
        pivot_longer(cols = -Row, names_to = "Column", values_to = "Sample_ID")
      
      merged_data <- left_join(absorbance_long, sample_ids_long, by = c("Row", "Column"))
      merged_data <- data.frame(merged_data)
      merged_data$timepoint <- timepoint
      full_merged <- rbind(full_merged, merged_data)
      timepoint <- timepoint + 10
    }
    
    return(full_merged)
  }
  
  observeEvent(input$file2, {
    req(input$file2)
    
    count <- uploaded_files$counter
    data2 <- process_file(input$file2$datapath)
    
    # Collect all existing Sample_IDs (from file1 and previous file2s)
    existing_ids <- c()
    if (!is.null(input$file)) {
      existing_ids <- unique(process_file(input$file$datapath)$Sample_ID)
    }
    existing_ids <- unique(c(existing_ids,
                             unlist(lapply(uploaded_files$extra_files, function(df) unique(df$Sample_ID)))))
    
    dup_ids <- intersect(existing_ids, unique(data2$Sample_ID))
    
    # Rename duplicates in data2 with _<count>
    for (id in dup_ids) {
      data2$Sample_ID[data2$Sample_ID == id] <- paste0(id, "_", count)
    }
    
    # Store data2 into list
    uploaded_files$extra_files[[paste0("file2_", count)]] <- data2
    uploaded_files$counter <- count + 1
  })
  
  
  processed_data <- reactive({
    req(input$file)
    
    data1 <- process_file(input$file$datapath)
    
    if (length(uploaded_files$extra_files) > 0) {
      data2_all <- do.call(bind_rows, uploaded_files$extra_files)
      combined <- bind_rows(data1, data2_all)
    } else {
      combined <- data1
    }
    
    combined$uniqueSampleID <- paste0(combined$Sample_ID, "_wellLoc", combined$Row, combined$Column)
    return(combined)
  })
  
  
  output$sample_selector <- renderUI({
    req(processed_data())
    vals <- unique(processed_data()$Sample_ID)
    prefix   <- str_extract(vals, "^[^_]+")
    treatment <- str_remove(vals, "^[^_]+_")
    compound <- str_extract(treatment, "^[A-Za-z]+")
    number   <- as.numeric(str_extract(treatment, "\\d+\\.?\\d*"))
    
    # Put non-numeric cases at the end
    number[is.na(number)] <- Inf
    
    # Order based on prefix, compound, then numeric value
    ord <- order(prefix, compound, number)
    sorted_vals <- vals[ord]
    sample_ids<-sorted_vals
    checkboxGroupInput("samples", "Select Samples to Include:", choices = sample_ids, selected = sample_ids)
  })
  
  plot_data <- reactive({
    req(filtered_data(), input$samples)
    data <- filtered_data()
    filtered_data <- data %>% dplyr::filter(Sample_ID %in% input$samples)
    
    if (input$toggle_average) {
      summary_data <- filtered_data %>%
        group_by(Sample_ID, timepoint) %>%
        summarise(
          mean_absorbance = mean(Absorbance, na.rm = TRUE),
          se_absorbance = sd(Absorbance, na.rm = TRUE) / sqrt(n()),
          .groups = "drop"
        )
      return(summary_data)
    } else {
      return(filtered_data)
    }
  })
  
  reactive_plot <- reactive({
    req(plot_data(), input$plot_button)
    
    # Use QurvE plotGrowthCurves if the user wants group means
    if (input$toggle_groupcurves && input$compare_button > 0) {
      ptab <- statmod_results()  # this runs compareGrowthCurves internally
      dat <- filtered_data() |> 
        filter(Sample_ID %in% input$samples) |> 
        mutate(replicate = paste(Sample_ID, Row, Column, sep = "_")) |> 
        filter(Sample_ID != "Empty")
      
      if (!is.null(input$max_timepoint) && !is.na(input$max_timepoint)) {
        dat <- dat %>% filter(timepoint <= input$max_timepoint)
      }
      
      wide <- dat |>
        select(replicate, timepoint, Absorbance) |>
        pivot_wider(names_from = timepoint, values_from = Absorbance)
      
      y_mat <- as.matrix(wide |> select(-replicate))
      ord <- order(as.numeric(colnames(y_mat)))
      y_mat <- y_mat[, ord]
      colnames(y_mat) <- colnames(y_mat)
      
      grp <- factor(stringr::str_remove(wide$replicate, "_[A-H]_[0-9]+$"))
      keep_levels <- names(which(table(grp) >= 2))
      y_mat <- y_mat[grp %in% keep_levels, , drop = FALSE]
      grp <- droplevels(grp[grp %in% keep_levels])
      
      # Create Rep IDs (if rownames not present)
      rep_ids <- paste0("Rep", seq_len(nrow(y_mat)))
      timepoints <- as.numeric(colnames(y_mat))
      
      gg_df <- as.data.frame(y_mat)
      gg_df$Group <- grp
      gg_df$Replicate <- rep_ids
      
      gg_long <- gg_df %>%
        pivot_longer(cols = -c(Group, Replicate), names_to = "Time", values_to = "Response") %>%
        mutate(Time = as.numeric(Time))
      
      plot <- ggplot(gg_long, aes(x = Time, y = Response, group = Replicate, color = Group)) +
        geom_line(alpha = 0.6) +
        labs(
          title = "Growth Curves by Group",
          x = "Time",
          y = "Response"
        ) +
        theme_minimal()
      
      return(plot)
    }    
    # Standard ggplot-based plotting
    data <- plot_data()
    
    if (input$toggle_average) {
      if (input$toggle_mean_only) {
        plot <- ggplot(data, aes(x = timepoint, y = mean_absorbance, color = Sample_ID)) +
          geom_point(size = input$point_size) +
          theme_minimal() +
          labs(
            title = "Absorbance over Time",
            x = "Timepoint",
            y = "Mean Absorbance",
            color = "Sample ID"
          )
      } else {
        plot <- ggplot(data, aes(x = timepoint, y = mean_absorbance, color = Sample_ID)) +
          geom_line(size = 1) +
          geom_point(size = input$point_size) +
          geom_errorbar(aes(ymin = mean_absorbance - se_absorbance, ymax = mean_absorbance + se_absorbance), width = 0.5) +
          theme_minimal() +
          labs(
            title = "Absorbance over Time",
            x = "Timepoint",
            y = "Mean Absorbance",
            color = "Sample ID"
          ) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right"
          )
      }
    } else {
      plot <- ggplot(data, aes(x = timepoint, y = Absorbance, color = Sample_ID, group = interaction(Row, Column))) +
        geom_line(alpha = 0.8) +
        geom_point(size = input$point_size, alpha = 0.8) +
        theme_minimal() +
        labs(
          title = "Absorbance over Time",
          x = "Timepoint",
          y = "Absorbance",
          color = "Sample ID"
        ) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right"
        )
    }
    
    plot <- plot + coord_cartesian(
      xlim = c(input$x_min, input$x_max),
      ylim = c(input$y_min, input$y_max)
    )
    
    return(plot)
  })  
  output$absorbance_plot <- renderPlotly({
    
    ggplotly(reactive_plot())
  })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("absorbance_plot", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = input$pdfWidth, height = input$pdfHeight, useDingbats = FALSE)
      print(
        reactive_plot()
      )
      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)
