# app.R
library(shiny)
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(plotly)
library(stringr)
library(statmod)
library(nlme)
library(splines)
library(broom)
library(purrr)
library(scales)

ui <- navbarPage(
  titlePanel("CLARIOstar Absorbance Visualization (Linear / Log visual only)"),
  tabPanel("Plotting",
           sidebarLayout(
             sidebarPanel(
               fileInput("file", "Upload Excel File:", accept = ".xlsx"),
               fileInput("file2", "Upload Additonal Excel File (repeatable):", accept = ".xlsx"),
               checkboxInput("use_blank_corrected", "Use Blank Corrected Values", value = TRUE),
               
               radioButtons("axis_scale", "Y-axis scale (visual only):",
                            choices = c("Linear" = "linear", "Log10 axis" = "log10", "Log2 axis" = "log2"),
                            selected = "linear"),
               helpText("Axis scale affects only the axis (data are NOT transformed)."),
               
               hr(),
               actionButton("uncheck_all", "Uncheck All Samples"),
               checkboxGroupInput("samples", "Select Samples to Include:", choices = NULL, selected = NULL),
               
               hr(),
               selectInput("stat_method", "Select Statistical Method:",
                           choices = list(
                             "Statmod permutation trajectory (original)" = "statmod",
                             "A) Linear mixed model (lme)" = "lme",
                             "B) Mixed model with splines/polynomials (lme_spline)" = "lme_spline",
                             "C) Repeated measures ANOVA (aov within-subject)" = "rm_anova",
                             "D) Per-timepoint tests (t-test / Wilcox) with correction" = "per_timepoint"
                           ), selected = "statmod"),
               conditionalPanel(
                 condition = "input.stat_method == 'lme' || input.stat_method == 'lme_spline'",
                 checkboxInput("lme_random_intercept_slope", "Include random slope for time (random intercept only otherwise)", value = FALSE)
               ),
               conditionalPanel(
                 condition = "input.stat_method == 'lme_spline'",
                 selectInput("spline_or_poly", "Spline or polynomial:", choices = c("natural_spline" = "ns", "polynomial" = "poly"), selected = "ns"),
                 conditionalPanel(condition = "input.spline_or_poly == 'poly'",
                                  numericInput("poly_degree", "Polynomial degree:", value = 2, min = 1, max = 5))
               ),
               conditionalPanel(
                 condition = "input.stat_method == 'per_timepoint'",
                 radioButtons("per_time_test", "Per-timepoint test:", choices = c("t.test" = "ttest", "wilcox.test" = "wilcox"), selected = "ttest"),
                 selectInput("per_time_adj", "Multiplicity correction:", choices = c("BH"="BH","holm"="holm","bonferroni"="bonferroni"), selected = "BH")
               ),
               
               hr(),
               numericInput("nsim_input", "Number of Permutations for Growth Curve Comparison (nsim)",
                            value = 200, min = 10, max = 5000),
               helpText("Statmod permutation method uses nsim; other methods ignore nsim."),
               actionButton("compare_button", "Run Pairwise Comparison")
             ),
             
             mainPanel(
               # ensure width="100%" so Shiny reports container width; the server will pass explicit pixel width to plotly
               plotlyOutput("absorbance_plot", height = "700px", width = "100%"),
               checkboxInput("toggle_average", "Show Average with Error Bars", TRUE),
               checkboxInput("toggle_mean_only", "Show Only Mean Values", FALSE),
               checkboxInput("toggle_groupcurves", "Show Growth Curve Group Means (Run Pairwise Comparison First)", FALSE),
               numericInput("max_timepoint", "Max Timepoint to Include:", value = NA, min = 0),
               actionButton("plot_button", "Generate Plot"),
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
                        actionButton("apply_ranges", "Apply Axis Ranges"),  # <- NEW button
                        sliderInput("point_size", "Point Size:", min = 0.25, max = 5, value = 1, step = 0.25)
                 ),
                 column(6,
                        h3("Pairwise Results"),
                        h5("Table will show pairwise comparisons and p-values depending on selected method"),
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
               "Using a large number will be more accurate but be much slower. Both growthcurve from ",
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
# compareGrowthCurves (unchanged) ...
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
  uploaded_files <- reactiveValues(extra_files = list(), counter = 2)
  
  observeEvent(input$uncheck_all, {
    updateCheckboxGroupInput(session, "samples", selected = character(0))
  })
  
  ## ---------- File processing ----------
  process_file <- function(file_path) {
    data <- readxl::read_excel(file_path, sheet = 1)
    data <- data[, 1:13]
    col2 <- as.character(unlist(data[[2]]))
    
    if (input$use_blank_corrected) {
      absorbance_start <- which(grepl("Blank corrected based on Raw Data", col2, ignore.case = TRUE))
    } else {
      absorbance_start <- which(grepl("Raw Data", col2, ignore.case = TRUE) &
                                  !grepl("blank", col2, ignore.case = TRUE))
    }
    
    sample_ids_start <- grep("Sample IDs", col2, ignore.case = TRUE)
    
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
    
    existing_ids <- c()
    if (!is.null(input$file)) {
      existing_ids <- unique(process_file(input$file$datapath)$Sample_ID)
    }
    existing_ids <- unique(c(existing_ids,
                             unlist(lapply(uploaded_files$extra_files, function(df) unique(df$Sample_ID)))))
    
    dup_ids <- intersect(existing_ids, unique(data2$Sample_ID))
    for (id in dup_ids) {
      data2$Sample_ID[data2$Sample_ID == id] <- paste0(id, "_", count)
    }
    uploaded_files$extra_files[[paste0("file2_", count)]] <- data2
    uploaded_files$counter <- count + 1
  })
  
  # ---------- processed_data reactive (NO transforms) ----------
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
    
    # KEEP RAW values only; do not apply any log transforms.
    combined
  })
  # -----------------------------------------------------------------------
  
  output$sampleRemover <- renderUI({
    req(processed_data())
    toSearch <- unique(processed_data()$uniqueSampleID)
    selectInput("samplesFilter", "Select SampleIDs to remove:", choices = toSearch, multiple = TRUE)
  })
  
  selectedSamples <- reactiveVal(NULL)
  observeEvent(input$applysampleRemover, {
    selectedSamples(input$samplesFilter)
  })
  
  # Update checkboxes choices, preserving previous explicit selection (including explicit all-unchecked).
  observeEvent(processed_data(), {
    pd <- processed_data()
    vals <- unique(pd$Sample_ID)
    
    prefix   <- str_extract(vals, "^[^_]+")
    treatment <- str_remove(vals, "^[^_]+_")
    compound <- str_extract(treatment, "^[A-Za-z]+")
    number   <- as.numeric(str_extract(treatment, "\\d+\\.?\\d*"))
    number[is.na(number)] <- Inf
    ord <- order(prefix, compound, number)
    sorted_vals <- vals[ord]
    
    prev <- isolate(input$samples)
    if (is.null(prev)) {
      selected_vals <- sorted_vals
    } else {
      selected_vals <- intersect(sorted_vals, prev)
    }
    
    updateCheckboxGroupInput(session, "samples", choices = sorted_vals, selected = selected_vals)
  }, ignoreNULL = FALSE, ignoreInit = FALSE)
  
  filtered_data <- reactive({
    req(processed_data())
    data <- processed_data()
    
    sampleRemove <- selectedSamples()
    if (!is.null(sampleRemove) && length(sampleRemove) > 0) {
      data <- data[!data$uniqueSampleID %in% sampleRemove, , drop = FALSE]
    }
    
    if (!is.null(input$samples)) {
      data <- data %>% filter(Sample_ID %in% input$samples)
    }
    
    if (!is.null(input$max_timepoint) && !is.na(input$max_timepoint)) {
      data <- data %>% filter(timepoint <= input$max_timepoint)
    }
    
    data <- data %>% filter(Sample_ID != "Empty")
    data
  })
  
  # ---------- axis ranges (linear) ----------
  display_axis_ranges <- reactive({
    req(filtered_data())
    data <- filtered_data()
    if (nrow(data) == 0) return(list(x_min = 0, x_max = 1, y_min = 0, y_max = 1))
    
    x_min <- min(data$timepoint, na.rm = TRUE)
    x_max <- max(data$timepoint, na.rm = TRUE)
    
    if (isTRUE(input$toggle_average)) {
      sd <- data %>%
        group_by(Sample_ID, timepoint) %>%
        summarise(
          mean_absorbance = mean(Absorbance, na.rm = TRUE),
          se_absorbance   = sd(Absorbance, na.rm = TRUE) / sqrt(n()),
          n = n(),
          .groups = "drop"
        )
      y_vec <- sd$mean_absorbance
    } else {
      y_vec <- data$Absorbance
    }
    
    y_vec <- y_vec[is.finite(y_vec) & !is.na(y_vec)]
    if (length(y_vec) == 0) {
      return(list(x_min = x_min, x_max = x_max, y_min = 0, y_max = 1))
    }
    
    y_min <- min(y_vec, na.rm = TRUE)
    y_max <- max(y_vec, na.rm = TRUE)
    
    pad <- (y_max - y_min) * 0.06
    if (!is.finite(pad)) pad <- 0.02
    y_min <- max(0, y_min - pad)
    y_max <- y_max + pad
    
    list(x_min = x_min, x_max = x_max, y_min = y_min, y_max = y_max)
  })
  
  output$x_range_inputs <- renderUI({
    req(display_axis_ranges())
    r <- display_axis_ranges()
    tagList(
      numericInput("x_min", "X-Axis Minimum:", value = r$x_min),
      numericInput("x_max", "X-Axis Maximum:", value = r$x_max)
    )
  })
  
  output$y_range_inputs <- renderUI({
    req(display_axis_ranges())
    r <- display_axis_ranges()
    tagList(
      numericInput("y_min", "Y-Axis Minimum:", value = r$y_min),
      numericInput("y_max", "Y-Axis Maximum:", value = r$y_max)
    )
  })
  
  # --- Track whether the user has explicitly applied ranges ---
  ranges_applied <- reactiveVal(FALSE)
  
  # Whenever data or toggle_average changes we programmatically update numericInputs
  # and reset ranges_applied() to FALSE (because those updates are NOT user intent).
  observeEvent(list(input$toggle_average, filtered_data()), {
    r <- display_axis_ranges()
    updateNumericInput(session, "y_min", value = r$y_min)
    updateNumericInput(session, "y_max", value = r$y_max)
    updateNumericInput(session, "x_min", value = r$x_min)
    updateNumericInput(session, "x_max", value = r$x_max)
    ranges_applied(FALSE)   # programmatic default ? user hasn't applied anything yet
  }, ignoreInit = FALSE)
  
  # When user clicks the Apply Axis Ranges button, mark ranges_applied TRUE.
  observeEvent(input$apply_ranges, {
    ranges_applied(TRUE)
    showNotification("Axis ranges applied", type = "message", duration = 2)
  })
  
  plot_data <- reactive({
    req(filtered_data())
    data <- filtered_data()
    if (input$toggle_average) {
      data <- data %>%
        group_by(Sample_ID, timepoint) %>%
        summarise(
          mean_absorbance = mean(Absorbance, na.rm = TRUE),
          se_absorbance   = sd(Absorbance, na.rm = TRUE) / sqrt(n()),
          n = n(),
          .groups = "drop"
        )
      data
    } else {
      data
    }
  })
  
  # helper to compute log ticks (majors & minors) given positive numeric range
  compute_log_ticks <- function(minv, maxv, base = 10) {
    # enforce strictly positive values
    minv <- max(minv, .Machine$double.eps)
    maxv <- max(maxv, minv * 1.001)
    
    if (base == 10) {
      # ---- LOG10 ----
      # clamp decades to 10^-3 .. 10^1 (0.001 .. 10)
      dmin <- floor(log10(minv))
      dmax <- ceiling(log10(maxv))
      dmin <- max(dmin, -3)
      dmax <- min(dmax, 1)
      
      majors <- 10 ^ seq(dmin, dmax)
      
      # classic decade minors: 2..9 × 10^d
      minors <- unlist(lapply(seq(dmin, dmax), function(d) (2:9) * 10^d))
      
      majors <- majors[majors >= 10^dmin & majors <= 10^dmax]
      minors <- minors[minors >= 10^dmin & minors <= 10^dmax]
      
      list(
        majors = as.numeric(sort(majors)),
        minors = as.numeric(sort(minors))
      )
      
    } else if (base == 2) {
      # ---- LOG2 ----
      # clamp to sensible powers of two
      # 2^-5 = 0.03125   2^4 = 16
      dmin <- floor(log(minv, base = 2))
      dmax <- ceiling(log(maxv, base = 2))
      dmin <- max(dmin, -5)
      dmax <- min(dmax, 4)
      
      # major ticks = powers of two
      majors <- 2 ^ seq(dmin, dmax)
      
      # minor ticks = halfway between powers of two (1.5 * 2^d)
      minors <- unlist(lapply(seq(dmin, dmax - 1), function(d) 1.5 * 2^d))
      
      majors <- majors[majors >= 2^dmin & majors <= 2^dmax]
      minors <- minors[minors >= 2^dmin & minors <= 2^dmax]
      
      list(
        majors = as.numeric(sort(majors)),
        minors = as.numeric(sort(minors))
      )
      
    } else {
      # fallback
      list(
        majors = c(minv, maxv),
        minors = numeric(0)
      )
    }
  }
  
  # ---------------- Plotly interactive plot (FORCED SAME PLOT-AREA PIXEL WIDTH) ----------------
  output$absorbance_plot <- renderPlotly({
    req(plot_data(), input$plot_button)
    pd <- plot_data()
    
    defaults <- display_axis_ranges()
    get_num <- function(x) {
      if (is.null(x)) return(NA_real_)
      v <- suppressWarnings(as.numeric(x))
      if (is.na(v) || !is.finite(v)) return(NA_real_)
      v
    }
    
    if (input$toggle_average) {
      sd <- pd %>% arrange(Sample_ID, timepoint)
      if (nrow(sd) == 0) return(NULL)
      sd$plot_y <- sd$mean_absorbance
      sd$plot_ymin <- sd$mean_absorbance - sd$se_absorbance
      sd$plot_ymax <- sd$mean_absorbance + sd$se_absorbance
    } else {
      sd <- pd %>% arrange(Sample_ID, Row, Column, timepoint)
      if (nrow(sd) == 0) return(NULL)
      sd$plot_y <- sd$Absorbance
    }
    
    # Note: only consider user-specified ranges if they've clicked Apply Axis Ranges.
    ui_ymin <- get_num(input$y_min)
    ui_ymax <- get_num(input$y_max)
    ui_xmin <- get_num(input$x_min)
    ui_xmax <- get_num(input$x_max)
    ui_ymin <- ifelse(input$axis_scale %in% c("log10", "log2") &&
                        is.finite(ui_ymin) && ui_ymin <= 0,
                      .Machine$double.eps,
                      ui_ymin)
    
    ui_specified <- ranges_applied() && is.finite(ui_ymin) && is.finite(ui_ymax) && (ui_ymax > ui_ymin)
    ui_x_specified <- ranges_applied() && is.finite(ui_xmin) && is.finite(ui_xmax) && (ui_xmax > ui_xmin)
    
    if (ui_specified) {
      use_ymin <- ui_ymin; use_ymax <- ui_ymax
    } else {
      use_ymin <- defaults$y_min; use_ymax <- defaults$y_max
    }
    
    if (ui_x_specified) {
      use_xmin <- ui_xmin; use_xmax <- ui_xmax
    } else {
      use_xmin <- defaults$x_min; use_xmax <- defaults$x_max
    }
    
    # get container width (pixels) from the client
    plot_w <- session$clientData$output_absorbance_plot_width
    if (is.null(plot_w) || !is.finite(plot_w) || plot_w <= 0) plot_w <- 900
    plot_h <- 760
    
    # fixed left margin (pixels) to accommodate y-axis labels consistently
    left_margin_px <- 120
    top_margin_px <- 110
    bottom_margin_px <- 80
    
    # Reserve a fixed pixel width for the legend area (to the right).
    # We'll compute a desired plot-area width in pixels and derive domain from that so it is identical
    # across axis types. The legend area width must be >= 140 px to show labels.
    reserved_legend_px <- 220  # space guaranteed for legend (adjust if you want more)
    
    # ensure reserved_legend_px isn't larger than available width
    reserved_legend_px <- min(reserved_legend_px, max(140, round(plot_w * 0.4)))
    
    # desired plot area width in pixels: use remaining width after left margin + reserved legend,
    # but limit to a sensible minimum/maximum so plot_area_px is stable.
    avail_for_plot <- plot_w - left_margin_px - reserved_legend_px
    # enforce minimum plot area (so we don't get too small)
    desired_plot_area_px <- max(420, round(avail_for_plot))
    # if container is narrow, reduce reserved legend to fit but keep plot area >= 420
    if (desired_plot_area_px + left_margin_px + reserved_legend_px > plot_w) {
      # shrink reserved_legend_px first
      reserved_legend_px <- max(140, plot_w - left_margin_px - 420)
      desired_plot_area_px <- max(420, plot_w - left_margin_px - reserved_legend_px)
    }
    
    # Compute domain for xaxis so the plot area occupies exactly desired_plot_area_px pixels
    domain_left <- left_margin_px / plot_w
    domain_right <- (left_margin_px + desired_plot_area_px) / plot_w
    domain_left <- max(0, min(domain_left, 0.9))
    domain_right <- max(domain_left + 0.2, min(domain_right, 0.99))
    
    # recompute right margin so total margins add up to container width
    # margin.r = pixels on the right (so that left + plot area + right = plot_w)
    right_margin_px <- max(10, plot_w - left_margin_px - desired_plot_area_px)
    
    # Start plotly with explicit width/height - this keeps pixel width stable
    p <- plot_ly(width = plot_w, height = plot_h)
    
    # add traces (unchanged)
    if (input$toggle_average) {
      for (g in unique(sd$Sample_ID)) {
        sg <- sd %>% filter(Sample_ID == g) %>% arrange(timepoint)
        if (nrow(sg) == 0) next
        if (!input$toggle_mean_only) {
          p <- p %>% add_ribbons(data = sg, x = ~timepoint, ymin = ~plot_ymin, ymax = ~plot_ymax,
                                 name = paste0(g, " CI"), inherit = FALSE, showlegend = FALSE, hoverinfo = "none", opacity = 0.22)
        }
        p <- p %>% add_lines(data = sg, x = ~timepoint, y = ~plot_y,
                             name = ~Sample_ID, inherit = FALSE, showlegend = TRUE,
                             hoverinfo = "text",
                             text = ~paste0("Sample: ", Sample_ID, "<br>Time: ", timepoint, "<br>Value: ", signif(plot_y, 4)))
        p <- p %>% add_markers(data = sg, x = ~timepoint, y = ~plot_y,
                               name = ~Sample_ID, inherit = FALSE, showlegend = FALSE,
                               marker = list(size = input$point_size * 3),
                               hoverinfo = "text",
                               text = ~paste0("Sample: ", Sample_ID, "<br>Time: ", timepoint, "<br>Value: ", signif(plot_y, 4)))
      }
      ytitle <- "Mean Absorbance"
    } else {
      reps <- unique(paste(sd$Sample_ID, sd$Row, sd$Column, sep = "_"))
      for (rep in reps) {
        rdf <- sd %>% filter(paste(Sample_ID, Row, Column, sep = "_") == rep)
        if (nrow(rdf) == 0) next
        lab <- unique(rdf$Sample_ID)
        p <- p %>% add_lines(data = rdf, x = ~timepoint, y = ~plot_y,
                             inherit = FALSE, name = paste0(lab, " (rep)"),
                             line = list(width = 0.8), showlegend = FALSE, opacity = 0.35,
                             hoverinfo = "text",
                             text = ~paste0("Sample: ", Sample_ID, "<br>Row: ", Row, "<br>Col: ", Column,
                                            "<br>Time: ", timepoint, "<br>Abs: ", signif(plot_y,4)))
      }
      p <- p %>% add_markers(data = sd, x = ~timepoint, y = ~plot_y,
                             inherit = FALSE, showlegend = FALSE,
                             marker = list(size = input$point_size * 2),
                             hoverinfo = "text",
                             text = ~paste0("Sample: ", Sample_ID, "<br>Row: ", Row, "<br>Col: ", Column,
                                            "<br>Time: ", timepoint, "<br>Abs: ", signif(plot_y,4)))
      ytitle <- "Absorbance"
    }
    
    axis_scale_choice <- input$axis_scale
    
    # Turn OFF automargin to prevent Plotly from changing domain to fit tick labels.
    # We'll provide margins explicitly so plot area in pixels stays identical.
    yaxis_list <- list(title = ytitle, showgrid = TRUE, zeroline = FALSE,
                       ticks = "outside", ticklen = 8, tickwidth = 1.2,
                       automargin = FALSE)
    
    if (axis_scale_choice == "linear") {
      yaxis_list <- modifyList(yaxis_list, list(range = c(use_ymin, use_ymax), autorange = FALSE))
    } else if (axis_scale_choice == "log10") {
      pos_min <- max(use_ymin, .Machine$double.eps)
      pos_max <- max(use_ymax, pos_min * 1.001)
      
      # DEFAULT start for log plots when user didn't set y_min: 0.04
      if (ui_specified) {
        pos_min_clamped <- max(pos_min, .Machine$double.eps)
        pos_max_clamped <- max(pos_max, pos_min_clamped * 1.001)
        axis_range_log <- c(log10(use_ymin), log10(use_ymax))
      } else {
        pos_min_clamped <- 0.04
        pos_max_clamped <- min(pos_max, 10)
        if (pos_max_clamped <= pos_min_clamped * 1.001) pos_max_clamped <- pos_min_clamped * 10
        axis_range_log <- c(log10(pos_min_clamped), log10(pos_max_clamped))
      }
      
      ticks <- compute_log_ticks(pos_min_clamped, pos_max_clamped, base = 10)
      tickvals <- sort(unique(c(ticks$majors, ticks$minors)))
      ticktext <- ifelse(tickvals %in% ticks$majors, format(tickvals, scientific = FALSE, trim = TRUE), "")
      
      yaxis_list <- modifyList(yaxis_list, list(
        type = "log",
        tickmode = "array",
        tickvals = tickvals,
        ticktext = ticktext,
        autorange = FALSE,
        range = axis_range_log
      ))
      
      zero_annotation <- list(xref = "paper", x = 0, xanchor = "right",
                              yref = "paper", y = 0, yanchor = "top",
                              text = "0", showarrow = FALSE,
                              font = list(size = 10, color = "black"),
                              bgcolor = "rgba(255,255,255,0.0)")
    } else if (axis_scale_choice == "log2") {
      pos_min <- max(use_ymin, .Machine$double.eps)
      pos_max <- max(use_ymax, pos_min * 1.001)
      
      if (ui_specified) {
        pos_min_clamped <- max(pos_min, .Machine$double.eps)
        pos_max_clamped <- max(pos_max, pos_min_clamped * 1.001)
        axis_range_log <- c(log10(use_ymin), log10(use_ymax))
      } else {
        pos_min_clamped <- 0.04
        pos_max_clamped <- min(pos_max, 16)
        if (pos_max_clamped <= pos_min_clamped * 1.001) pos_max_clamped <- pos_min_clamped * 8
        axis_range_log <- c(log10(pos_min_clamped), log10(pos_max_clamped))
      }
      
      ticks <- compute_log_ticks(pos_min_clamped, pos_max_clamped, base = 2)
      majors <- ticks$majors
      minors <- ticks$minors
      
      tickvals <- sort(unique(c(majors, minors)))
      ticktext <- ifelse(tickvals %in% majors, as.character(tickvals), "")
      
      yaxis_list <- modifyList(yaxis_list, list(
        type = "log",
        tickmode = "array",
        tickvals = tickvals,
        ticktext = ticktext,
        autorange = FALSE,
        range = axis_range_log
      ))
      
      zero_annotation <- list(xref = "paper", x = 0, xanchor = "right",
                              yref = "paper", y = 0, yanchor = "top",
                              text = "0", showarrow = FALSE,
                              font = list(size = 10, color = "black"),
                              bgcolor = "rgba(255,255,255,0.0)")
    }
    
    # x-axis domain forced to pixel-precise region computed above
    xaxis_layout <- list(
      title = "Timepoint",
      range = if (ui_x_specified) c(use_xmin, use_xmax) else c(defaults$x_min, defaults$x_max),
      domain = c(domain_left, domain_right),
      automargin = FALSE
    )
    
    legend_layout <- list(
      orientation = "v",
      x = domain_right + 0.02,
      xanchor = "left",
      y = 0.5,
      yanchor = "middle",
      bgcolor = "rgba(255,255,255,0.0)",
      borderwidth = 0
    )
    
    # Final layout: use explicit width/height and margins (left/right in px)
    p <- p %>% layout(title = if (input$toggle_average) "Absorbance over Time (summary)" else "Absorbance over Time (raw replicates)",
                      xaxis = xaxis_layout,
                      yaxis = yaxis_list,
                      margin = list(l = left_margin_px, r = right_margin_px, b = bottom_margin_px, t = top_margin_px),
                      legend = legend_layout,
                      autosize = FALSE,
                      width = plot_w,
                      height = plot_h)
    
    if (exists("zero_annotation")) p <- p %>% layout(annotations = list(zero_annotation))
    
    # Ensure Plotly doesn't try to be responsive and change the sizing after render
    p <- p %>% config(responsive = FALSE)
    
    p
  })
  # ---------------------------------------------------------------------------
  
  # ---------------- PDF ggplot (respects axis_scale visual choice) ----------------
  reactive_plot <- reactive({
    req(plot_data(), input$plot_button)
    data <- plot_data()
    
    common_theme <- theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
    
    if (input$toggle_average) {
      data <- data %>% arrange(Sample_ID, timepoint)
      data$plot_y <- data$mean_absorbance
      data$plot_ymin <- data$mean_absorbance - data$se_absorbance
      data$plot_ymax <- data$mean_absorbance + data$se_absorbance
      
      ytitle <- "Mean Absorbance"
      if (input$toggle_mean_only) {
        p <- ggplot(data, aes(x = timepoint, y = plot_y, color = Sample_ID)) +
          geom_point(size = input$point_size) + common_theme +
          labs(title = "Absorbance over Time", x = "Timepoint", y = ytitle, color = "Sample ID")
      } else {
        p <- ggplot(data, aes(x = timepoint, y = plot_y, color = Sample_ID)) +
          geom_line(size = 1) +
          geom_point(size = input$point_size) +
          geom_errorbar(aes(ymin = plot_ymin, ymax = plot_ymax), width = 0.5) +
          common_theme +
          labs(title = "Absorbance over Time", x = "Timepoint", y = ytitle, color = "Sample ID")
      }
      
    } else {
      data$plot_y <- data$Absorbance
      p <- ggplot(data, aes(x = timepoint, y = plot_y, color = Sample_ID, group = interaction(Row, Column))) +
        geom_line(alpha = 0.8) +
        geom_point(size = input$point_size, alpha = 0.8) +
        common_theme +
        labs(title = "Absorbance over Time", x = "Timepoint", y = "Absorbance", color = "Sample ID")
    }
    
    defaults <- display_axis_ranges()
    get_num <- function(x) {
      if (is.null(x)) return(NA_real_)
      v <- suppressWarnings(as.numeric(x))
      if (is.na(v) || !is.finite(v)) return(NA_real_)
      v
    }
    ui_ymin <- get_num(input$y_min)
    ui_ymax <- get_num(input$y_max)
    ui_xmin <- get_num(input$x_min)
    ui_xmax <- get_num(input$x_max)
    ui_ymin <- ifelse(input$axis_scale %in% c("log10", "log2") &&
                        is.finite(ui_ymin) && ui_ymin <= 0,
                      .Machine$double.eps,
                      ui_ymin)
    
    ui_specified <- ranges_applied() && is.finite(ui_ymin) && is.finite(ui_ymax) && (ui_ymax > ui_ymin)
    ui_x_specified <- ranges_applied() && is.finite(ui_xmin) && is.finite(ui_xmax) && (ui_xmax > ui_xmin)
    
    use_xmin <- if (ui_x_specified) ui_xmin else defaults$x_min
    use_xmax <- if (ui_x_specified) ui_xmax else defaults$x_max
    use_ymin <- if (ui_specified) ui_ymin else defaults$y_min
    use_ymax <- if (ui_specified) ui_ymax else defaults$y_max
    
    axis_scale_choice <- input$axis_scale
    
    if (axis_scale_choice == "linear") {
      p <- p + coord_cartesian(xlim = c(use_xmin, use_xmax), ylim = c(use_ymin, use_ymax))
    } else if (axis_scale_choice == "log10") {
      pos_min <- max(use_ymin, .Machine$double.eps)
      pos_max <- max(use_ymax, pos_min * 1.001)
      
      if (ui_specified) {
        pos_min_clamped <- max(pos_min, .Machine$double.eps)
        pos_max_clamped <- max(pos_max, pos_min_clamped * 1.001)
        limits_vals <- c(use_ymin, use_ymax)
      } else {
        pos_min_clamped <- 0.04
        pos_max_clamped <- min(pos_max, 10)
        if (pos_max_clamped <= pos_min_clamped * 1.001) pos_max_clamped <- pos_min_clamped * 10
        limits_vals <- c(pos_min_clamped, pos_max_clamped)
      }
      
      ticks <- compute_log_ticks(pos_min_clamped, pos_max_clamped, base = 10)
      majors <- ticks$majors; minors <- ticks$minors
      
      p <- p + scale_y_log10(limits = limits_vals,
                             breaks = majors,
                             minor_breaks = minors) +
        coord_cartesian(xlim = c(use_xmin, use_xmax))
    } else if (axis_scale_choice == "log2") {
      pos_min <- max(use_ymin, .Machine$double.eps)
      pos_max <- max(use_ymax, pos_min * 1.001)
      
      if (ui_specified) {
        pos_min_clamped <- max(pos_min, .Machine$double.eps)
        pos_max_clamped <- max(pos_max, pos_min_clamped * 1.001)
        limits_vals <- c(use_ymin, use_ymax)
      } else {
        pos_min_clamped <- 0.04
        pos_max_clamped <- min(pos_max, 16)
        if (pos_max_clamped <= pos_min_clamped * 1.001) pos_max_clamped <- pos_min_clamped * 8
        limits_vals <- c(pos_min_clamped, pos_max_clamped)
      }
      
      ticks <- compute_log_ticks(pos_min_clamped, pos_max_clamped, base = 2)
      majors <- ticks$majors
      minors <- ticks$minors
      
      maj_labels <- ifelse(majors < 1,
                           formatC(majors, format = "f", digits = 3),
                           formatC(majors, format = "f", digits = 0))
      maj_labels <- sub("0+\\.$", "0", maj_labels)
      
      p <- p +
        scale_y_continuous(
          trans = log2_trans(),
          limits = limits_vals,
          breaks = majors,
          minor_breaks = minors,
          labels = maj_labels
        ) +
        coord_cartesian(xlim = c(use_xmin, use_xmax)) +
        theme(
          panel.grid.minor = element_line(color = "gray90", size = 0.25),
          panel.grid.major = element_line(color = "gray85"),
          axis.ticks.length = unit(4, "pt")
        )
    }
    
    
    p
  })
  
  # Force-evaluate reactive_plot and log to console + UI when 'Generate Plot' clicked
  observeEvent(input$plot_button, {
    # small UI feedback
    showNotification("Plot button clicked: computing ticks...", type = "message", duration = 2)
    
    # Force evaluation of reactive_plot (keeps behavior you already added)
    p <- reactive_plot()
    
    # compute the same use_ymin/use_ymax logic that your plotting code uses
    defaults <- display_axis_ranges()
    get_num <- function(x) {
      if (is.null(x)) return(NA_real_)
      v <- suppressWarnings(as.numeric(x))
      if (is.na(v) || !is.finite(v)) return(NA_real_)
      v
    }
    ui_ymin <- get_num(input$y_min)
    ui_ymax <- get_num(input$y_max)
    ui_specified <- ranges_applied() && is.finite(ui_ymin) && is.finite(ui_ymax) && (ui_ymax > ui_ymin)
    
    if (ui_specified) {
      use_ymin <- ui_ymin
      use_ymax <- ui_ymax
    } else {
      use_ymin <- defaults$y_min
      use_ymax <- defaults$y_max
    }
    
    axis_scale_choice <- isolate(input$axis_scale)
    
    cat("=== Tick debug ===\n")
    cat("Axis choice:", axis_scale_choice, "\n")
    cat("Computed use_ymin:", use_ymin, " use_ymax:", use_ymax, "\n")
    
    if (axis_scale_choice == "linear") {
      pretty_ticks <- pretty(c(use_ymin, use_ymax), n = 5)
      cat("Linear major ticks (pretty):", paste(pretty_ticks, collapse = ", "), "\n")
      showNotification(paste0("Linear ticks: ", paste(pretty_ticks, collapse = ", ")), type = "message", duration = 4)
      
    } else if (axis_scale_choice == "log10") {
      pos_min <- max(use_ymin, .Machine$double.eps)
      pos_max <- max(use_ymax, pos_min * 1.001)
      
      if (!ui_specified) {
        pos_min_clamped <- max(pos_min, 0.04)
      } else {
        pos_min_clamped <- max(pos_min, .Machine$double.eps)
      }
      pos_max_clamped <- min(pos_max, 10)
      if (pos_max_clamped <= pos_min_clamped * 1.001) pos_max_clamped <- pos_min_clamped * 10
      
      ticks <- compute_log_ticks(pos_min_clamped, pos_max_clamped, base = 10)
      majors <- ticks$majors
      minors <- ticks$minors
      
      cat("Log10 majors:", paste(majors, collapse = ", "), "\n")
      cat("Log10 minors (subset):", if (length(minors)>0) paste(head(minors, 20), collapse = ", ") else "(none)", "\n")
      showNotification(paste0("Log10 majors: ", paste(majors, collapse = ", ")), type = "message", duration = 5)
      
    } else if (axis_scale_choice == "log2") {
      pos_min <- max(use_ymin, .Machine$double.eps)
      pos_max <- max(use_ymax, pos_min * 1.001)
      
      if (!ui_specified) {
        pos_min_clamped <- max(pos_min, 0.04)
      } else {
        pos_min_clamped <- max(pos_min, .Machine$double.eps)
      }
      pos_max_clamped <- min(pos_max, 16)
      if (pos_max_clamped <= pos_min_clamped * 1.001) pos_max_clamped <- pos_min_clamped * 8
      
      ticks <- compute_log_ticks(pos_min_clamped, pos_max_clamped, base = 2)
      majors <- ticks$majors
      minors <- ticks$minors
      
      cat("Log2 pos_min_clamped:", pos_min_clamped, " pos_max_clamped:", pos_max_clamped, "\n")
      cat("Log2 majors (powers of 2):", paste(majors, collapse = ", "), "\n")
      if (length(minors)) cat("Log2 minors:", paste(minors, collapse = ", "), "\n") else cat("Log2 minors: (none)\n")
      
      showNotification(paste0("Log2 majors: ", paste(majors, collapse = ", ")), type = "message", duration = 6)
    }
    
    if (inherits(p, "ggplot")) {
      cat("reactive_plot is a ggplot: layers=", length(p$layers), " scales=", length(p$scales$scales), "\n")
      if (length(p$scales$scales) > 0) {
        cat("ggplot scale objects:\n")
        for (i in seq_along(p$scales$scales)) {
          sc <- p$scales$scales[[i]]
          cat(" - scale", i, "class:", paste(class(sc), collapse = ","), "\n")
        }
      }
    } else {
      cat("reactive_plot not a ggplot (class):", paste(class(p), collapse = ", "), "\n")
    }
    
    flush.console()
  })
  
  
  output$download_plot <- downloadHandler(
    filename = function() paste0("absorbance_plot_", Sys.Date(), ".pdf"),
    content = function(file) {
      pdf(file, width = input$pdfWidth, height = input$pdfHeight, useDingbats = FALSE)
      print(reactive_plot())
      dev.off()
    }
  )
  
  # ... rest of server (pairwise functions etc) unchanged
  build_pair_table <- function(df) {
    rep_id <- paste(df$Sample_ID, df$Row, df$Column, sep = "_")
    wide <- df %>% mutate(replicate = rep_id) %>%
      select(replicate, timepoint, Absorbance) %>%
      pivot_wider(names_from = timepoint, values_from = Absorbance)
    
    if (nrow(wide) == 0) return(data.frame())
    y_mat <- as.matrix(wide %>% select(-replicate))
    ord <- order(as.numeric(colnames(y_mat)))
    y_mat <- y_mat[, ord, drop = FALSE]
    grp <- factor(stringr::str_remove(wide$replicate, "_[A-H]_[0-9]+$"))
    keep_levels <- names(which(table(grp) >= 2))
    if (length(keep_levels) < 2) return(data.frame())
    
    g1 <- g2 <- character()
    for (i in 1:(length(keep_levels)-1)) {
      for (j in (i+1):length(keep_levels)) {
        g1 <- c(g1, keep_levels[i])
        g2 <- c(g2, keep_levels[j])
      }
    }
    data.frame(Sample1 = g1, Sample2 = g2, stringsAsFactors = FALSE)
  }
  
  stat_results <- eventReactive(input$compare_button, {
    req(filtered_data())
    df <- filtered_data()
    pair_tab <- build_pair_table(df)
    if (nrow(pair_tab) == 0) {
      showNotification("Not enough groups with >= 2 replicates for pairwise comparison", type = "error")
      return(data.frame())
    }
    
    out_rows <- vector("list", nrow(pair_tab))
    for (i in seq_len(nrow(pair_tab))) {
      s1 <- pair_tab$Sample1[i]; s2 <- pair_tab$Sample2[i]
      subdat <- df %>% filter(Sample_ID %in% c(s1, s2)) %>%
        mutate(replicate = paste(Sample_ID, Row, Column, sep = "_"))
      subdat$Sample_ID <- factor(subdat$Sample_ID)
      subdat$replicate <- factor(subdat$replicate)
      
      result_row <- list(Sample1 = s1, Sample2 = s2,
                         method = input$stat_method, p_value = NA_real_,
                         stat = NA_real_, extra = NA_character_)
      
      tryCatch({
        if (input$stat_method == "statmod") {
          wide <- subdat %>% select(replicate, timepoint, Absorbance) %>%
            pivot_wider(names_from = timepoint, values_from = Absorbance)
          y_mat <- as.matrix(wide %>% select(-replicate))
          ord <- order(as.numeric(colnames(y_mat)))
          y_mat <- y_mat[, ord, drop = FALSE]
          grp <- factor(stringr::str_remove(wide$replicate, "_[A-H]_[0-9]+$"))
          res <- compareGrowthCurves(group = grp, y = y_mat, nsim = input$nsim_input, adjust = "none", verbose = FALSE)
          hit <- res %>% filter(Group1 == s1 & Group2 == s2)
          if (nrow(hit) == 0) hit <- res %>% filter(Group1 == s2 & Group2 == s1)
          if (nrow(hit) > 0) { result_row$p_value <- hit$P.Value[1]; result_row$stat <- hit$Stat[1] }
        } else if (input$stat_method == "lme" || input$stat_method == "lme_spline") {
          subdat <- subdat %>% mutate(timepoint = as.numeric(timepoint))
          random_formula <- if (isTRUE(input$lme_random_intercept_slope)) as.formula("~ timepoint | replicate") else as.formula("~ 1 | replicate")
          
          if (input$stat_method == "lme_spline") {
            if (input$spline_or_poly == "ns") fixed_formula <- as.formula("Absorbance ~ ns(timepoint, df = 3) * Sample_ID")
            else { deg <- as.integer(input$poly_degree); fixed_formula <- as.formula(paste0("Absorbance ~ poly(timepoint, ", deg, ", raw = TRUE) * Sample_ID")) }
          } else fixed_formula <- as.formula("Absorbance ~ timepoint * Sample_ID")
          
          fit_full <- tryCatch(lme(fixed = fixed_formula, random = random_formula, data = subdat, na.action = na.omit,
                                   control = lmeControl(opt = "optim")), error = function(e) NULL)
          
          if (!is.null(fit_full)) {
            if (input$stat_method == "lme_spline") {
              if (input$spline_or_poly == "ns") reduced_formula <- as.formula("Absorbance ~ ns(timepoint, df = 3) + Sample_ID")
              else { deg <- as.integer(input$poly_degree); reduced_formula <- as.formula(paste0("Absorbance ~ poly(timepoint, ", deg, ", raw = TRUE) + Sample_ID")) }
            } else reduced_formula <- as.formula("Absorbance ~ timepoint + Sample_ID")
            
            fit_reduced <- tryCatch(lme(fixed = reduced_formula, random = random_formula, data = subdat, na.action = na.omit,
                                        control = lmeControl(opt = "optim")), error = function(e) NULL)
            
            if (!is.null(fit_reduced)) {
              anov <- tryCatch(anova(fit_reduced, fit_full), error = function(e) NULL)
              if (!is.null(anov) && "p-value" %in% colnames(anov)) {
                result_row$p_value <- anov$"p-value"[2]
                result_row$extra <- paste0("ModelDFs:", paste(anov$df, collapse = ";"))
              }
            }
          }
        } else if (input$stat_method == "rm_anova") {
          subdat <- subdat %>% mutate(timepoint = as.factor(timepoint))
          res_aov <- tryCatch(aov(Absorbance ~ timepoint * Sample_ID + Error(replicate/timepoint), data = subdat),
                              error = function(e) NULL)
          if (!is.null(res_aov)) {
            s <- summary(res_aov)
            pval <- NA_real_
            for (el in s) {
              if (is.list(el) && length(el) >= 1) {
                if (is.matrix(el[[1]]) || is.data.frame(el[[1]])) {
                  tb <- as.data.frame(el[[1]])
                  if ("timepoint:Sample_ID" %in% rownames(tb)) {
                    pval <- tb["timepoint:Sample_ID", "Pr(>F)"]
                  }
                }
              }
            }
            result_row$p_value <- pval
          }
        } else if (input$stat_method == "per_timepoint") {
          subdat <- subdat %>% mutate(timepoint = as.numeric(timepoint))
          tps <- sort(unique(subdat$timepoint))
          pvals <- numeric(length(tps))
          for (k in seq_along(tps)) {
            tp <- tps[k]
            x <- subdat %>% filter(timepoint == tp & Sample_ID == s1) %>% pull(Absorbance)
            y <- subdat %>% filter(timepoint == tp & Sample_ID == s2) %>% pull(Absorbance)
            if (length(x) < 2 || length(y) < 2) pvals[k] <- NA_real_ else {
              if (input$per_time_test == "ttest") {
                tt <- tryCatch(t.test(x, y), error = function(e) NULL)
                pvals[k] <- if (!is.null(tt)) tt$p.value else NA_real_
              } else {
                wt <- tryCatch(wilcox.test(x, y), error = function(e) NULL)
                pvals[k] <- if (!is.null(wt)) wt$p.value else NA_real_
              }
            }
          }
          adj <- p.adjust(pvals, method = input$per_time_adj)
          result_row$p_value <- min(adj, na.rm = TRUE)
          result_row$extra <- paste0("n_timepoints=", length(pvals), "; min_adj_p=", signif(result_row$p_value,4))
        }
      }, error = function(e) {
        result_row$extra <- paste0("ERROR: ", conditionMessage(e))
      })
      
      out_rows[[i]] <- result_row
    }
    
    res_df <- bind_rows(out_rows)
    res_df <- res_df %>% mutate(p_value_numeric = as.numeric(p_value))
    res_df <- res_df %>% mutate(adj_p = ifelse(!is.na(p_value_numeric), p.adjust(p_value_numeric, method = "BH"), NA_real_))
    res_df <- res_df %>% select(Sample1, Sample2, method, p_value = p_value_numeric, adj_p, stat, extra)
    res_df <- res_df %>% mutate(across(c(p_value, adj_p), ~ signif(.x, 4)))
    res_df
  })
  
  output$pairwise_pvalues <- renderDataTable({
    req(input$compare_button)
    req(stat_results())
    stat_results()
  }, options = list(pageLength = 10, scrollX = TRUE))
  
} # server

shinyApp(ui = ui, server = server)
