
library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RMySQL)
library(shinyjs)
library(DT)
library(RMariaDB)
library(stringr)
library(plotly)

p_values_df<-read.csv("data/pathwaypValues.csv")
path_celltypes<-unique(p_values_df$Celltype)
print(path_celltypes)

library(RMariaDB)
SSLloc<-"data/usertrust-intermediate-root.crt"
conditions<-c("JH01","JH02","JH03","JH04")
geneNames<-read.delim("data/geneNames_3.tsv",header = FALSE)[,1]
cellData2<-readRDS("data/neutrophil_cellData3.RDS")
rownames(cellData2)<-gsub("-", ".", rownames(cellData2))
print(dim(cellData2))
pathwayDF<-read.table("data/pathwayGenes.csv",sep = ",",header = TRUE)                          

database2<-"********"
query <- "SELECT DISTINCT gene_name FROM gene_expression;"
tryCatch(
  {
    con <-  dbConnect(RMariaDB::MariaDB(),user = '********',password ='********',host = '********',
                      port = 3306,dbname =database2,ssl.ca=SSLloc,mysql='true')
  },
  error = function(e) {
    # Attempt to reconnect or handle the error appropriately
    dbDisconnect(con); con <-  dbConnect(RMariaDB::MariaDB(),user = '********',password ='********',host = '********',
                                         port = 3306,dbname =database2,ssl.ca=SSLloc,mysql='true')
  }
)
geneNames_db <- dbGetQuery(con, query)
dbDisconnect(con)
geneNames_db<-geneNames_db[,1]


#############################
query <- paste0("SELECT * FROM cell_metadata;")
tryCatch(
  {
    con <-  dbConnect(RMariaDB::MariaDB(),user = '********',password ='********',host = '********',
                      port = 3306,dbname ='********',ssl.ca=SSLloc,mysql='true')
  },
  error = function(e) {
    # Attempt to reconnect or handle the error appropriately
    dbDisconnect(con); con <-  dbConnect(RMariaDB::MariaDB(),user = '********',password ='********',host = '********',
                                         port = 3306,dbname ='********',ssl.ca=SSLloc,mysql='true')
  }
)
cellData <- dbGetQuery(con, query)
dbDisconnect(con)
if(!grepl("\\.",cellData[1,1])){
  cellData[,1]<-gsub("-", ".", cellData[,1])
}
rownames(cellData)<-cellData[,1]
cellData<-cellData[,c(2,3,4,5,1)]
cellData<-cellData[cellData$orig_ident %in% conditions,]


############################
rownames(cellData)<-gsub("-",".",rownames(cellData))
cellTypes<-as.character(sort(unique(cellData$celltype)))
cellTypes2<-as.character(sort(unique(cellData2$celltype)))
Samples<-c("JH01","JH02","JH03","JH04") #unique(cellData$obj.orig.ident)

Samples<-conditions

#########################################################################################################################################
PathwayexpressionDot <- function(inputdf, pathway_name, subset,splits,xtext,ytext) {
  #inputdf # currently semicolon separated column 2
  # Extract unique genes from the input dataframe
  print(Sys.time())
  path_name<-pathway_name
  unique_genes <- unique(inputdf$gene_name)
  print(unique_genes)

  if (subset) {
    df <- expand.grid(
      cellType = unique(cellData2$celltype),
      gene = unique_genes,
      objOrigIdent = unique(cellData2$orig_ident)
    )
    
    df$avg.exp <- 0
    df$pct.exp <- 0
    df$pct.exp.SIZENORM <- 0
    
    # Calculate average expression and percentage for each celltype and gene
    for (i in 1:nrow(df)) {
      cells <- rownames(cellData2[cellData2$orig_ident == df$objOrigIdent[i] & 
                                    cellData2$celltype == df$cellType[i], ])
      
      gene_data <- inputdf[inputdf$gene_name == df$gene[i], ]
      
      expr_values <- gene_data[rownames(cellData2) %in% cells, 3]
      
      df$avg.exp[i] <- mean(as.numeric(expr_values), na.rm = TRUE)
      df$pct.exp[i] <- mean(as.numeric(expr_values) > 0, na.rm = TRUE) * 100
      df$pct.exp.SIZENORM[i] <- df$pct.exp[i] * (length(cells) / nrow(cellData2))
    }
  } else {
    if("All_celltypes" %in% splits){
      cellData_unique<-cellData
      inputdf2<-inputdf
    }else{
      cellData_unique<-cellData
      cellData_unique<-cellData_unique[cellData_unique$celltype %in% splits,]
      inputdf2<-inputdf[inputdf$cell_name %in% cellData_unique$cell_name,]
    }
    
    df <- expand.grid(
      cellType = unique(cellData_unique$celltype),
      gene = unique_genes,
      objOrigIdent = unique(cellData_unique$orig_ident)
    )
    
    df$avg.exp <- 0
    df$pct.exp <- 0
    df$pct.exp.SIZENORM <- 0
    
    merged_data <- inputdf2 %>%
      left_join(cellData_unique, by = "cell_name")
    
    df <- merged_data %>%
      group_by(orig_ident, celltype, gene_name) %>%
      summarize(
        avg.exp = mean(Expression, na.rm = TRUE),
        pct.exp = mean(Expression > 0, na.rm = TRUE) * 100,
        cells_in_group = n()) %>%
      ungroup()
    
    df <- df %>%
      left_join(
        merged_data %>%
          group_by(orig_ident, celltype) %>%
          summarize(total_cells = n()) %>%
          ungroup(),
        by = c("orig_ident", "celltype")
      ) %>%
      mutate(
        pct.exp.SIZENORM = pct.exp * (cells_in_group / dim(cellData)[1])
      )
  }
  print(Sys.time())
  df$celltype_orig<-paste(df$celltype,df$orig_ident,sep = "_")
  df<-df[df$orig_ident %in% splits,]


  ggplot(df, aes(x = gene_name, y = celltype_orig, size = pct.exp.SIZENORM, color = avg.exp)) +
    geom_point() +
    scale_size_continuous(range = c(1, 6)) +
    scale_color_gradient(low = "turquoise", high = "red", 
                         limits = c(min(-0.5, min(df$avg.exp, na.rm = TRUE)), 
                                    max(df$avg.exp, na.rm = TRUE))) +
    labs(title = paste0(path_name, " Expression Across Cell Types"), 
         x = "Genes", 
         y = "Cell Types") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = xtext, angle = 45, hjust = 1),  # Reactive X-axis text size
      axis.text.y = element_text(size = ytext),  # Reactive Y-axis text size
      plot.title = element_text(hjust = 0.5, size = 16),
      panel.grid = element_blank(),
      aspect.ratio = NULL)
}

######################################################################################################################################
expressionUmap<- function(inputdf,gene_name,pt_size,cells,sampleName,subset){


  tmpgene<-inputdf
  tmpgene<-tmpgene[,3]
  if(subset){
    tmpgene<-data.frame(cellData2[,1:4],tmpgene)
    tmpgene<-tmpgene[which(rownames(cellData2) %in% cells),]
  }else{
    tmpgene<-data.frame(cellData[,1:4],tmpgene)
    tmpgene<-tmpgene[which(rownames(cellData) %in% cells),]
  }
  colnames(tmpgene)[5]<-"expression"
  tmpgene <- tmpgene[order(tmpgene$expression,decreasing = FALSE), ]
  ggplot(tmpgene, aes(x = umap_1, y = umap_2, color = expression)) +
    geom_point(size = pt_size) +
    scale_color_gradient(low = "turquoise", high = "red") +
    ggtitle(paste0(gene_name," expression in\n",sampleName," cells")) +
    theme_minimal() +  
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      plot.margin = unit(rep(0, 4), "cm"),
      plot.title = element_text(hjust = 0.5),  # Center the ggtitle
      legend.key.width = unit(0.5, "cm"),  # Adjust the width of the legend key
      legend.key.height = unit(0.5, "cm")  # Adjust the height of the legend key
    ) +
    coord_fixed(ratio = 1)
}


#################################
### Function for gene DotPlots

expressionDot<-function(inputdf,gene_name,subset){
  # Subset the relevant data
  # Create a data frame for ggplot
  tmpgene<-inputdf
  tmpgene<-tmpgene[,3]
  if(subset){
    df <- data.frame(
      cellType = rep(unique(cellData2$celltype),length(unique(cellData2$orig_ident))),
      

      objOrigIdent = lapply(unique(cellData2$orig_ident), function(x) rep(x, each = length(unique(cellData2$celltype)))) %>% unlist(),
      
      avg.exp = 0,
      pct.exp = 0,
      pct.exp.SIZENORM=0
    )
    for(row in 1:dim(df)[1]){
      cells2<-rownames(cellData2[(cellData2$orig_ident==df$objOrigIdent[row])&(cellData2$celltype==df$cellType[row]),])
      
      df$avg.exp[row]<-mean(as.numeric(tmpgene[which(rownames(cellData2) %in% cells2)]))
      df$pct.exp[row]<-mean(as.numeric(tmpgene[which(rownames(cellData2) %in% cells2)]) > 0) * 100
      df$pct.exp.SIZENORM[row]<-(df$pct.exp[row]*(length(cells2)/dim(cellData2)[1]))
    }
  }else{
    df <- data.frame(
      cellType = rep(unique(cellData$celltype),length(unique(cellData$orig_ident))),
      

      objOrigIdent = lapply(unique(cellData$orig_ident), function(x) rep(x, each = length(unique(cellData$celltype)))) %>% unlist(),
      
      avg.exp = 0,
      pct.exp = 0,
      pct.exp.SIZENORM=0
    )
    for(row in 1:dim(df)[1]){
      cells2<-rownames(cellData[(cellData$orig_ident==df$objOrigIdent[row])&(cellData$celltype==df$cellType[row]),])
      
      df$avg.exp[row]<-mean(as.numeric(tmpgene[which(rownames(cellData) %in% cells2)]))
      df$pct.exp[row]<-mean(as.numeric(tmpgene[which(rownames(cellData) %in% cells2)]) > 0) * 100
      df$pct.exp.SIZENORM[row]<-(df$pct.exp[row]*(length(cells2)/dim(cellData)[1]))
    }
  }
  
  df$objOrigIdent<-gsub("JH01","JH01_MRSA_Infected_Kidney",df$objOrigIdent)
  df$objOrigIdent<-gsub("JH02","JH02_MRSA_Infected_Kidney",df$objOrigIdent)
  df$objOrigIdent<-gsub("JH03","JH03_MRSA_Infected_Kidney",df$objOrigIdent)
  df$objOrigIdent<-gsub("JH04","JH04_Control_noMRSA",df$objOrigIdent)

  # Plot
  ggplot(df, aes(x = objOrigIdent, y = cellType, size = pct.exp.SIZENORM, color = avg.exp)) +
    geom_point() +
    scale_size_continuous(range = c(1, 6)) +
    scale_color_gradient(low = "turquoise", high = "red", limits = c(min(-0.5,min(df$avg.exp, na.rm = TRUE)), max(df$avg.exp, na.rm = TRUE))) +
    labs(title = paste0(gene_name, " expression")) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      axis.title.x = element_blank(),  # Remove x-axis title
      axis.title.y = element_blank()   # Remove y-axis title
    ) +
    coord_fixed(ratio = 0.5) 
  
}

#################################
### Function for activity UMAPs

activityUmap<-function(inputdf,gene_name,pt_size,cells,sampleName,subset){
  tmpgene<-inputdf
  tmpgene<-tmpgene[,4]
  if(subset){
    tmpgene<-data.frame(cellData2[,1:4],tmpgene)
    tmpgene<-tmpgene[which(rownames(cellData2) %in% cells),]
  }else{
    tmpgene<-data.frame(cellData[,1:4],tmpgene)
    tmpgene<-tmpgene[which(rownames(cellData) %in% cells),]
  }
  colnames(tmpgene)[5]<-"activity"
  tmpgene <- tmpgene[order(tmpgene$activity,decreasing = FALSE), ]
  ggplot(tmpgene, aes(x = umap_1, y = umap_2, color = activity)) +
    geom_point(size = pt_size) +
    scale_color_gradient(low = "grey", high = "blue") +
    ggtitle(paste0(gene_name," activity in\n",sampleName," cells")) +
    theme_minimal() +  
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      plot.margin = unit(rep(0, 4), "cm"),
      plot.title = element_text(hjust = 0.5),  # Center the ggtitle
      legend.key.width = unit(0.5, "cm"),  # Adjust the width of the legend key
      legend.key.height = unit(0.5, "cm")  # Adjust the height of the legend key
      
    ) +
    coord_fixed(ratio = 1)
  
}

#################################
### Function for activity DotPlots

activityDot<-function(inputdf,gene_name,subset){
  # Subset the relevant data
  tmpgene<-inputdf
  tmpgene<-as.numeric(tmpgene[,4])
  if(subset){
    df <- data.frame(
      cellType = rep(unique(cellData2$celltype),length(unique(cellData2$orig_ident))),
      

      objOrigIdent = lapply(unique(cellData2$orig_ident), function(x) rep(x, each = length(unique(cellData2$celltype)))) %>% unlist(),
      
      avg.exp = 0,
      pct.exp = 0,
      pct.exp.SIZENORM=0
    )
    for(row in 1:dim(df)[1]){
      cells<-rownames(cellData2[(cellData2$orig_ident==df$objOrigIdent[row])&(cellData2$celltype==df$cellType[row]),])
      
      df$avg.exp[row]<-mean(as.numeric(tmpgene[which(rownames(cellData2) %in% cells)]))
      df$pct.exp[row]<-mean(as.numeric(tmpgene[which(rownames(cellData2) %in% cells)]) > 0) * 100
      df$pct.exp.SIZENORM[row]<-(df$pct.exp[row]*(length(cells)/dim(cellData2)[1]))
    }
  }else{
    df <- data.frame(
      cellType = rep(unique(cellData$celltype),length(unique(cellData$orig_ident))),
      
      objOrigIdent = lapply(unique(cellData$orig_ident), function(x) rep(x, each = length(unique(cellData$celltype)))) %>% unlist(),
      
      avg.exp = 0,
      pct.exp = 0,
      pct.exp.SIZENORM=0
    )
    for(row in 1:dim(df)[1]){
      cells<-rownames(cellData[(cellData$orig_ident==df$objOrigIdent[row])&(cellData$celltype==df$cellType[row]),])
      
      df$avg.exp[row]<-mean(as.numeric(tmpgene[which(rownames(cellData) %in% cells)]))
      df$pct.exp[row]<-mean(as.numeric(tmpgene[which(rownames(cellData) %in% cells)]) > 0) * 100
      df$pct.exp.SIZENORM[row]<-(df$pct.exp[row]*(length(cells)/dim(cellData)[1]))
    }
  }
  
  df$objOrigIdent<-gsub("JH01","JH01_MRSA_Infected_Kidney",df$objOrigIdent)
  df$objOrigIdent<-gsub("JH02","JH02_MRSA_Infected_Kidney",df$objOrigIdent)
  df$objOrigIdent<-gsub("JH03","JH03_MRSA_Infected_Kidney",df$objOrigIdent)
  df$objOrigIdent<-gsub("JH04","JH04_Control_noMRSA",df$objOrigIdent)
  # Plot
  ggplot(df, aes(x = objOrigIdent, y = cellType, size = pct.exp.SIZENORM, color = avg.exp)) +
    geom_point() +
    scale_size_continuous(range = c(1, 6)) +
    scale_color_gradient(low = "grey", high = "blue", limits = c(min(df$avg.exp, na.rm = TRUE), max(df$avg.exp, na.rm = TRUE))) +
    labs(title = paste0(gene_name," activity")) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      axis.title.x = element_blank(),  # Remove x-axis title
      axis.title.y = element_blank()   # Remove y-axis title
    ) +
    coord_fixed(ratio = 0.5) 
  
}


####### Find Markers Function for full object
FindMarkersWilcox <- function(
    data, 
    cells1 = NULL, 
    cells2 = NULL, 
    features = NULL, 
    logfc.threshold = 0.1, 
    min.pct = 0.01,
    p.adj = TRUE,
    ...) {
  
  wilcox_results <- lapply(features, function(feature) {
    group1 <- data[feature, cells1, drop = FALSE]  # Subset data for group 1 using row names
    group2 <- data[feature, cells2, drop = FALSE]  # Subset data for group 2 using row names
    
    # Perform Wilcoxon rank sum test for the specified feature
    wresult <- wilcox.test(group1, group2)
    
    # Calculate log2 fold change
    log2FC <- log2(mean(group1) + 1) - log2(mean(group2) + 1)
    
    p_val_adj <- wresult$p.value
    
    return(list(feature = feature, 
                p_value = wresult$p.value, 
                p_value_adj = p_val_adj,
                statistic = wresult$statistic,
                log2FC = log2FC))
  })
  
  # Convert results to data frame
  wilcox_results <- do.call(rbind, wilcox_results)
  wilcox_results <- as.data.frame(wilcox_results)
  
  if (p.adj) {
    p_val_adj <- p.adjust(wilcox_results$p_value, method = "bonferroni")
    wilcox_results$p_value_adj<-p_val_adj
  } else {}
  # Apply thresholds
  wilcox_results <- subset(wilcox_results, p_value_adj < min.pct)
  wilcox_results <- subset(wilcox_results, abs(as.numeric(log2FC)) > logfc.threshold)
  
  # Order results by p-value
  wilcox_results <- wilcox_results[order(wilcox_results$p_value_adj), ]
  
  return(wilcox_results[,-4])
}
#####################################################################################################################################


# Define UI for application
ui <- navbarPage(
  useShinyjs(),
  tags$hr(),  # Horizontal grey bar
  div(
    style = "text-align:center; background-color: #f0f0f0;",
    h4("MRSA Mouse Kidney Multiome")
  ),
  
  # Sidebar with a text input for gene
  tabPanel(
    "Plots",
    tags$hr(),  # Horizontal grey bar
    div(
      style = "text-align:center; background-color: #f0f0f0;",
      h4("Cell Clustering")
    ),
    fluidRow(
      splitLayout(
        cellWidths = c("50%","50%"),
        img(src = "umap.png", width = "auto", height = "750px", alt = "Samples", style = "object-fit: cover; margin-left: auto;"),
        img(src = "umap_number.png", width = "auto", height = "750px", alt = "clustNumber", style = "object-fit: cover; margin-left: auto;")
      )
      
    ),
    fluidRow(
      column(
        width = 3,
        selectizeInput(
          "gene_name",
          label = "Gene or pathway of interest:",
          choices = as.list(geneNames),
          multiple = FALSE,
          selected = "Pdk4"
        )
      ),
      column(
        width = 3,
        actionButton(
          inputId = "submit_loc",
          label = "Submit gene"
        )
      ),
      column(
        width = 3,
        sliderInput(
          inputId = "pt_size",
          label = "Point size",
          min = 0.0,
          max = 1,
          value = 0.30
        )
      )
    ),
    tags$hr(),  # Horizontal grey bar
    div(
      style = "text-align:center; background-color: #f0f0f0;",
      h4("Gene Expression")
    ),
    fluidRow(
      splitLayout(
        cellWidths = c("50%","50%"),
        plotOutput("plot1"),
        plotOutput("plot2")
      )
    ),
    tags$hr(),  # Horizontal grey bar
    div(
      style = "text-align:center; background-color: #f0f0f0;",
      h4("Gene Activity and Expression Vln")
    ),
    fluidRow(
      splitLayout(
        cellWidths = c("33%","33%","33%"),
        plotOutput("plot3"),
        plotOutput("plot4"),
        plotOutput("plotvln")
      )
    ),
    ######################################################################################
    ### conditional split, vln by cluster
    # Gene Expression by Condition
    tags$hr(),  # Horizontal grey bar
    div(
      style = "text-align:center; background-color: #f0f0f0;",
      h4("Split Vln plot")
    ),
    fluidRow(
      checkboxGroupInput("vln_split", "Select groups", choices = c(Samples, cellTypes), selected = c(Samples[1:4], cellTypes[8:11]),inline = TRUE)
      
    ),
    fluidRow(
      splitLayout(
        cellWidths = c("90%"),
        plotOutput("plotvln2")
      )
    ),# Gene Expression by Condition
    tags$hr(),  # Horizontal grey bar
    div(
      style = "text-align:center; background-color: #f0f0f0;",
      h4("Significant expression differences between selected groups")
    ),
    fluidRow(
      mainPanel(
        checkboxInput("filter_padj", "Filter Adjusted P-Values > 0.05", value = FALSE),
        DTOutput("statTable")
      )),
    ######################################################################################
    
    actionButton("toggle_splits", "Show Conditional Splits"), #############CHANGE FOR CAROLINE 4 SAMPLES
    
    # Checkbox inputs
    fluidRow(
      column(4, checkboxInput("check_JH01", "JH01_MRSA_infected_kidney", value = TRUE)),
      column(4, checkboxInput("check_JH02", "JH02_MRSA_infected_kidney", value = TRUE)),
      column(4, checkboxInput("check_JH03", "JH03_MRSA_infected_kidney", value = TRUE)),
      column(4, checkboxInput("check_JH04", "JH04_Control_noMRSA", value = TRUE))
    ),
    
    # Conditionally show/hide based on toggle button
    conditionalPanel(
      condition = "input.toggle_splits % 2 == 1",
      
      # Gene Expression by Condition
      tags$hr(),  # Horizontal grey bar
      div(
        style = "text-align:center; background-color: #f0f0f0;",
        h4("Gene Expression by Condition")
      ),
      fluidRow(
        splitLayout(
          cellWidths = c("24%", "24%", "24%", "24%"),
          plotOutput("plot5"),
          plotOutput("plot6"),
          plotOutput("plot7"),
          plotOutput("plot8")
        )
      ),
      
      # Horizontal grey bar
      tags$hr(),
      
      # Gene Activity by Condition
      div(
        style = "text-align:center; background-color: #f0f0f0;",
        h4("Gene Activity by Condition")
      ),
      fluidRow(
        splitLayout(
          cellWidths = c("24%", "24%", "24%", "24%"),
          plotOutput("plot9"),
          plotOutput("plot10"),
          plotOutput("plot11"),
          plotOutput("plot12")
        )
      )
    )
    
    
  ), ## end of tab panel
  tabPanel("Pathway Figures",
           mainPanel(
             fluidRow(
               column(3, selectInput(
                 inputId = "pathway_type",
                 label = "Select Pathway Type:",
                 choices = c("KEGG", "HALLMARK"),
                 selected = "HALLMARK"
               )),
               column(3, selectInput(
                 inputId = "celltypePath",
                 label = "Select Cell Type:",
                 choices = as.list(unique(p_values_df$Celltype)),
                 selected = "All"
               )),
               column(3, selectInput(
                 inputId = "comp1",
                 label = "Select Comparison 1:",
                 choices = as.list(conditions),
                 selected = conditions[1]
               )),
               column(3, selectInput(
                 inputId = "comp2",
                 label = "Select Comparison 2:",
                 choices = as.list(conditions),
                 selected = conditions[4]
               )),
               column(3, textInput(
                 inputId = "top_n",
                 label = "Number of Top Pathways (Type 'All' to show all):",
                 value = "35"
               )),
               column(3, sliderInput(
                 inputId = "plot_height",
                 label = "Adjust Plot Height:",
                 min = 300,
                 max = 1000,
                 value = 750,
                 step = 50
               ))
             ),
             plotlyOutput("pathway_plot", height = "auto")
           )
  ),
  
  tabPanel(
    "Pathway Views",
    tags$hr(),  # Horizontal grey bar
    div(
      style = "text-align:center; background-color: #f0f0f0;",
      h4("Pathway Expression Levels By Individual Gene for Each Celltype")
    ),
    
    sidebarLayout(
      sidebarPanel(
        width = 3,  # Sidebar width (adjustable)
        
        # Pathway selection
        selectizeInput(
          "path_name",
          label = "HALLMARK and mmu pathway:",
          choices = as.list(gsub("1$","",geneNames[grep("(^mmu|^HALLMARK)",geneNames)])),
          multiple = FALSE,
          selected = "mmu04210_Apoptosis"
        ),
        
        # Submit button
        actionButton(
          inputId = "submit_loc_path",
          label = "Submit Path"
        ),
        
        # Checkbox for selecting groups
        checkboxGroupInput(
          "path_split", "Select groups", 
          choices = c(Samples, cellTypes, "All_celltypes"), 
          selected = c(Samples[c(1,2)], cellTypes[c(15,16,17)]),
          inline = TRUE
        ),
        
        # Text size inputs
        numericInput("x_text_size", "X-axis (Gene) Text Size:", value = 12, min = 6, max = 24, step = 1),
        numericInput("y_text_size", "Y-axis (Cell Type) Text Size:", value = 14, min = 6, max = 24, step = 1),
        
        # Plot height slider
        sliderInput("plot_height", "Adjust Plot Height", min = 400, max = 1200, value = 600, step = 50),
        textInput("pdf_name", "PDF File Name:", value = "pathway_plot"),
        numericInput("pdf_width", "PDF Width:", value = 12, min = 5, max = 20, step = 1),
        numericInput("pdf_height", "PDF Height:", value = 6, min = 3, max = 15, step = 1),
        downloadButton("downloadPlot", "Download Plot as PDF")
      ),
      
      # Main panel for plot
      mainPanel(
        width = 9,  # Expands to fill remaining space
        div(
          style = "text-align:center; background-color: #f0f0f0;",
          h4("Pathway Expression"),
          p("May take over 10 seconds for larger pathways to appear")
        ),
        
        # Plot output with dynamic height
        plotOutput("plotPath")
      )
    )
  ),
  
  ####################### Find Markers tab ###########################################
  ####################################################################################
  tabPanel(
    "FindMarkers",
    titlePanel("Differential Expression Analysis"),
    sidebarLayout(
      sidebarPanel(
        column(6,
               checkboxGroupInput("group1", "Select Groups for Set 1:", choices = c(Samples, cellTypes), selected = c(Samples[1], cellTypes[1:2]))
        ),
        column(6,
               checkboxGroupInput("group2", "Select Groups for Set 2:", choices = c(Samples, cellTypes), selected = c(Samples[2], cellTypes[3:4]))
        ),

        actionButton("submit_DE_btn", "Initilize Data and Find Markers")
      ),
      mainPanel(
        dataTableOutput("top_de_table")  # Output for the data table
      )
    )
  ),
  ######################## end of Find Markers tab
  
  ################################### Neutrophil subset panel ##################################################################################
  tabPanel(
    "Neutrophil Plots",
    tags$hr(),  # Horizontal grey bar
    div(
      style = "text-align:center; background-color: #f0f0f0;",
      h4("Sub Clustering")
    ),
    fluidRow(
      splitLayout(
        cellWidths = c("50%","50%"),
        img(src = "neutrophil_cl_umap3.png", width = "auto", height = "750px", alt = "Samples", style = "object-fit: cover; margin-left: auto;"),
        img(src = "neutrophil_cl_umap3.png", width = "auto", height = "750px", alt = "clustNumber", style = "object-fit: cover; margin-left: auto;")
      )
      
    ),
    fluidRow(
      column(
        width = 3,
        selectizeInput(
          "gene_name2",
          label = "Gene of interest:",
          choices = as.list(geneNames),
          multiple = FALSE,
          selected = "Hmox1"
        )
      ),
      column(
        width = 3,
        actionButton(
          inputId = "submit_loc2",
          label = "Submit gene"
        )
      ),
      column(
        width = 3,
        sliderInput(
          inputId = "pt_size2",
          label = "Point size",
          min = 0.0,
          max = 1,
          value = 0.90
        )
      )
    ),
    tags$hr(),  # Horizontal grey bar
    div(
      style = "text-align:center; background-color: #f0f0f0;",
      h4("Gene Expression")
    ),
    fluidRow(
      splitLayout(
        cellWidths = c("50%","50%"),
        plotOutput("plot13"),
        plotOutput("plot14")
      )
    ),
    tags$hr(),  # Horizontal grey bar
    div(
      style = "text-align:center; background-color: #f0f0f0;",
      h4("Gene Activity and Expression Vln")
    ),
    fluidRow(
      splitLayout(
        cellWidths = c("33%","33%","33%"),
        plotOutput("plot15"),
        plotOutput("plot16"),
        plotOutput("plotvln3")
      )
    ),
    ######################################################################################
    ### conditional split, vln by cluster
    # Gene Expression by Condition
    tags$hr(),  # Horizontal grey bar
    div(
      style = "text-align:center; background-color: #f0f0f0;",
      h4("Split Vln plot")
    ),
    fluidRow(
      checkboxGroupInput("vln_split2", "Select groups", choices = c(Samples, cellTypes2), selected = c(Samples[1:4], cellTypes2[2:4]),inline = TRUE)
      
    ),
    fluidRow(
      splitLayout(
        cellWidths = c("90%"),
        plotOutput("plotvln4")
      )
    ),
    ######################################################################################
    
    actionButton("toggle_splits2", "Show Conditional Splits"), #############CHANGE FOR CAROLINE 4 SAMPLES
    
    # Checkbox inputs
    fluidRow(
      column(4, checkboxInput("check_JH012", "JH01_MRSA_infected_kidney", value = TRUE)),
      column(4, checkboxInput("check_JH022", "JH02_MRSA_infected_kidney", value = TRUE)),
      column(4, checkboxInput("check_JH032", "JH03_MRSA_infected_kidney", value = TRUE)),
      column(4, checkboxInput("check_JH042", "JH04_Control_noMRSA", value = TRUE))
    ),
    
    # Conditionally show/hide based on toggle button
    conditionalPanel(
      condition = "input.toggle_splits2 % 2 == 1",
      
      # Gene Expression by Condition
      tags$hr(),  # Horizontal grey bar
      div(
        style = "text-align:center; background-color: #f0f0f0;",
        h4("Gene Expression by Condition")
      ),
      fluidRow(
        splitLayout(
          cellWidths = c("24%", "24%", "24%", "24%"),
          plotOutput("plot17"),
          plotOutput("plot18"),
          plotOutput("plot19"),
          plotOutput("plot20")
        )
      ),
      
      # Horizontal grey bar
      tags$hr(),
      
      # Gene Activity by Condition
      div(
        style = "text-align:center; background-color: #f0f0f0;",
        h4("Gene Activity by Condition")
      ),
      fluidRow(
        splitLayout(
          cellWidths = c("24%", "24%", "24%", "24%"),
          plotOutput("plot21"),
          plotOutput("plot22"),
          plotOutput("plot23"),
          plotOutput("plot24")
        )
      )
    )
    
    
  ),
  ##############################################################################################################################################
  
  ########## FIND MARKERS FOR SUBSET
  tabPanel(
    "Neutrophil FindMarkers",
    titlePanel("Differential Expression Analysis for Subset"),
    sidebarLayout(
      sidebarPanel(
        column(6,
               checkboxGroupInput("group1_2", "Select Groups for Set 1:", choices = c(Samples, cellTypes2), selected = c(Samples[1], cellTypes2[1:2]))
        ),
        column(6,
               checkboxGroupInput("group2_2", "Select Groups for Set 2:", choices = c(Samples, cellTypes2), selected = c(Samples[2], cellTypes2[3:4]))
        ),
        #numericInput("num_top_de", "Number of Top DE Genes to Display:", value = 25, min = 1, max = 100),
        actionButton("submit_DE_btn2", "Initilize Data and Find Markers")
      ),
      mainPanel(
        dataTableOutput("top_de_table2")  # Output for the data table
      )
    )
  ),
  ################################################################################################################################################
  tabPanel(
    "Info Tab",
    fluidRow(
      column(
        width = 6,
        h4("Mouse Multiome dataset JH"),
        p("JH_01 MRSA infected kidney"),
        p("JH_02 MRSA infected kidney"),
        p("JH_03 MRSA infected kidney"),
        p("JH_04 Control kidney without MRSA infection")
      )
    )
  ),
  # Footer.
  hr(),
  fluidRow(
    column(width = 1, align = "center", img(src="https://www.iu.edu/images/brand/brand-expression/iu-trident-promo.jpg", width='100%')),
    column(width = 11,
           p( "Jered Myslinski", a("jmyslins@iu.edu",href="mailto:jmyslins@iu.edu"),
              br(),
              "Takashi Hato", a("thato@iu.edu", href="mailto:thato@iu.edu")
           )
    )
  )
)


server <- function(input, output) {
  filtered_data <- reactive({
    df <- p_values_df
    
    # Filter pathway type
    if (input$pathway_type == "HALLMARK") {
      df <- df %>% dplyr::filter(grepl("^HALLMARK", Pathway))
    } else {
      df <- df %>% dplyr::filter(grepl("^mmu", Pathway))
    }
    
    # Filter by cell type
    df <- df %>% dplyr::filter(Celltype == input$celltypePath)
    
    # Filter by selected comparisons
    df <- df[df$Comp1==as.character(input$comp1)|df$Comp1==as.character(input$comp2),]
    df <- df[df$Comp2==as.character(input$comp1)|df$Comp2==as.character(input$comp2),]
    #df <- df %>% dplyr::filter(Comp1 == input$comp1, Comp2 == input$comp2)
    
    
    # df <- df[df$Comp1==as.character(input$comp1),]
    # Limit to top N pathways if not "All"
    if (tolower(input$top_n) != "all") {
      top_n <- as.numeric(input$top_n)
      if (!is.na(top_n) && nrow(df) > top_n) {
        df <- df %>% arrange(desc(Padj)) %>% head(top_n)
      }
    }
    
    print(dim(df))
    df
  })
  
  output$pathway_plot <- renderPlotly({
    df <- filtered_data()
    
    # Convert Padj to -log10 scale
    df$Padj <- as.numeric(df$Padj)
    df <- df[!is.na(df$Padj), ]
    df$Log10PadjValue <- -log10(df$Padj)
    
    # Plot with grouped comparisons
    p <- ggplot(df, aes(x = reorder(Pathway, Log10PadjValue), y = Log10PadjValue, color = Celltype)) +
      geom_point(size = 3) +
      coord_flip() +
      labs(
        title = paste0("Significant Pathways: ", input$comp1, " vs ", input$comp2, " - ", input$pathway_type),
        x = "Pathways",
        y = "-log10(P-Value)"
      ) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 10))
    
    ggplotly(p) %>% layout(height = input$plot_height)
  })
  
  
  #get gene expression and gene activity data for select gene for plots (not findmarkers table)
  
  
  
  database2<-"jmysapps_caroline_mouse_multiome"
  genesql2 <- eventReactive(input[["submit_loc_path"]],
                            {
                              print(paste("genesql2pushed"))
                              matches <- geneNames_db[sapply(geneNames_db, function(x) grepl(paste0("^", x), input$path_name))]
                              print(input$path_name)
                              
                              inputPath<-as.character(input$path_name)
                              genes_string <- pathwayDF[pathwayDF$Pathway == gsub("\\.","-",inputPath), "Genes"]
                              genes_list <- unlist(strsplit(genes_string, ";"))  # Split genes by semicolon
                              genes_in_clause <- paste0("'", paste(genes_list, collapse = "', '"), "'")
                              
                              query <- paste0("SELECT * FROM gene_expression WHERE gene_name IN (", genes_in_clause, ");")
                              
                              tryCatch(
                                {
                                  con <-  dbConnect(RMariaDB::MariaDB(),user = '********',password ='********',host = '********',
                                                    port = 3306,dbname =database2,ssl.ca=SSLloc,mysql='true')
                                },
                                error = function(e) {
                                  # Attempt to reconnect or handle the error appropriately
                                  dbDisconnect(con); con <-  dbConnect(RMariaDB::MariaDB(),user = '********',password ='********',host = '********',
                                                                       port = 3306,dbname =database2,ssl.ca=SSLloc,mysql='true')
                                }
                              )
                              resultdf <- dbGetQuery(con, query)
                              dbDisconnect(con)
                              resultdf<-resultdf[,c(1,2)]
                              print(dim(resultdf))
                              cellnames<-rownames(cellData)
                            
                              
                              # Loop through each gene and process separately
                              split_values_list <- lapply(unique(resultdf$gene_name), function(gene) {
                                resultdf %>%
                                  dplyr::filter(gene_name == gene) %>%
                                  separate_rows(expression_level, sep = ";") %>%
                                  mutate(cell_name = rownames(cellData)) %>%
                                  mutate(scaled_expression = as.vector(scale(as.numeric(expression_level)))) %>%
                                  dplyr::select(gene_name, cell_name, scaled_expression)  # Explicit dplyr::select()
                              })
                              
                              # Combine all processed genes into one dataframe
                              split_values_df <- bind_rows(split_values_list)
                              split_values_df<-rbind(c(inputPath,inputPath,0),split_values_df)
                              colnames(split_values_df)[3]<-"Expression"
                              split_values_df
                              
                            }
  )
  
  
  #get gene expression and gene activity data for select gene for plots (not findmarkers table)
  
  #database2<-"jmysapps_basile_rat_snRNA"
  genesql <- eventReactive(input[["submit_loc"]],
                           {
                             if(length(grep("^(HALLMARK|mmu)",input$gene_name))==0){
                               matches <- input$gene_name
                             }else{
                               matches <- geneNames_db[sapply(geneNames_db, function(x) grepl(paste0("^", x), input$gene_name))]
                             }
                             query <- paste0("SELECT * FROM gene_expression WHERE gene_name = '", matches, "';")
                             tryCatch(
                               {
                                 con <-  dbConnect(RMariaDB::MariaDB(),user = '********',password ='********',host = '********',
                                                   port = 3306,dbname =database2,ssl.ca=SSLloc,mysql='true')
                               },
                               error = function(e) {
                                 # Attempt to reconnect or handle the error appropriately
                                 dbDisconnect(con); con <-  dbConnect(RMariaDB::MariaDB(),user = '********',password ='********',host = '********',
                                                                      port = 3306,dbname =database2,ssl.ca=SSLloc,mysql='true')
                               }
                             )
                             resultdf <- dbGetQuery(con, query)
                             dbDisconnect(con)
                             cellnames<-rownames(cellData)
                             split_values_X<- lapply(strsplit(as.character(resultdf), split = ";", fixed = TRUE), as.numeric)[[2]]
                             split_values_A<- lapply(strsplit(as.character(resultdf), split = ";", fixed = TRUE), as.numeric)[[3]]
                             
                             split_values_df <- t(rbind(rep(resultdf[,1],length(split_values_X)),cellnames,as.numeric(split_values_X),as.numeric(split_values_A)))
                             rownames(split_values_df) <- split_values_df[, 2]
                             colnames(split_values_df)[1:2]<-c("gene_name","cell_name")
                             colnames(split_values_df)[3]<-gsub(" ","",paste(input$gene_name,c("_Expression")))
                             colnames(split_values_df)[4]<-gsub(" ","",paste(input$gene_name,c("_Activity")))
                             
                             split_values_df<-as.data.frame(split_values_df)
                             split_values_df[,3]<-as.numeric(split_values_df[,3])
                             split_values_df[,3]<-scale(split_values_df[,3])
                             split_values_df[,4]<-as.numeric(split_values_df[,4])
                             split_values_df[,4]<-scale(split_values_df[,4])
                             split_values_df
                           }
  )
  
  plot_reactive <- reactive({
    inputdf <- data.frame(genesql2())
    print("pathwaysdf")
    print(head(inputdf))
    
    path_name <- inputdf[1, 1]
    print(path_name)
    
    inputdf <- inputdf[-1,]
    inputdf <- as.data.frame(inputdf)
    inputdf[, 3] <- as.numeric(inputdf[, 3])
    
    print(head(inputdf))
    
    choiceTypes <- input$path_split
    choiceTypes <- as.character(choiceTypes)
    print(choiceTypes)
    
    xtext <- input$x_text_size
    ytext <- input$y_text_size
    
    p <- PathwayexpressionDot(inputdf = inputdf, pathway_name = path_name, subset = FALSE, splits = choiceTypes, xtext = xtext, ytext = ytext)
    
    return(p)  # Store the plot
  })
  
  output$plotPath <- renderPlot({
    plot_reactive()
  }, height = reactive({input$plot_height}))
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0(input$pdf_name, ".pdf")  # Use user-defined name
    },
    
    content = function(file) {
      pdf(file, width = input$pdf_width, height = input$pdf_height, useDingbats = FALSE)  # Save as PDF
      print(plot_reactive())  # Print the stored plot
      dev.off()
    }
  )
  
  stat_results <- reactive({
    cellsChosen <- as.character(input$vln_split)
    gene_expression_df <- data.frame(genesql())
    gene_expression_df$cellTypes <- cellData$celltype
    gene_expression_df$obj.orig.ident <- cellData$orig_ident
    gene_expression_df$Cell_Groups <- paste0(gene_expression_df$cellTypes, " ", gene_expression_df$obj.orig.ident)
    
    # Filter based on input
    gene_expression_df <- gene_expression_df %>%
      dplyr::filter(obj.orig.ident %in% cellsChosen, cellTypes %in% cellsChosen)
    
    colnames(gene_expression_df)[3] <- "expression"
    
    # Compute mean expression per group
    mean_expression <- gene_expression_df %>%
      group_by(Cell_Groups) %>%
      summarize(avg_expression = mean(expression, na.rm = TRUE), .groups = 'drop')
    
    # Generate pairwise comparisons
    comparisons <- combn(unique(gene_expression_df$Cell_Groups), 2, simplify = FALSE)
    
    # Perform statistical tests
    results_list <- lapply(comparisons, function(pair) {
      group1 <- pair[1]
      group2 <- pair[2]
      
      expr1 <- gene_expression_df$expression[gene_expression_df$Cell_Groups == group1]
      expr2 <- gene_expression_df$expression[gene_expression_df$Cell_Groups == group2]
      
      if (length(expr1) > 1 && length(expr2) > 1) {
        p_value <- wilcox.test(expr1, expr2)$p.value  # Wilcoxon rank-sum test
      } else {
        p_value <- NA  # Avoid errors when sample size is too small
      }
      
      data.frame(Group1 = group1, Group2 = group2, p_value = p_value)
    })
    
    results_df <- do.call(rbind, results_list)
    results_df <- results_df[!is.na(results_df$p_value), ]
    # Adjust p-values using Benjamini-Hochberg correction
    results_df$p_adj <- p.adjust(results_df$p_value, method = "BH")
    print(head(results_df))
    print(head(mean_expression))
    # Merge with mean expression values
    
    # Ensure column names are clean
    colnames(results_df) <- trimws(colnames(results_df))
    
    # Print column names to check
    print(colnames(results_df))
    
    # Rename safely
    
    results_df <- results_df %>%
      left_join(mean_expression, by = c("Group1" = "Cell_Groups")) %>%
      rename(avg_exp_Group1 = avg_expression) %>%
      left_join(mean_expression, by = c("Group2" = "Cell_Groups")) %>%
      rename(avg_exp_Group2 = avg_expression)
    
    return(results_df)
  })
  
  # Render DataTable
  output$statTable <- renderDT({
    data <- stat_results()
    
    # Apply filtering if checkbox is checked
    if (input$filter_padj) {
      data <- data %>% dplyr::filter(p_adj <= 0.05)
    }
    
    data <- data %>%
      mutate(
        avg_exp_Group1 = round(avg_exp_Group1, 3),
        avg_exp_Group2 = round(avg_exp_Group2, 3)
      )
    
    datatable(data, options = list(pageLength = 10), rownames = FALSE) %>%
      formatSignif(columns = c("p_value", "p_adj"), digits = 3)  # Display in scientific notation while keeping numeric sorting
  })
  
  #########################
  genesql3 <- eventReactive(input[["submit_loc2"]],
                            {query <- paste0("SELECT * FROM gene_expression WHERE gene_name = '", input$gene_name2, "';")
                            tryCatch(
                              {
                                con <-  dbConnect(RMariaDB::MariaDB(),user = '********',password ='********',host = '********',
                                                  port = 3306,dbname ='********',ssl.ca=SSLloc,mysql='true')
                              },
                              error = function(e) {
                                # Attempt to reconnect or handle the error appropriately
                                dbDisconnect(con); con <-  dbConnect(RMariaDB::MariaDB(),user = '********',password ='********',host = '********',
                                                                     port = 3306,dbname ='********',ssl.ca=SSLloc,mysql='true')
                              }
                            )
                            resultdf2 <- dbGetQuery(con, query)
                            dbDisconnect(con)
                            colnames(resultdf2)[3:4]<-gsub(" ","",paste(input$gene_name,c("_Expression","_Activity")))
                            resultdf2[,3]<-scale(resultdf2[,3])
                            resultdf2[,4]<-scale(resultdf2[,4])
                            resultdf2
                            }
  )
  ########################################################
  ###Find Markers TESTS
  ########################################################
  rv <- reactiveValues(data = NULL) #stores table once it is downloaded from SQL server
  wresults_reactive <- reactiveVal(NULL)  # Initialize a reactive variable
  
  observeEvent(input$submit_DE_btn, {
    if (is.null(rv$data)) {  # Check if data is not already loaded
      
      query <- paste0("SELECT gene, geneexpression
    FROM caroline_expr_activity
    WHERE geneexpression NOT LIKE '0;0;0;0;0;%';") #gets all geneexpressions that are not all zeroes
      
      tryCatch(
        {
          con <-  dbConnect(RMariaDB::MariaDB(),user = '********',password ='********',host = '********',
                            port = 3306,dbname ='********',ssl.ca=SSLloc,mysql='true')
        },
        error = function(e) {
          # Attempt to reconnect or handle the error appropriately
          dbDisconnect(con); con <-  dbConnect(RMariaDB::MariaDB(),user = '********',password ='********',host = '********',
                                               port = 3306,dbname ='********',ssl.ca=SSLloc,mysql='true')
        }
      )
      withProgress(message = "Extracting Data from SQL database...", value = 0, {
        result <- dbGetQuery(con, query)
        dbDisconnect(con)
        split_values <- lapply(strsplit(as.character(result$geneexpression), split = ";", fixed = TRUE), as.numeric)
        split_values_df <- do.call(rbind, split_values)
        rownames(split_values_df) <- result[, 1]
        # Assign data to reactiveValues
        rv$data <- split_values_df #table of all genes and expressions (columns are all cells in order)
      })
    }
    # Filter data based on selected groups
    ################
    
    
    ###############
    withProgress(message = "Wilcox Tests running, finding Markers..", value = 0, {
      
      tryCatch({
        subset1 <- unique(which((cellData$celltype %in% as.character(input$group1)) & 
                                  (cellData$orig_ident %in% as.character(input$group1))))
        subset2 <- unique(which((cellData$celltype %in% as.character(input$group2)) & 
                                  (cellData$orig_ident %in% as.character(input$group2))))
        
        if (is.null(subset1) || is.null(subset2)) {
          showNotification("Error: One or both groups are empty.", type = "warning")
          return(NULL)
        } else if ((all((as.character(input$group1) %in% cellData$celltype)==FALSE) | 
                    all((as.character(input$group1) %in% cellData$orig_ident)==FALSE)) ||
                   (all((as.character(input$group2) %in% cellData$celltype)==FALSE) | 
                    all((as.character(input$group2) %in% cellData$orig_ident)==FALSE))) {
          showNotification("Error: One or both groups do not contain valid elements.", type = "warning")
          return(NULL)
        }
        
        # Run Wilcoxon test
        withProgress(message = "Wilcox Tests running, finding Markers..", value = 0, {
          features <- rownames(rv$data)
          wresults <- FindMarkersWilcox(
            data = rv$data,
            cells1 = subset1,
            cells2 = subset2,
            features = features,
            logfc.threshold = 0.1,
            min.pct = 0.01,
            only.pos = FALSE
          )
        })
        wresults[] <- lapply(wresults, function(x) if (is.list(x)) unlist(x) else x)
        
        # Ensure log2FC is numeric (if applicable)
        if ("log2FC" %in% colnames(wresults)) {
          wresults$log2FC <- as.numeric(wresults$log2FC)
        }
        
        # Store cleaned data
        wresults_reactive(wresults)
      }, error = function(e) {
        showNotification(paste("An error occurred:", e$message), type = "error")
      })
      
      incProgress(1, detail = "Analysis completed.")
      
      output$top_de_table <- renderDataTable({
        req(wresults_reactive())
        
        datatable(
          wresults_reactive(), 
          options = list(pageLength = 25),
          rownames = FALSE
        ) %>% formatRound(columns = "log2FC", digits = 3)  # Ensures numeric format
      })
    })
  })
  
  
  ########################################################
  ##############################################################################################################################
  
  ##Find markers subset ##
  
  observeEvent(input$submit_DE_btn2, {
    if (is.null(rv$data)) {  # Check if data is not already loaded
      
      query <- paste0("SELECT gene, geneexpression
    FROM caroline_expr_activity
    WHERE geneexpression NOT LIKE '0;0;0;0;0;%';") #gets all geneexpressions that are not all zeroes
      
      tryCatch(
        {
          con <-  dbConnect(RMariaDB::MariaDB(),user = '********',password ='********',host = '********',
                            port = 3306,dbname ='********',ssl.ca=SSLloc,mysql='true')
        },
        error = function(e) {
          # Attempt to reconnect or handle the error appropriately
          dbDisconnect(con); con <-  dbConnect(RMariaDB::MariaDB(),user = '********',password ='********',host = '********',
                                               port = 3306,dbname ='********',ssl.ca=SSLloc,mysql='true')
        }
      )
      withProgress(message = "Extracting Data from SQL database...", value = 0, {
        result <- dbGetQuery(con, query)
        dbDisconnect(con)
        split_values <- lapply(strsplit(as.character(result$geneexpression), split = ";", fixed = TRUE), as.numeric)
        split_values_df <- do.call(rbind, split_values)
        rownames(split_values_df) <- result[, 1]
        # Assign data to reactiveValues
        rv$data <- split_values_df #table of all genes and expressions (columns are all cells in order)
      })
    }
    # Filter data based on selected groups
    ################
    
    
    ###############
    withProgress(message = "Wilcox Tests running, finding Markers..", value = 0, {
      
      tryCatch({
        subset1 <- rownames(cellData2)[unique(which((cellData2$celltype %in% as.character(input$group1_2)) & 
                                                      (cellData2$orig_ident %in% as.character(input$group1_2))))]
        subset2 <- rownames(cellData2)[unique(which((cellData2$celltype %in% as.character(input$group2_2)) & 
                                                      (cellData2$orig_ident %in% as.character(input$group2_2))))]
        subset1 <- which(rownames(cellData) %in% subset1)
        subset2 <- which(rownames(cellData) %in% subset2)
        if (is.null(subset1) || is.null(subset2)) {
          showNotification("Error: One or both groups are empty.", type = "warning")
          return(NULL)
        } else if ((all((as.character(input$group1_2) %in% cellData2$celltype)==FALSE) | all((as.character(input$group1_2) %in% cellData2$orig_ident)==FALSE)) ||
                   (all((as.character(input$group2_2) %in% cellData2$celltype)==FALSE) | all((as.character(input$group2_2) %in% cellData2$orig_ident)==FALSE))) {
          showNotification("Error: One or both groups do not contain valid elements.", type = "warning")
          return(NULL)
        }
        
        # Define a function to find markers using Wilcoxon rank sum test
        withProgress(message = "Wilcox Tests running, finding Markers..", value = 0, {
          features <- rownames(rv$data)
          wresults <- FindMarkersWilcox(data = rv$data,
                                        cells1 = subset1,
                                        cells2 = subset2,
                                        features = features,
                                        logfc.threshold = 0.1,
                                        min.pct = 0.01,
                                        only.pos = FALSE)
        })
      }, error = function(e) {
        showNotification("An error occurred: " + e$message, type = "error")
      })
      incProgress(1, detail = "Analysis completed.")
      
      output$top_de_table2 <- renderDataTable({
        
        # Convert the data frame to a DataTable object
        datatable(wresults,
                  options = list(pageLength = 25))  # Sets the default number of rows displayed
      })
    })
  })
  
  
  ###################################################################################################################################################
  
  ##############################################For 'Plots'##################################################  
  # Generate a reactive list of plots that changes when the gene changes.
  
  output$plot1 <- renderPlot({
    # print(input$vln_celltypes)
    
    inputdf<-data.frame(genesql())
    gene_name<-inputdf[1,1]
    pt_size<-input$pt_size
    cells<-inputdf[,2]
    p<-expressionUmap(inputdf=inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = "all",subset=FALSE)
    print(p)
  })
  
  output$plot2 <- renderPlot({
    inputdf<-data.frame(genesql())
    gene_name<-inputdf[1,1]
    
    p<-expressionDot(inputdf = inputdf,gene_name = gene_name,subset=FALSE)
    print(p)
    
  })
  
  
  output$plot3 <- renderPlot({
    inputdf<-data.frame(genesql())
    gene_name<-inputdf[1,1]
    pt_size<-input$pt_size
    cells<-inputdf[,2]
    p<-activityUmap(inputdf = inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = "all",subset=FALSE)
    print(p)
  })
  
  output$plot4 <- renderPlot({
    inputdf<-data.frame(genesql())
    gene_name<-inputdf[1,1]
    p<-activityDot(inputdf=inputdf,gene_name = gene_name,subset=FALSE)
    print(p)
    
  })
  
  output$plotvln <- renderPlot({
    
    
    gene_expression_df<-data.frame(genesql())
    gene<-gene_expression_df[1,1]
    gene_expression_df$cellTypes <- cellData$celltype
    gene_expression_df$obj.orig.ident <- cellData$orig_ident
    gene_expression_df$Cell_Groups <- paste0(gene_expression_df$cellTypes," ",gene_expression_df$obj.orig.ident)
    set.seed(123)  # Set seed for reproducibility
    noise <- rnorm(n = nrow(gene_expression_df)) / 100000
    gene_expression_df[,3] <- gene_expression_df[,3] + noise
    colnames(gene_expression_df)[3]<-"expression"
    
    color_palette <- c("red", "blue", "green", "purple","orange","grey")
    
    p<-ggplot(gene_expression_df, aes(x = obj.orig.ident, y = expression, fill = obj.orig.ident)) +
      geom_violin(trim = TRUE, scale = "width") +
      geom_jitter(size=input$pt_size) +
      scale_fill_manual(values = setNames(color_palette, unique(gene_expression_df$obj.orig.ident))) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      ggtitle(paste0(gene," expression")) +
      labs(x = "Condition Groups")+
      labs(y = paste0(gene))+
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        plot.margin = unit(rep(0, 4), "cm"),
        plot.title = element_text(hjust = 0.5),  # Center the ggtitle
        legend.key.width = unit(0.5, "cm"),  # Adjust the width of the legend key
        legend.key.height = unit(0.5, "cm")  # Adjust the height of the legend key
      )
    
    print(p)
  })
  
  output$plotvln2 <- renderPlot({
    cellsChosen<-as.character(input$vln_split)
    gene_expression_df<-data.frame(genesql())
    gene<-gene_expression_df[1,1]
    gene_expression_df$cellTypes <- cellData$celltype
    gene_expression_df$obj.orig.ident <- cellData$orig_ident
    gene_expression_df$Cell_Groups <- paste0(gene_expression_df$cellTypes," ",gene_expression_df$obj.orig.ident)
    gene_expression_df<-gene_expression_df[gene_expression_df$obj.orig.ident %in% cellsChosen,]
    gene_expression_df<-gene_expression_df[gene_expression_df$cellTypes %in% cellsChosen,]
    set.seed(123)  # Set seed for reproducibility
    noise <- rnorm(n = nrow(gene_expression_df)) / 100000
    gene_expression_df[,3] <- gene_expression_df[,3] + noise
    colnames(gene_expression_df)[3]<-"expression"
    
    color_palette <- c("red", "blue", "green", "purple","orange","grey")
    gene_expression_df$Cell_Groups<-factor(gene_expression_df$Cell_Groups,levels = sort(unique(gene_expression_df$Cell_Groups)))
    p<-ggplot(gene_expression_df, aes(x = Cell_Groups, y = expression, fill = obj.orig.ident)) +
      geom_violin(trim = TRUE, scale = "width") +
      geom_jitter(size=input$pt_size) +
      scale_fill_manual(values = setNames(color_palette, unique(gene_expression_df$obj.orig.ident))) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      ggtitle(paste0(gene," expression")) +
      labs(x = "Condition Groups")+
      labs(y = paste0(gene))+
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        plot.margin = unit(rep(0, 4), "cm"),
        plot.title = element_text(hjust = 0.5),  # Center the ggtitle
        legend.key.width = unit(0.5, "cm"),  # Adjust the width of the legend key
        legend.key.height = unit(0.5, "cm")  # Adjust the height of the legend key
      )
    
    print(p)
  })
  #expression subsets
  output$plot5 <- renderPlot({
    if(input$check_JH01){
      inputdf<-data.frame(genesql())
      gene_name<-inputdf[1,1]
      pt_size<-input$pt_size
      cells<-rownames(cellData)[cellData$orig_ident==conditions[1]]
      p<-expressionUmap(inputdf = inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = conditions[1],subset=FALSE)
      print(p)}
  })
  output$plot6 <- renderPlot({
    if(input$check_JH02){
      inputdf<-data.frame(genesql())
      gene_name<-inputdf[1,1]
      pt_size<-input$pt_size
      cells<-rownames(cellData)[cellData$orig_ident==conditions[2]]
      p<-expressionUmap(inputdf = inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = conditions[2],subset=FALSE)
      print(p)}
  })
  output$plot7 <- renderPlot({
    if(input$check_JH03){
      inputdf<-data.frame(genesql())
      gene_name<-inputdf[1,1]
      pt_size<-input$pt_size
      cells<-rownames(cellData)[cellData$orig_ident==conditions[3]]
      p<-expressionUmap(inputdf = inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = conditions[3],subset=FALSE)
      print(p)}
  })
  output$plot8 <- renderPlot({
    if(input$check_JH04){
      inputdf<-data.frame(genesql())
      gene_name<-inputdf[1,1]
      pt_size<-input$pt_size
      cells<-rownames(cellData)[cellData$orig_ident==conditions[4]]
      p<-expressionUmap(inputdf = inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = conditions[4],subset=FALSE)
      print(p)}
  })
  
  
  #activity subsets
  output$plot9 <- renderPlot({
    if(input$check_JH01){
      inputdf<-data.frame(genesql())
      gene_name<-inputdf[1,1]
      pt_size<-input$pt_size
      cells<-rownames(cellData)[cellData$orig_ident==conditions[1]]
      p<-activityUmap(inputdf = inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = conditions[1],subset=FALSE)
      print(p)}
  })
  output$plot10 <- renderPlot({
    if(input$check_JH02){
      inputdf<-data.frame(genesql())
      gene_name<-inputdf[1,1]
      pt_size<-input$pt_size
      cells<-rownames(cellData)[cellData$orig_ident==conditions[2]]
      p<-activityUmap(inputdf = inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = conditions[2],subset=FALSE)
      print(p)}
  })
  output$plot11 <- renderPlot({
    if(input$check_JH03){
      inputdf<-data.frame(genesql())
      gene_name<-inputdf[1,1]
      pt_size<-input$pt_size
      cells<-rownames(cellData)[cellData$orig_ident==conditions[3]]
      p<-activityUmap(inputdf = inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = conditions[3],subset=FALSE)
      print(p)}
  })
  output$plot12 <- renderPlot({
    if(input$check_JH04){
      inputdf<-data.frame(genesql())
      gene_name<-inputdf[1,1]
      pt_size<-input$pt_size
      cells<-rownames(cellData)[cellData$orig_ident==conditions[4]]
      p<-activityUmap(inputdf = inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = conditions[4],subset=FALSE)
      print(p)}
  })
  
  ############################ end 'Plots' ######################
  
  ##############################################For 'Neutrophil Plots'##################################################  
  # Generate a reactive list of plots that changes when the gene changes.
  
  output$plot13 <- renderPlot({
    # print(input$vln_celltypes)
    inputdf<-data.frame(genesql3())
    gene_name<-inputdf[1,1]
    inputdf<-inputdf[inputdf$cell_name %in% rownames(cellData2),]
    print(dim(inputdf))
    pt_size<-input$pt_size2
    cells<-inputdf[,2]
    p<-expressionUmap(inputdf=inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = "all",subset=TRUE)
    print(p)
  })
  
  output$plot14 <- renderPlot({
    inputdf<-data.frame(genesql3())
    gene_name<-inputdf[1,1]
    inputdf<-inputdf[inputdf$cell_name %in% rownames(cellData2),]
    p<-expressionDot(inputdf = inputdf,gene_name = gene_name,subset=TRUE)
    print(p)
    
  })
  
  
  output$plot15 <- renderPlot({
    inputdf<-data.frame(genesql3())
    gene_name<-inputdf[1,1]
    inputdf<-inputdf[inputdf$cell_name %in% rownames(cellData2),]
    pt_size<-input$pt_size2
    cells<-inputdf[,2]
    p<-activityUmap(inputdf = inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = "all",subset=TRUE)
    print(p)
  })
  
  output$plot16 <- renderPlot({
    inputdf<-data.frame(genesql3())
    gene_name<-inputdf[1,1]
    inputdf<-inputdf[inputdf$cell_name %in% rownames(cellData2),]
    p<-activityDot(inputdf=inputdf,gene_name = gene_name,subset=TRUE)
    print(p)
    
  })
  
  output$plotvln3 <- renderPlot({
    
    
    gene_expression_df<-data.frame(genesql3())
    gene<-gene_expression_df[1,1]
    gene_expression_df<-gene_expression_df[gene_expression_df$cell_name %in% rownames(cellData2),]
    gene_expression_df$cellTypes <- cellData2$celltype
    gene_expression_df$obj.orig.ident <- cellData2$orig_ident
    gene_expression_df$Cell_Groups <- paste0(gene_expression_df$cellTypes," ",gene_expression_df$obj.orig.ident)
    set.seed(123)  # Set seed for reproducibility
    noise <- rnorm(n = nrow(gene_expression_df)) / 100000
    gene_expression_df[,3] <- gene_expression_df[,3] + noise
    colnames(gene_expression_df)[3]<-"expression"
    
    color_palette <- c("red", "blue", "green", "purple","orange","grey")
    
    p<-ggplot(gene_expression_df, aes(x = obj.orig.ident, y = expression, fill = obj.orig.ident)) +
      geom_violin(trim = TRUE, scale = "width") +
      geom_jitter(size=input$pt_size) +
      scale_fill_manual(values = setNames(color_palette, unique(gene_expression_df$obj.orig.ident))) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      ggtitle(paste0(gene," expression")) +
      labs(x = "Condition Groups")+
      labs(y = paste0(gene))+
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        plot.margin = unit(rep(0, 4), "cm"),
        plot.title = element_text(hjust = 0.5),  # Center the ggtitle
        legend.key.width = unit(0.5, "cm"),  # Adjust the width of the legend key
        legend.key.height = unit(0.5, "cm")  # Adjust the height of the legend key
      )
    
    print(p)
  })
  
  output$plotvln4 <- renderPlot({
    cellsChosen<-as.character(input$vln_split2)
    
    gene_expression_df<-data.frame(genesql3())
    gene<-gene_expression_df[1,1]
    gene_expression_df<-gene_expression_df[gene_expression_df$cell_name %in% rownames(cellData2),]
    gene_expression_df$cellTypes <- cellData2$celltype
    gene_expression_df$obj.orig.ident <- cellData2$orig_ident
    gene_expression_df$Cell_Groups <- paste0(gene_expression_df$cellTypes," ",gene_expression_df$obj.orig.ident)
    gene_expression_df<-gene_expression_df[gene_expression_df$obj.orig.ident %in% cellsChosen,]
    gene_expression_df<-gene_expression_df[gene_expression_df$cellTypes %in% cellsChosen,]
    set.seed(123)  # Set seed for reproducibility
    noise <- rnorm(n = nrow(gene_expression_df)) / 100000
    gene_expression_df[,3] <- gene_expression_df[,4] + noise
    colnames(gene_expression_df)[3]<-"expression"
    
    color_palette <- c("red", "blue", "green", "purple","orange","grey")
    gene_expression_df$Cell_Groups<-factor(gene_expression_df$Cell_Groups,levels = sort(unique(gene_expression_df$Cell_Groups)))
    p<-ggplot(gene_expression_df, aes(x = Cell_Groups, y = expression, fill = obj.orig.ident)) +
      geom_violin(trim = TRUE, scale = "width") +
      geom_jitter(size=input$pt_size) +
      scale_fill_manual(values = setNames(color_palette, unique(gene_expression_df$obj.orig.ident))) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      ggtitle(paste0(gene," expression")) +
      labs(x = "Condition Groups")+
      labs(y = paste0(gene))+
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        plot.margin = unit(rep(0, 4), "cm"),
        plot.title = element_text(hjust = 0.5),  # Center the ggtitle
        legend.key.width = unit(0.5, "cm"),  # Adjust the width of the legend key
        legend.key.height = unit(0.5, "cm")  # Adjust the height of the legend key
      )
    
    print(p)
  })
  #expression subsets
  output$plot17 <- renderPlot({
    if(input$check_JH012){
      inputdf<-data.frame(genesql3())
      gene_name<-inputdf[1,1]
      inputdf<-inputdf[inputdf$cell_name %in% rownames(cellData2),]
      pt_size<-input$pt_size2
      cells<-rownames(cellData2)[cellData2$orig_ident==conditions[1]]
      p<-expressionUmap(inputdf = inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = conditions[1],subset=TRUE)
      print(p)}
  })
  output$plot18 <- renderPlot({
    if(input$check_JH022){
      inputdf<-data.frame(genesql3())
      gene_name<-inputdf[1,1]
      inputdf<-inputdf[inputdf$cell_name %in% rownames(cellData2),]
      pt_size<-input$pt_size2
      cells<-rownames(cellData)[cellData2$orig_ident==conditions[2]]
      p<-expressionUmap(inputdf = inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = conditions[2],subset=TRUE)
      print(p)}
  })
  output$plot19 <- renderPlot({
    if(input$check_JH032){
      inputdf<-data.frame(genesql3())
      gene_name<-inputdf[1,1]
      inputdf<-inputdf[inputdf$cell_name %in% rownames(cellData2),]
      pt_size<-input$pt_size2
      cells<-rownames(cellData)[cellData2$orig_ident==conditions[3]]
      p<-expressionUmap(inputdf = inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = conditions[3],subset=TRUE)
      print(p)}
  })
  output$plot20 <- renderPlot({
    if(input$check_JH042){
      inputdf<-data.frame(genesql3())
      gene_name<-inputdf[1,1]
      inputdf<-inputdf[inputdf$cell_name %in% rownames(cellData2),]
      pt_size<-input$pt_size2
      cells<-rownames(cellData)[cellData2$orig_ident==conditions[4]]
      p<-expressionUmap(inputdf = inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = conditions[4],subset=TRUE)
      print(p)}
  })
  
  
  #activity subsets
  output$plot21 <- renderPlot({
    if(input$check_JH012){
      inputdf<-data.frame(genesql3())
      gene_name<-inputdf[1,1]
      inputdf<-inputdf[inputdf$cell_name %in% rownames(cellData2),]
      pt_size<-input$pt_size2
      cells<-rownames(cellData)[cellData2$orig_ident==conditions[1]]
      p<-activityUmap(inputdf = inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = conditions[1],subset=TRUE)
      print(p)}
  })
  output$plot22 <- renderPlot({
    if(input$check_JH022){
      inputdf<-data.frame(genesql3())
      gene_name<-inputdf[1,1]
      inputdf<-inputdf[inputdf$cell_name %in% rownames(cellData2),]
      pt_size<-input$pt_size2
      cells<-rownames(cellData)[cellData2$orig_ident==conditions[2]]
      p<-activityUmap(inputdf = inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = conditions[2],subset=TRUE)
      print(p)}
  })
  output$plot23 <- renderPlot({
    if(input$check_JH032){
      inputdf<-data.frame(genesql3())
      gene_name<-inputdf[1,1]
      inputdf<-inputdf[inputdf$cell_name %in% rownames(cellData2),]
      pt_size<-input$pt_size2
      cells<-rownames(cellData)[cellData2$orig_ident==conditions[3]]
      p<-activityUmap(inputdf = inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = conditions[3],subset=TRUE)
      print(p)}
  })
  output$plot24 <- renderPlot({
    if(input$check_JH042){
      inputdf<-data.frame(genesql3())
      gene_name<-inputdf[1,1]
      inputdf<-inputdf[inputdf$cell_name %in% rownames(cellData2),]
      pt_size<-input$pt_size2
      cells<-rownames(cellData)[cellData2$orig_ident==conditions[4]]
      p<-activityUmap(inputdf = inputdf,gene_name = gene_name,pt_size = pt_size,cells = cells,sampleName = conditions[4],subset=TRUE)
      print(p)}
  })
  
  ############################ end 'Neutrophil Plots' ######################
}

# Run the application
shinyApp(ui = ui, server = server)
