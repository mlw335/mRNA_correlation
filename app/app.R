######
# Gene Correlation Explorer
# Data from Tjaden, B. (2023)
######

# ---------------------------
# Basic libraries
# ---------------------------
library(shiny)
library(tidyverse)
library(pheatmap)
library(DT)
library(svglite)
library(htmltools)
library(zip)

# ---------------------------
# Resolve project paths
# (assumes app is run from project root)
# ---------------------------
project_root <- getwd()
data_dir <- file.path(project_root, "../data")

cor_mat_path <- file.path(
  data_dir,
  "correlation_matrix_all_genes.rds"
)

supp_table_path <- file.path(
  data_dir,
  "Supplementary_Table_2.txt"
)

# ---------------------------
# Fail fast if data missing
# ---------------------------
required_files <- c(cor_mat_path, supp_table_path)
missing <- required_files[!file.exists(required_files)]

if (length(missing) > 0) {
  stop(
    paste(
      "Missing required data files:",
      paste(basename(missing), collapse = ", "),
      "\n\nExpected location:",
      data_dir,
      "\n\nPlease download the data before running the app.",
      sep = "\n"
    ),
    call. = FALSE
  )
}

message("Using user-downloaded local data")

# ---------------------------
# Load small metadata only
# ---------------------------
df <- read.delim(
  supp_table_path,
  header = TRUE,
  fill = TRUE
)

# ---------------------------
# Large matrix will be loaded lazily in server()
# ---------------------------
# UI
# ---------------------------
ui <- fluidPage(
  titlePanel("Gene Correlation Explorer"),
  
  tags$head(
    tags$script(src = "https://unpkg.com/@panzoom/panzoom/dist/panzoom.min.js"),
    tags$style(HTML("
      #heatmap_viewport {
        width: 100%;
        height: 900px;
        border: 1px solid #ddd;
        overflow: hidden;
        position: relative;
        background: white;
      }
      #heatmap_viewport svg {
        width: 100%;
        height: 100%;
        cursor: grab;
      }
      #heatmap_viewport svg:active {
        cursor: grabbing;
      }
      .zoom_controls {
        margin: 8px 0;
        display: flex;
        gap: 8px;
        align-items: center;
      }
    "))
  ),
  
  sidebarLayout(
    sidebarPanel(
      checkboxInput("batch", "Batch input mode", value = FALSE),
      
      conditionalPanel(
        condition = "!input.batch",
        textInput("geneX", "Enter Gene of Interest:", value = "geneX")
      ),
      
      conditionalPanel(
        condition = "input.batch",
        fileInput(
          "batch_file",
          "Upload gene list (one gene per line)",
          accept = c(".txt", ".csv")
        )
      ),
      
      textInput("GOIs", "Additional Genes (comma-separated):", value = ""),
      
      sliderInput(
        "cor_threshold",
        "Correlation cutoff:",
        min = 0, max = 1, value = 0.5, step = 0.05
      ),
      
      numericInput(
        "n_threshold",
        "Minimum n-value (overlaps):",
        value = 100, min = 1
      ),
      
      numericInput(
        "primary_n",
        "Top primary hits:",
        value = 25, min = 1
      ),
      
      checkboxInput(
        "heatmap_num",
        "Include values in heatmap?",
        value = TRUE
      ),
      
      actionButton("run", "Run Analysis")
    ),
    
    mainPanel(
      h6(
        "Data from Tjaden, B. (2023)... ",
        "https://doi.org/10.1080/15476286.2023.2189331"
      ),
      
      conditionalPanel(
        condition = "!input.batch",
        h3("Correlation Heatmap"),
        uiOutput("heatmap_svg_ui"),
        downloadButton("download_heatmap", "Download Heatmap (SVG)"),
        h3("Top Correlated Genes"),
        DTOutput("corTable"),
        downloadButton("download_table", "Download Table (CSV)")
      ),
      
      conditionalPanel(
        condition = "input.batch",
        h3("Batch analysis"),
        downloadButton("download_batch", "Download batch results (ZIP)")
      )
    )
  )
)
  
run_single_gene <- function(
    geneX,
    cor_mat_all,
    df,
    GOIs,
    cor_threshold,
    n_threshold,
    primary_n,
    show_numbers = TRUE
) {
  
  geneX <- trimws(geneX)
  
  # ---- validation ----
  if (!geneX %in% rownames(cor_mat_all)) {
    stop(paste0("Gene '", geneX, "' not found in correlation matrix"))
  }
  
  # ---- reshape expression table ----
  df_long <- df %>%
    pivot_longer(
      cols = -(1:7),
      names_to = "Condition",
      values_to = "Expression"
    ) %>%
    select(Name, Condition, Expression)
  
  df_wide <- df_long %>%
    pivot_wider(names_from = Name, values_from = Expression)
  
  # ---- compute overlaps ----
  df_wide_clean <- df_wide %>%
    filter(!is.na(.data[[geneX]]))
  
  numeric_cols <- df_wide_clean %>%
    select(-Condition) %>%
    mutate(across(
      everything(),
      ~ suppressWarnings(as.numeric(.x))
    ))
  
  cor_with_geneX <- cor_mat_all[geneX, ]
  cor_with_geneX <- cor_with_geneX[names(cor_with_geneX) != geneX]
  
  n_with_geneX <- colSums(
    !is.na(numeric_cols) & !is.na(numeric_cols[[geneX]])
  )
  n_with_geneX <- n_with_geneX[names(n_with_geneX) != geneX]
  
  common_genes <- intersect(
    names(cor_with_geneX),
    names(n_with_geneX)
  )
  
  correlation_df <- tibble(
    Gene = common_genes,
    Correlation_with_geneX = as.numeric(cor_with_geneX[common_genes]),
    n_with_geneX = as.integer(n_with_geneX[common_genes])
  )
  
  # ---- filter significant hits ----
  significant_hits <- correlation_df %>%
    filter(!is.na(Correlation_with_geneX)) %>%
    filter(abs(Correlation_with_geneX) >= cor_threshold) %>%
    filter(n_with_geneX >= n_threshold) %>%
    arrange(desc(abs(Correlation_with_geneX)))
  
  # ---- build heatmap matrix ----
  heatmap_genes <- unique(
    c(geneX, head(significant_hits$Gene, primary_n), GOIs)
  )
  
  heatmap_genes <- heatmap_genes[
    heatmap_genes %in% rownames(cor_mat_all)
  ]
  
  if (length(heatmap_genes) < 2) {
    stop("Not enough genes to construct heatmap")
  }
  
  cor_mat <- cor_mat_all[
    heatmap_genes,
    heatmap_genes,
    drop = FALSE
  ]
  
  # ---- return results ----
  list(
    table = significant_hits,
    heatmap_matrix = cor_mat
  )
}


# ---------------------------
# SERVER
# ---------------------------
server <- function(input, output, session) {
  
  # ---------------------------
  # Reactive storage
  # ---------------------------
  cor_mat_rv <- reactiveVal(NULL)
  heatmap_svg <- reactiveVal(NULL)
  significant_hits_rv <- reactiveVal(NULL)
  batch_results_rv <- reactiveVal(NULL)
  
  # ---------------------------
  # Lazy loader for large matrix
  # ---------------------------
  load_matrix <- function() {
    if (is.null(cor_mat_rv())) {
      message("Loading user-downloaded correlation matrix...")
      cor_mat_rv(readRDS(cor_mat_path))
    }
    cor_mat_rv()
  }
  
  observeEvent(input$run, {
    
    withProgress(message = "Running analysis...", value = 0, {
      
      cor_mat_all <- load_matrix()
      
      GOIs <- strsplit(input$GOIs, ",")[[1]] |> trimws()
      GOIs <- GOIs[GOIs != ""]
      
      params <- list(
        cor_threshold = input$cor_threshold,
        n_threshold   = input$n_threshold,
        primary_n     = input$primary_n,
        show_numbers = input$heatmap_num
      )
      
      if (!input$batch) {
        # --------------------
        # Single-gene mode
        # --------------------
        
        geneX <- trimws(input$geneX)
        
        res <- tryCatch(
          run_single_gene(
            geneX,
            cor_mat_all,
            df,
            GOIs,
            params$cor_threshold,
            params$n_threshold,
            params$primary_n,
            params$show_numbers
          ),
          error = function(e) {
            showNotification(e$message, type = "error")
            return(NULL)
          }
        )
        
        if (is.null(res)) return()
        
        significant_hits_rv(res$table)
        
        # render SVG
        svg_file <- file.path(tempdir(), paste0("heatmap_", Sys.getpid(), ".svg"))
        svglite::svglite(svg_file, width = 18, height = 18)
        pheatmap::pheatmap(
          res$heatmap_matrix,
          clustering_distance_rows = "euclidean",
          clustering_method = "ward.D2",
          display_numbers = isTRUE(params$show_numbers),
          main = paste0("Neighborhood of ", geneX)
        )
        dev.off()
        
        heatmap_svg(paste(readLines(svg_file), collapse = "\n"))
        
      } else {
        # --------------------
        # Batch mode
        # --------------------
        
        req(input$batch_file)
        
        genes <- readLines(input$batch_file$datapath, warn = FALSE)
        genes <- trimws(genes)
        genes <- genes[genes != ""]
        genes <- unique(genes)
        
        n_genes <- length(genes)
        results <- vector("list", n_genes)
        names(results) <- genes
        
        for (i in seq_along(genes)) {
          g <- genes[i]
          
          incProgress(
            1 / n_genes,
            detail = paste("Processing", g)
          )
          
          results[[g]] <- tryCatch(
            run_single_gene(
              g,
              cor_mat_all,
              df,
              GOIs,
              params$cor_threshold,
              params$n_threshold,
              params$primary_n,
              params$show_numbers
            ),
            error = function(e) NULL
          )
        }
        
        batch_results_rv(results)
      }
    })
  })
  
  # ---------------------------
  # Outputs
  # ---------------------------
  output$corTable <- DT::renderDT({
    req(significant_hits_rv())
    DT::datatable(
      significant_hits_rv(),
      filter = "top",
      options = list(pageLength = 25, scrollX = TRUE)
    )
  })
  
  output$heatmap_svg_ui <- renderUI({
    req(heatmap_svg())
    tagList(
      div(id = "heatmap_viewport", HTML(heatmap_svg()))
    )
  })
  
  output$download_heatmap <- downloadHandler(
    filename = function() {
      paste0(
        "correlation_heatmap_",
        trimws(input$geneX),
        "_threshold_", input$cor_threshold,
        ".svg"
      )
    },
    content = function(file) {
      req(heatmap_svg())
      writeLines(heatmap_svg(), file, useBytes = TRUE)
    }
  )
  
  output$download_table <- downloadHandler(
    filename = function() {
      paste0(
        "correlated_genes_",
        trimws(input$geneX),
        "_threshold_", input$cor_threshold,
        ".csv"
      )
    },
    content = function(file) {
      req(significant_hits_rv())
      write.csv(significant_hits_rv(),
                file, row.names = FALSE)
    }
  )
  
  output$download_batch <- downloadHandler(
    filename = function() {
      paste0("GeneCorrelationExplorer_batch_", Sys.Date(), ".zip")
    },
    content = function(file) {
      
      results <- batch_results_rv()
      req(results)
      
      tmpdir <- tempdir()
      tables_dir   <- file.path(tmpdir, "tables")
      heatmaps_dir <- file.path(tmpdir, "heatmaps")
      
      dir.create(tables_dir, showWarnings = FALSE)
      dir.create(heatmaps_dir, showWarnings = FALSE)
      
      for (g in names(results)) {
        res <- results[[g]]
        if (is.null(res)) next
        
        write.csv(
          res$table,
          file.path(tables_dir, paste0(g, "_correlations.csv")),
          row.names = FALSE
        )
        
        svglite::svglite(
          file.path(heatmaps_dir, paste0(g, "_heatmap.svg")),
          width = 8,
          height = 8
        )
        pheatmap::pheatmap(res$heatmap_matrix)
        dev.off()
      }
      
      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(tmpdir)
      
      zip::zip(
        zipfile = file,
        files = c(
          file.path("tables", list.files("tables")),
          file.path("heatmaps", list.files("heatmaps"))
        )
      )
    }
  )}

shinyApp(ui = ui, server = server)
