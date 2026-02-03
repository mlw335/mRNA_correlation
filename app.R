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
      textInput("geneX", "Enter Gene of Interest:", value = "geneX"),
      textInput("GOIs", "Additional Genes (comma-separated):", value = ""),
      sliderInput("cor_threshold", "Correlation cutoff:",
                  min = 0, max = 1, value = 0.5, step = 0.05),
      numericInput("n_threshold", "Minimum n-value (overlaps):",
                   value = 100, min = 1),
      numericInput("primary_n", "Top primary hits:",
                   value = 25, min = 1),
      checkboxInput("heatmap_num", "Include Values in Heatmap?",
                    value = TRUE),
      
      actionButton("run", "Run Analysis"),
      
      div(class = "zoom_controls",
          actionButton("zoom_in", "+"),
          actionButton("zoom_out", "−"),
          actionButton("zoom_reset", "reset")
      )
    ),
    
    mainPanel(
      h6("Data from Tjaden, B. (2023)... https://doi.org/10.1080/15476286.2023.2189331"),
      h3("Correlation Heatmap (pheatmap → SVG, zoomable)"),
      uiOutput("heatmap_svg_ui"),
      downloadButton("download_heatmap", "Download Heatmap (SVG)"),
      h3("Top Correlated Genes"),
      DTOutput("corTable"),
      downloadButton("download_table", "Download Table (CSV)")
    )
  )
)

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
      
      # ---- load matrix only when needed ----
      cor_mat_all <- load_matrix()
      
      incProgress(0.1, detail = "Reshaping expression table")
      
      df_long <- df %>%
        pivot_longer(
          cols = -(1:7),
          names_to = "Condition",
          values_to = "Expression"
        ) %>%
        select(Name, Condition, Expression)
      
      df_wide <- df_long %>%
        pivot_wider(names_from = Name, values_from = Expression)
      
      geneX <- trimws(input$geneX)
      GOIs  <- strsplit(input$GOIs, ",")[[1]] |> trimws()
      GOIs  <- GOIs[GOIs != ""]
      
      cor_threshold <- input$cor_threshold
      n_threshold   <- input$n_threshold
      primary_n     <- input$primary_n
      
      # ---- validation ----
      if (!geneX %in% rownames(cor_mat_all)) {
        showNotification(
          paste0("geneX '", geneX, "' not found in correlation matrix."),
          type = "error"
        )
        heatmap_svg(NULL)
        significant_hits_rv(NULL)
        return(NULL)
      }
      
      missing_gois <- setdiff(GOIs, rownames(cor_mat_all))
      if (length(missing_gois) > 0) {
        showNotification(
          paste(
            "Some GOIs not found and will be ignored:",
            paste(missing_gois, collapse = ", ")
          ),
          type = "warning"
        )
        GOIs <- intersect(GOIs, rownames(cor_mat_all))
      }
      
      incProgress(0.2, detail = "Computing n-overlaps")
      
      df_wide_clean <- df_wide %>%
        filter(!is.na(.data[[geneX]]))
      
      numeric_cols <- df_wide_clean %>%
        select(-Condition) %>%
        mutate(across(everything(),
                      ~ suppressWarnings(as.numeric(.x))))
      
      cor_with_geneX_all <- cor_mat_all[geneX, ]
      cor_with_geneX_all <- cor_with_geneX_all[
        names(cor_with_geneX_all) != geneX
      ]
      
      n_with_geneX <- colSums(
        !is.na(numeric_cols) & !is.na(numeric_cols[[geneX]])
      )
      n_with_geneX <- n_with_geneX[
        names(n_with_geneX) != geneX
      ]
      
      common_genes <- intersect(
        names(cor_with_geneX_all),
        names(n_with_geneX)
      )
      
      correlation_df <- tibble(
        Gene = common_genes,
        Correlation_with_geneX = as.numeric(
          cor_with_geneX_all[common_genes]
        ),
        n_with_geneX = as.integer(
          n_with_geneX[common_genes]
        )
      )
      
      incProgress(0.2, detail = "Filtering significant hits")
      
      significant_hits <- correlation_df %>%
        filter(!is.na(Correlation_with_geneX)) %>%
        filter(abs(Correlation_with_geneX) >= cor_threshold) %>%
        filter(n_with_geneX >= n_threshold) %>%
        arrange(desc(abs(Correlation_with_geneX)))
      
      significant_hits_rv(significant_hits)
      
      incProgress(0.2, detail = "Building heatmap matrix")
      
      heatmap_genes <- unique(
        c(geneX, head(significant_hits$Gene, primary_n), GOIs)
      )
      heatmap_genes <- heatmap_genes[
        heatmap_genes %in% rownames(cor_mat_all)
      ]
      
      if (length(heatmap_genes) < 2) {
        heatmap_svg(NULL)
        showNotification("Not enough genes to plot heatmap.",
                         type = "warning")
        return(NULL)
      }
      
      cor_mat <- cor_mat_all[
        heatmap_genes,
        heatmap_genes,
        drop = FALSE
      ]
      
      incProgress(0.2, detail = "Rendering pheatmap to SVG")
      
      svg_file <- file.path(
        tempdir(),
        paste0("heatmap_", Sys.getpid(), ".svg")
      )
      
      svglite::svglite(svg_file, width = 18, height = 18)
      
      pheatmap::pheatmap(
        cor_mat,
        clustering_distance_rows = "euclidean",
        clustering_method = "ward.D2",
        display_numbers = isTRUE(input$heatmap_num),
        main = paste0(
          "Neighborhood of ", geneX,
          "\nCorrelation Threshold = ", cor_threshold
        ),
        fontsize = 10,
        fontsize_number = 8,
        border_color = NA
      )
      
      dev.off()
      
      heatmap_svg(
        paste(readLines(svg_file, warn = FALSE),
              collapse = "\n")
      )
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
}

shinyApp(ui = ui, server = server)
