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
library(ggrepel)
library(forcats)

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

uniprot_path <- file.path(
  data_dir,
  "uniprotkb_proteome_UP000000625_2026_03_19.tsv"
)

# ---------------------------
# Fail fast if data missing
# ---------------------------
required_files <- c(cor_mat_path, supp_table_path, uniprot_path)
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

uniprot_to_function <- readr::read_tsv(
  uniprot_path,
  show_col_types = FALSE
)

# ---------------------------
# Build GO lookup table once
# ---------------------------
unique_GO <- uniprot_to_function %>%
  tidyr::separate_rows(`Gene Ontology (biological process)`, sep = ";\\s*") %>%
  dplyr::pull(`Gene Ontology (biological process)`) %>%
  unique()

go_df <- tibble(go_term = unique_GO) %>%
  filter(!is.na(go_term), go_term != "") %>%
  mutate(
    go_term = str_remove(go_term, "\\s*\\[GO:\\d+\\]$"),
    go_group = case_when(
      
      # --- central carbon metabolism ---
      str_detect(go_term, "glycolysis|gluconeogenesis|pentose phosphate|pentose-phosphate|TCA|tricarboxylic acid|fermentation|carbon fixation|pyruvate|acetate|acetyl-CoA|propionate|succinyl-CoA|2-oxoglutarate|glyoxylate|succinate|citrate|malate|fumarate|lactate|ethanol|acetaldehyde|glycerol|generation of precursor metabolites") ~ "Central carbon metabolism",
      
      # --- secondary carbon utilization ---
      str_detect(go_term, "carbon utilization|glucose|galactose|ribose|fucose|arabinose|maltose|lactose|mannitol|fructose|carbohydrate|hexose|monosaccharide|oligosaccharide|polysaccharide|glucarate|glucuronate|galactarate|galactonate|galacturonate|trehalose|cellobiose|xylan|pectin|gluconate|glucan|glycogen|mannose|rhamnose|sorbitol|melibiose|inositol|allantoin|ascorbic acid|tartrate|galactitol|allose|glycolate|ketone|carboxylic acid|organic acid|dicarboxylic acid|monocarboxylic acid|aldonic acid|organic phosphonate|aminophosphonate|N-acetylglucosamine|N-acetylmuramic acid|N-acetylmannosamine|N-acetylneuraminate|diacetylchitobiose|chitin|xylose|lyxose|xylulose|indole|methylglyoxal|formaldehyde|oxalate|phenylpropanoid|bile acid|carnitine|ureide|acetoacetic acid|butanol|dimethyl sulfoxide|amino sugar|cyanate") ~ "Secondary carbon utilization",
      
      # --- nucleotide metabolism ---
      str_detect(go_term, "purine|pyrimidine|nucleotide|nucleoside|nucleobase|adenosine|guanosine|inosine|xanthosine|cytosine|uracil|uridine|thymine|thymidine|adenine|guanine|hypoxanthine|UMP|GMP|IMP|XMP|CMP|TMP|TTP|dATP|dTTP|UTP|CTP|ADP|dTMP|salvage|deoxycytidine|UDP|TDP|GTP|AMP") ~ "Nucleotide metabolism",
      
      # --- protein homeostasis ---
      str_detect(go_term, "protein catabolic|polymerization|stabilization|destabilization|folding|refolding|assembly|targeting|insertion|secretion|trimerization|tetramerization|oligomerization|lipoylation|glycosylation|methylation|adenylylation|proteolysis|peptide catabolic|peptide metabolic|chaperone|zymogen activation|oligomerization|trimerization|tetramerization|hexamerization|protein complex assembly|protein processing|signal peptide processing|glycosylation|lipoylation|protein repair|protein maturation|protein unfolding|protein unfolding|protein quality|protein autoprocessing|protein localization|denaturation|unfolded protein|protein de") ~ "Protein homeostasis",
      
      # --- energy metabolism ---
      str_detect(go_term, "electron transport chain|respiration|oxidative phosphorylation|ATP synthesis|hydrogen metabolic process|respiratory chain|generation of precursor metabolites and energy") ~ "Energy metabolism",
      
      # --- cofactor / vitamin metabolism ---
      str_detect(go_term, "heme|porphyrin|cytochrome|quinone|ubiquinone|menaquinone|NAD|NADP|FAD|FMN|lipoate|molybdopterin|vitamin|cofactor|prosthetic group|siderophore|enterobactin|thiamine|folic acid|tetrahydrofolate|dihydrofolate|pteridine|cobalamin|riboflavin|pyridox|coenzyme A|pantothenate|isoprenoid|terpenoid|polyprenol|ferredoxin|glutathione") ~ "Cofactor / vitamin metabolism",
      
      # --- ion / homeostasis ---
      str_detect(go_term, "magnesium|iron|zinc|copper|manganese|nickel|cadmium|cobalt|selenium|molybdate|potassium|phosphate|phosphorus|sulfate|sulfur compound|ion homeostasis|membrane potential|cation|homeostasis") ~ "Ion homeostasis / ion response",
      
      # --- transport ---
      str_detect(go_term, "transport|import|export") ~ "Transport",
      
      # --- transcription ---
      str_detect(go_term, "transcription|gene expression") ~ "Transcription",
      
      # --- DNA maintenance ---
      str_detect(go_term, "DNA replication|DNA repair|recombination|translesion|SOS response|chromosome|DNA damage|DNA integration|transposition|transformation|base-excision repair|double-strand break repair|interstrand cross-link repair|mismatch|plasmid partitioning|sister chromatid cohesion|DNA protection|CRISPR|DNA|photoreactive repair|methylation|plasmid") ~ "DNA maintenance",
      
      # --- translation / RNA biology ---
      str_detect(go_term, "translation|ribosome|ribosomal|tRNA|rRNA|mRNA|RNA processing|RNA modification|RNA methylation|RNA capping|ncRNA|pseudouridine|endoribonuclease|RNA metabolic process|RNA") ~ "Translation / RNA biology",
      
      # --- motility ---
      str_detect(go_term, "flagell|chemotaxis|motility|aerotaxis|thermotaxis|pilus") ~ "Motility / taxis",
      
      # --- cell envelope ---
      str_detect(go_term, "cell wall|peptidoglycan|outer membrane|lipopolysaccharide|capsule|membrane assembly|membrane organization|periplasmic space organization|lipoprotein|cellulose|biofilm|cell shape|cell morphogenesis|cell growth|cell size|division|septum|cytokinesis|curli|amyloid fibril|O antigen|common antigen|adhesion|cell killing|autolysis|ring assembly|growth|osmosensory|cytoskeleton|invagination|cell tip|colanic acid") ~ "Cell envelope",
      
      # --- lipid metabolism ---
      str_detect(go_term, "fatty acid|lipid|phospholipid|cardiolipin|lipid A|phosphatid|acyl-CoA|butyryl-CoA|poly-hydroxybutyrate|steroid") ~ "Lipid metabolism",
      
      # --- amino acid metabolism ---
      str_detect(go_term, "amino acid|arginine|lysine|methionine|threonine|serine|glycine|tryptophan|tyrosine|histidine|cysteine|glutamate|aspartate|proline|alanine|valine|leucine|isoleucine|phenylalanine|ornithine|putrescine|spermidine|spermine|taurine|GABA|asparagine|glutamine|amine metabolic|biogenic amine|sulfur utilization|gamma-aminobutyric acid|nitrate|nitrite|ammonia|urea|amine|nitrogen utilization") ~ "Amino acid metabolism",
      
      
      # --- stress / response ---
      str_detect(go_term, "stress|response to|detoxification|starvation|oxidative|heat|cold|temperature|osmotic|hypotonic|hyperosmotic|acidic pH|alkaline pH|radiation|antibiotic|xenobiotic|toxic substance|nitric oxide|hydrogen peroxide|superoxide|defense|phage shock|dormancy|programmed cell death|symbiotic interaction|virus|viral|toxin|stringent") ~ "Stress / environmental response",
      
      # --- signaling ---
      str_detect(go_term, "signal transduction|phosphorelay|phosphorylation|autophosphorylation|cAMP|quorum sensing|sensory|detection") ~ "Signal transduction / signaling",
      
      
      TRUE ~ "Other"
    )
  )

# ---------------------------
# Large matrix will be loaded lazily in server()
# ---------------------------

# ---------------------------
# GO helper functions
# ---------------------------
build_go_analysis_df <- function(corr_df, uniprot_to_function, go_df) {
  
  uniprot_small <- uniprot_to_function %>%
    transmute(
      Gene = `Gene Names (primary)`,
      go_raw = `Gene Ontology (biological process)`
    ) %>%
    filter(!is.na(Gene), Gene != "")
  
  corr_df %>%
    left_join(uniprot_small, by = "Gene") %>%
    separate_rows(go_raw, sep = ";\\s*") %>%
    mutate(
      go_term = stringr::str_trim(stringr::str_remove(go_raw, "\\s*\\[GO:\\d+\\]$")),
      go_term = na_if(go_term, "")
    ) %>%
    left_join(go_df, by = "go_term") %>%
    mutate(go_group = dplyr::coalesce(go_group, "Other"))
}

annotate_with_go <- function(df, uniprot_to_function, go_df) {
  df %>%
    left_join(
      uniprot_to_function %>%
        transmute(
          Gene = `Gene Names (primary)`,
          go_raw = `Gene Ontology (biological process)`
        ),
      by = "Gene"
    ) %>%
    separate_rows(go_raw, sep = ";\\s*") %>%
    mutate(
      go_term = stringr::str_remove(go_raw, "\\s*\\[GO:\\d+\\]$"),
      go_term = stringr::str_trim(go_term)
    ) %>%
    left_join(go_df, by = "go_term")
}

goAnalysisModuleUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(
        3,
        selectInput(
          ns("analysis_mode"),
          "Analysis mode",
          choices = c(
            "Continuous" = "continuous",
            "Threshold enrichment" = "threshold"
          ),
          selected = "continuous"
        )
      ),
      column(
        3,
        numericInput(
          ns("threshold"),
          "GO threshold",
          value = 0.25,
          step = 0.01
        )
      ),
      column(
        3,
        numericInput(
          ns("min_group_size"),
          "Minimum GO group size",
          value = 5,
          min = 1,
          step = 1
        )
      )),
    fluidRow(
      column(
        3,
        checkboxInput(ns("hide_translation"), "Hide Translation / RNA biology", value = FALSE)
      ),
      column(
        3,
        checkboxInput(ns("drop_other"), "Hide Other", value = FALSE)
      )
    ),
  
  tabsetPanel(
      tabPanel("GO summary", DTOutput(ns("summary_table"))),
      tabPanel("GO distribution", plotOutput(ns("distribution_plot"), height = "650px")),
      tabPanel("GO main", plotOutput(ns("main_plot"), height = "650px")),
      tabPanel("GO quadrant", plotOutput(ns("quadrant_plot"), height = "650px"))
    )
  )
}

goAnalysisModuleServer <- function(id, go_df, corr_df, uniprot_to_function) {
  moduleServer(id, function(input, output, session) {
    
    annotated_df <- reactive({
      req(go_df(), corr_df(), uniprot_to_function())
      build_go_analysis_df(
        corr_df = corr_df(),
        uniprot_to_function = uniprot_to_function(),
        go_df = go_df()
      )
    })
    
    filtered_df <- reactive({
      df <- annotated_df() %>%
        filter(!is.na(Correlation_with_geneX), !is.na(go_group), go_group != "")
      
      if (isTRUE(input$drop_other)) {
        df <- df %>% filter(go_group != "Other")
      }
      
      if (isTRUE(input$hide_translation)) {
        df <- df %>% filter(go_group != "Translation / RNA biology")
      }
      
      df
    })
    
    summary_df <- reactive({
      df <- filtered_df()
      req(nrow(df) > 0)
      
      global_median <- median(df$Correlation_with_geneX, na.rm = TRUE)
      global_prop_high <- mean(df$Correlation_with_geneX > input$threshold, na.rm = TRUE)
      
      df %>%
        group_by(go_group) %>%
        summarise(
          n_total = n(),
          mean_corr = mean(Correlation_with_geneX, na.rm = TRUE),
          median_corr = median(Correlation_with_geneX, na.rm = TRUE),
          n_high = sum(Correlation_with_geneX > input$threshold, na.rm = TRUE),
          prop_high = mean(Correlation_with_geneX > input$threshold, na.rm = TRUE),
          expected_high = n_total * global_prop_high,
          enrichment = ifelse(expected_high > 0, n_high / expected_high, NA_real_),
          delta_median = median_corr - global_median,
          .groups = "drop"
        ) %>%
        filter(n_total >= input$min_group_size)
    })
    
    output$summary_table <- DT::renderDT({
      req(summary_df())
      
      DT::datatable(
        summary_df() %>% arrange(desc(delta_median)),
        filter = "top",
        options = list(pageLength = 15, scrollX = TRUE)
      )
    })
    
    output$distribution_plot <- renderPlot({
      df <- filtered_df() %>%
        group_by(go_group) %>%
        mutate(n_group = n()) %>%
        ungroup() %>%
        filter(n_group >= input$min_group_size)
      
      req(nrow(df) > 0)
      
      ggplot(
        df,
        aes(
          x = forcats::fct_reorder(go_group, Correlation_with_geneX, .fun = median, na.rm = TRUE),
          y = Correlation_with_geneX
        )
      ) +
        geom_boxplot(outlier.alpha = 0.25) +
        coord_flip() +
        theme_minimal() +
        labs(x = "GO group", y = "Correlation")
    })
    
    output$main_plot <- renderPlot({
      df <- summary_df()
      req(nrow(df) > 0)
      
      if (input$analysis_mode == "continuous") {
        ggplot(df, aes(x = reorder(go_group, delta_median), y = delta_median, fill = delta_median)) +
          geom_col() +
          coord_flip() +
          geom_hline(yintercept = 0, linetype = 2) +
          scale_fill_gradient2(midpoint = 0) +
          theme_minimal() +
          labs(
            x = "GO group",
            y = "Median correlation - global median",
            fill = "Shift"
          )
      } else {
        ggplot(df, aes(x = reorder(go_group, enrichment), y = enrichment, fill = enrichment)) +
          geom_col() +
          coord_flip() +
          geom_hline(yintercept = 1, linetype = 2) +
          scale_fill_gradient2(midpoint = 1) +
          theme_minimal() +
          labs(
            x = "GO group",
            y = "Observed / Expected",
            fill = "Enrichment"
          )
      }
    })
    
    output$quadrant_plot <- renderPlot({
      df <- summary_df()
      req(nrow(df) > 0)
      
      ggplot(df, aes(x = enrichment, y = delta_median, size = n_total, label = go_group)) +
        geom_point(alpha = 0.7) +
        geom_vline(xintercept = 1, linetype = 2) +
        geom_hline(yintercept = 0, linetype = 2) +
        ggrepel::geom_text_repel(max.overlaps = 25) +
        theme_minimal() +
        labs(
          x = "Enrichment (Observed / Expected)",
          y = "Median correlation - global median",
          size = "Group size"
        )
    })
    
    list(
      annotated_df = annotated_df,
      summary_df = summary_df
    )
  })
}

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
        "RNAseq Data from Supplementary Data 2 of Tjaden, B. (2023)... ",
        "https://doi.org/10.1080/15476286.2023.2189331"
      ),
      h6(
        "GO terms from Uniprot GO annotation, sorted into 17 subcategories (see README)"
      ),
      
      conditionalPanel(
        condition = "!input.batch",
        h3("Correlation Heatmap"),
        uiOutput("heatmap_svg_ui"),
        downloadButton("download_heatmap", "Download Heatmap (SVG)"),
        h3("Top Correlated Genes"),
        DTOutput("corTable"),
        downloadButton("download_table", "Download Table (CSV)"),
        h3("GO Group Analysis"),
        goAnalysisModuleUI("go_mod")
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
    full_table = correlation_df,
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
  full_correlation_rv <- reactiveVal(NULL)
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
        show_numbers  = input$heatmap_num
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
        full_correlation_rv(res$full_table)
        
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
  # GO module wiring
  # ---------------------------
  go_df_r <- reactive({ go_df })
  uniprot_r <- reactive({ uniprot_to_function })
  corr_for_go_r <- reactive({ full_correlation_rv() })
  
  goAnalysisModuleServer(
    id = "go_mod",
    go_df = go_df_r,
    corr_df = corr_for_go_r,
    uniprot_to_function = uniprot_r
  )
  
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
      write.csv(significant_hits_rv(), file, row.names = FALSE)
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
        
        # ---- add GO annotation ----
        annotated <- annotate_with_go(
          res$table,
          uniprot_to_function,
          go_df
        )
        
        write.csv(
          annotated,
          file.path(tables_dir, paste0(g, "_correlations_GO.csv")),
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
  )
}

shinyApp(ui = ui, server = server)
