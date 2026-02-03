# -----------------------------
# Gene Correlation Explorer
# One-click launcher
# -----------------------------


# ---------------------------
# Package bootstrap
# ---------------------------
required_packages <- c(
  "shiny",
  "tidyverse",
  "pheatmap",
  "DT",
  "svglite",
  "htmltools"
)

installed <- rownames(installed.packages())

missing <- setdiff(required_packages, installed)

if (length(missing) > 0) {
  message("Installing missing packages: ",
          paste(missing, collapse = ", "))
  install.packages(missing, dependencies = TRUE)
}

invisible(lapply(required_packages, library, character.only = TRUE))

message("Launching Gene Correlation Explorer...")

# Fail early if app directory is missing
if (!dir.exists("app")) {
  stop(
    "Could not find 'app/' directory.\n",
    "Please run this script from the project root.",
    call. = FALSE
  )
}

# Launch Shiny app

shiny::runApp("app", launch.browser = TRUE)
