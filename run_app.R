# -----------------------------
# Gene Correlation Explorer
# One-click launcher
# -----------------------------

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