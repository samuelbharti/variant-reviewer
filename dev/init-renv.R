# dev/init-renv.R
# Helper script to initialize renv for a new project based on this template.

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

# Initialize renv if not already
if (!file.exists("renv.lock")) {
  message("Initializing renv and taking initial snapshot...")
  renv::init(bare = TRUE)
  # Install a few recommended packages used by the template
  renv::install(c(
    "shiny",
    "bslib",
    "brand.yml",
    "ggplot2",
    "httr2",
    "reactable",
    "jsonlite",
    "shinycssloaders",
    "shinytest2",
    "testthat"
  ))
  renv::snapshot()
  message(
    "renv initialized and renv.lock written. Review renv.lock before committing."
  )
} else {
  message("renv.lock already exists. Use renv::restore() to restore packages.")
}
