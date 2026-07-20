# Load libraries and source files
library(shiny)
library(bslib)

# bs_theme(brand = TRUE) in ui.R reads _brand.yml via the {brand.yml} package.
# This reference keeps it visible to renv (it is otherwise a phantom dependency)
# and fails fast with a clear message if it is missing.
if (!requireNamespace("brand.yml", quietly = TRUE)) {
  stop("The 'brand.yml' package is required (used by bs_theme(brand = TRUE)).")
}

# Optionally theme base/ggplot/lattice output to match the app theme. Activates
# only if the {thematic} package is installed, so it adds no hard dependency.
if (requireNamespace("thematic", quietly = TRUE)) {
  thematic::thematic_shiny(font = "auto")
}

# Theme reactable tables to match the bslib theme. reactable renders its own
# (light) theme by default, which clashes with the app palette; binding its
# colors to Bootstrap CSS variables makes every table adopt the theme (borders,
# text, header, hover) and track it if the palette changes.
options(
  reactable.theme = reactable::reactableTheme(
    color = "var(--bs-body-color)",
    backgroundColor = "transparent",
    borderColor = "var(--bs-border-color)",
    stripedColor = "var(--bs-tertiary-bg)",
    highlightColor = "var(--bs-secondary-bg)",
    cellPadding = "6px 8px",
    headerStyle = list(
      backgroundColor = "var(--bs-tertiary-bg)",
      borderColor = "var(--bs-border-color)"
    )
  )
)

source("R/load_components.R")

# Load data/connections
# Example: app_data <- readRDS("data/app_data.rds")

# Preprocess small data
