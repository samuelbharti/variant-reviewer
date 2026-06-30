page_fluid(
  title = "Variant Reviewer",
  # Apply branding from _brand.yml (colors, fonts). brand = TRUE requires the
  # file to exist; switch to bslib::bs_theme() to make it optional.
  theme = bslib::bs_theme(brand = TRUE),
  tags$h2("Variant Reviewer", class = "mt-3"),
  tags$p(
    class = "text-muted",
    "A lightweight gene and variant interpretation companion."
  ),
  dashboard_page
)
