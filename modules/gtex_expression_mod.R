# GTEx tissue expression card. Plots median TPM across tissues for the gene.

gtex_expression_ui <- function(id) {
  ns <- NS(id)
  card(
    card_header("Tissue expression (GTEx)"),
    card_body(shinycssloaders::withSpinner(
      uiOutput(ns("content")),
      proxy.height = "300px"
    ))
  )
}

# resolved: reactive() -> mygene_resolve() result (uses the gene symbol).
gtex_expression_server <- function(id, resolved) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    expression <- reactive({
      res <- resolved()
      if (is.null(res) || !isTRUE(res$ok)) {
        return(NULL)
      }
      gtex_median_expression(res$symbol)
    })

    output$plot <- renderPlot({
      res <- expression()
      req(res, isTRUE(res$ok))
      df <- res$data
      df$tissue <- factor(df$tissue, levels = df$tissue[order(df$median_tpm)])
      ggplot2::ggplot(df, ggplot2::aes(x = median_tpm, y = tissue)) +
        ggplot2::geom_col(fill = "#0066cc") +
        ggplot2::labs(x = "Median TPM", y = NULL) +
        ggplot2::theme_minimal(base_size = 12)
    })

    output$content <- renderUI({
      res <- expression()
      if (is.null(res)) {
        return(vr_empty("Search for a gene to see tissue expression."))
      }
      if (!isTRUE(res$ok)) {
        return(vr_error(res$error))
      }
      # Scale height to the number of tissues so labels stay legible.
      height <- max(300, nrow(res$data) * 16)
      plotOutput(ns("plot"), height = paste0(height, "px"))
    })
  })
}
