# Gene summary card (MyGene). Renders the shared resolved-gene reactive.

gene_summary_ui <- function(id) {
  ns <- NS(id)
  card(
    card_header("Gene"),
    card_body(shinycssloaders::withSpinner(
      uiOutput(ns("content")),
      proxy.height = "120px"
    ))
  )
}

# resolved: reactive() returning the mygene_resolve() result (or NULL before a
# search has run).
gene_summary_server <- function(id, resolved) {
  moduleServer(id, function(input, output, session) {
    output$content <- renderUI({
      res <- resolved()
      if (is.null(res)) {
        return(vr_empty("Search for a gene to see its summary."))
      }
      if (!isTRUE(res$ok)) {
        return(vr_error(res$error))
      }
      tagList(
        tags$h5(
          res$symbol,
          tags$small(class = "text-muted", paste0("  ", res$name))
        ),
        vr_field("Type", res$type_of_gene),
        vr_field("Entrez", res$entrez),
        vr_field("Ensembl", res$ensembl_gene),
        vr_field("UniProt", res$uniprot),
        if (!is_blank(res$summary)) {
          tags$p(class = "mt-2 mb-0", res$summary)
        }
      )
    })
  })
}
