# Variant summary card (MyVariant). Annotates the optional variant from the
# search query.

variant_summary_ui <- function(id) {
  ns <- NS(id)
  card(
    card_header("Variant"),
    card_body(shinycssloaders::withSpinner(
      uiOutput(ns("content")),
      proxy.height = "120px"
    ))
  )
}

# search: reactive() returning list(gene, variant) or NULL.
variant_summary_server <- function(id, search) {
  moduleServer(id, function(input, output, session) {
    annotation <- reactive({
      query <- search()
      if (is.null(query) || is_blank(query$variant)) {
        return(NULL)
      }
      myvariant_annotate(query$variant)
    })

    output$content <- renderUI({
      res <- annotation()
      if (is.null(res)) {
        return(vr_empty("Enter a variant (rsID or HGVS) to annotate it."))
      }
      if (!isTRUE(res$ok)) {
        return(vr_error(res$error))
      }
      tagList(
        vr_field("Variant", res$id),
        vr_field("dbSNP", res$rsid),
        vr_field("Gene", res$gene),
        vr_field("Protein change", res$hgvsp),
        vr_field(
          "CADD (phred)",
          if (is_blank(res$cadd_phred)) NULL else vr_num(res$cadd_phred, 1)
        ),
        vr_field("ClinVar", res$clinvar_significance)
      )
    })
  })
}
