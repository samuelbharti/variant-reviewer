# Gene-first search box. Returns a reactive carrying the submitted query so the
# parent can fan it out to the result modules.

gene_search_ui <- function(id) {
  ns <- NS(id)

  card(
    card_body(
      layout_columns(
        col_widths = c(5, 4, 3),
        textInput(
          ns("gene"),
          label = "Gene symbol",
          placeholder = "e.g. TP53",
          width = "100%"
        ),
        textInput(
          ns("variant"),
          label = "Variant (optional)",
          placeholder = "e.g. R175H or rs113488022",
          width = "100%"
        ),
        div(
          class = "d-grid align-self-end mb-1",
          actionButton(
            ns("submit"),
            label = "Review",
            icon = icon("magnifying-glass"),
            class = "btn-primary"
          )
        )
      )
    )
  )
}

# Returns a reactive carrying list(gene = <chr>, variant = <chr|NULL>), or NULL
# before the first submit (and when the gene is blank). A reactiveVal is used
# instead of eventReactive so reading it before any submit yields NULL rather
# than a silent error, which lets the result cards show their initial
# placeholder messages.
gene_search_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    query <- reactiveVal(NULL)
    observeEvent(input$submit, {
      gene <- trimws(input$gene %||% "")
      if (gene == "") {
        query(NULL)
        return()
      }
      variant <- trimws(input$variant %||% "")
      query(list(
        gene = gene,
        variant = if (variant == "") NULL else variant
      ))
    })
    query
  })
}
