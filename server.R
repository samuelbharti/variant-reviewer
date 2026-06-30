# Shiny Server
function(input, output, session) {
  # Submitted search query: reactive(list(gene, variant)) or NULL.
  search <- gene_search_server("search")

  # Resolve the gene once and share its identifiers with every gene-level
  # module (summary, GTEx, STRING, protein, links).
  resolved <- reactive({
    query <- search()
    if (is.null(query)) {
      return(NULL)
    }
    mygene_resolve(query$gene)
  })

  # Annotate the variant once and share it with the variant and protein cards
  # (both need it), so MyVariant is only queried a single time per search.
  variant_annotation <- reactive({
    query <- search()
    if (is.null(query) || is_blank(query$variant)) {
      return(NULL)
    }
    myvariant_annotate(query$variant)
  })

  gene_summary_server("gene_summary", resolved)
  variant_summary_server("variant_summary", variant_annotation)
  protein_summary_server(
    "protein_summary",
    resolved,
    search,
    variant_annotation
  )
  gtex_expression_server("gtex", resolved)
  string_ppi_server("string_ppi", resolved)
  external_links_server("links", resolved)
}
