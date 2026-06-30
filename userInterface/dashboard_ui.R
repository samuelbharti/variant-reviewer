# Single-page dashboard: search bar on top, then a grid of result cards.

dashboard_page <- tagList(
  gene_search_ui("search"),
  layout_columns(
    col_widths = c(4, 4, 4),
    gene_summary_ui("gene_summary"),
    variant_summary_ui("variant_summary"),
    protein_summary_ui("protein_summary")
  ),
  layout_columns(
    col_widths = c(6, 6),
    clinvar_ui("clinvar"),
    gnomad_ui("gnomad")
  ),
  gtex_expression_ui("gtex"),
  layout_columns(
    col_widths = c(6, 6),
    string_ppi_ui("string_ppi"),
    opentargets_ui("opentargets")
  ),
  external_links_ui("links")
)
