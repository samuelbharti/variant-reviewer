# Home layout building blocks: the scrollable stack of result cards
# (`dashboard_page`, left column) and the chat panel (`chat_panel`, right
# column). ui.R composes them side by side.

# The result cards the user can show/hide, id -> display label. The id is both
# the module id and the value tracked by the `visible_cards` checkbox group.
.dashboard_cards <- c(
  gene_summary = "Gene summary",
  variant_summary = "Variant",
  predictions = "In-silico predictions",
  protein_summary = "Protein (ProtVar)",
  domains = "Protein domains",
  structure = "3D structure",
  clinvar = "ClinVar",
  gnomad = "gnomAD",
  constraint = "Gene constraint",
  ensembl = "Ensembl VEP",
  gtex = "GTEx expression",
  string_ppi = "STRING",
  opentargets = "Open Targets",
  links = "External links"
)

# Show a card only when its id is ticked in the `visible_cards` control. The
# JS runs client-side, so toggling is instant and needs no server round-trip.
.card_when_shown <- function(id, ui) {
  conditionalPanel(
    condition = sprintf(
      "input.visible_cards && input.visible_cards.indexOf('%s') > -1",
      id
    ),
    ui
  )
}

# Popover to pick which result cards are shown (all ticked by default).
.cards_toggle <- div(
  class = "d-flex justify-content-end mb-2",
  popover(
    actionButton(
      "cards_btn",
      "Cards",
      icon = icon("table-cells"),
      class = "btn-outline-secondary btn-sm"
    ),
    checkboxGroupInput(
      "visible_cards",
      "Show cards",
      choices = stats::setNames(
        names(.dashboard_cards),
        unname(.dashboard_cards)
      ),
      selected = names(.dashboard_cards)
    ),
    title = "Show cards",
    placement = "bottom"
  )
)

# What the annotations are scoped to. Shown under the search so users know the
# species/assembly/transcript basis before reading the cards.
.annotation_note <- tags$p(
  class = "text-muted small mb-2",
  icon("circle-info"),
  paste(
    " Annotating human (Homo sapiens), genome assembly GRCh38 (hg38);",
    "variant consequences are reported per Ensembl transcript (VEP)."
  )
)

# Left column: search bar, a scope note and any gene/variant warning, the card
# picker, then the grid of result cards.
dashboard_page <- tagList(
  gene_search_ui("search"),
  .annotation_note,
  uiOutput("search_notice"),
  .cards_toggle,
  layout_columns(
    col_widths = c(4, 4, 4),
    .card_when_shown("gene_summary", gene_summary_ui("gene_summary")),
    .card_when_shown("variant_summary", variant_summary_ui("variant_summary")),
    .card_when_shown("protein_summary", protein_summary_ui("protein_summary"))
  ),
  layout_columns(
    col_widths = c(6, 6),
    .card_when_shown("clinvar", clinvar_ui("clinvar")),
    .card_when_shown("gnomad", gnomad_ui("gnomad"))
  ),
  layout_columns(
    col_widths = c(6, 6),
    .card_when_shown("predictions", predictions_ui("predictions")),
    .card_when_shown("constraint", gene_constraint_ui("constraint"))
  ),
  .card_when_shown("domains", protein_domains_ui("domains")),
  .card_when_shown("structure", protein_structure_ui("structure")),
  .card_when_shown("ensembl", ensembl_ui("ensembl")),
  .card_when_shown("gtex", gtex_expression_ui("gtex")),
  layout_columns(
    col_widths = c(6, 6),
    .card_when_shown("string_ppi", string_ppi_ui("string_ppi")),
    .card_when_shown("opentargets", opentargets_ui("opentargets"))
  ),
  .card_when_shown("links", external_links_ui("links"))
)

# Right column: bring-your-own-key chat, grounded in the current search via a
# context tool (see server.R). It is pinned in view while the left column
# scrolls (see the .vr-chat-col rule in www/css/app.css); the chat height is
# viewport-relative so it fills the pinned column.
chat_panel <- byok_chat_ui(
  "chat",
  title = "Ask the assistant",
  height = "calc(100vh - 9rem)",
  # Widen the Model & key drawer so the provider/model controls have room.
  sidebar_width = 470,
  greeting = paste(
    "Hi! Open **Model & key** (the gear button), choose a model, and click",
    "**Connect** — a key set in the environment is used automatically; otherwise",
    "paste your own. Then ask me about the gene or variant you're reviewing."
  ),
  # Clickable example-prompt chips below the input (fill-to-edit on click).
  suggestions = VR_CHAT_SUGGESTIONS
)
