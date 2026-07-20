# 3D structure card. Shows the AlphaFold-predicted structure for the gene's
# protein (by UniProt accession) with the variant residue highlighted, via the
# r3dmol viewer. The coordinate file is downloaded only when the viewer renders.

protein_structure_ui <- function(id) {
  ns <- NS(id)
  card(
    full_screen = TRUE,
    vr_card_header("3D structure (AlphaFold)", ns),
    card_body(shinycssloaders::withSpinner(
      uiOutput(ns("content")),
      proxy.height = "360px"
    ))
  )
}

# resolved:   reactive() -> mygene_resolve() result (uses the UniProt accession)
# search:     reactive() -> list(gene, variant)
# annotation: reactive() -> myvariant_annotate() result (for the HGVS position)
protein_structure_server <- function(id, resolved, search, annotation) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Lightweight metadata only (model URL + residue) -- no coordinate download,
    # so this is safe to read from the assistant snapshot on every search.
    meta <- reactive({
      res <- resolved()
      if (is.null(res) || !isTRUE(res$ok)) {
        return(NULL)
      }
      if (is_blank(res$uniprot)) {
        return(list(ok = FALSE, error = "No UniProt accession for this gene."))
      }
      out <- alphafold_model(res$uniprot)
      if (isTRUE(out$ok)) {
        out$accession <- res$uniprot
        query <- search()
        out$position <- if (is.null(query) || is_blank(query$variant)) {
          NULL
        } else {
          protein_resolve_position(query$variant, annotation())
        }
      }
      out
    })

    output$viewer <- r3dmol::renderR3dmol({
      info <- meta()
      req(info, isTRUE(info$ok))
      pdb <- alphafold_pdb_text(info$pdb_url)
      req(isTRUE(pdb$ok))
      viewer <- r3dmol::r3dmol(backgroundColor = "#faf8f3") |>
        r3dmol::m_add_model(data = pdb$text, format = "pdb") |>
        r3dmol::m_set_style(
          style = r3dmol::m_style_cartoon(color = "#7a7468")
        ) |>
        r3dmol::m_zoom_to()
      pos <- suppressWarnings(as.integer(info$position %||% NA))
      if (!is.na(pos)) {
        viewer <- viewer |>
          r3dmol::m_add_style(
            sel = r3dmol::m_sel(resi = pos),
            style = r3dmol::m_style_stick(color = "#c07a52", radius = 0.3)
          ) |>
          r3dmol::m_add_style(
            sel = r3dmol::m_sel(resi = pos),
            style = r3dmol::m_style_sphere(color = "#c07a52", scale = 0.5)
          ) |>
          r3dmol::m_zoom_to(sel = r3dmol::m_sel(resi = pos))
      }
      viewer
    })

    output$source <- renderUI({
      info <- meta()
      req(!is.null(info), isTRUE(info$ok))
      vr_source_link(src_alphafold(info$accession), "AlphaFold")
    })

    output$content <- renderUI({
      info <- meta()
      if (is.null(info)) {
        return(vr_empty("Search for a gene to see its predicted 3D structure."))
      }
      if (!isTRUE(info$ok)) {
        return(vr_error(info$error))
      }
      pos <- info$position
      note <- if (!is.null(pos) && !is.na(suppressWarnings(as.integer(pos)))) {
        sprintf(
          "AlphaFold model for %s; variant residue %s highlighted in clay.",
          info$accession,
          pos
        )
      } else {
        sprintf("AlphaFold model for %s.", info$accession)
      }
      tagList(
        r3dmol::r3dmolOutput(ns("viewer"), height = "360px"),
        tags$p(
          class = "text-muted small mb-0 mt-2",
          note,
          " Predicted (computed) structure from AlphaFold DB, not experimental."
        )
      )
    })

    # Returned so the parent can surface this card's metadata to the assistant.
    meta
  })
}
