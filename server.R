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

  # The dbSNP rsID drives the gnomAD and ClinVar lookups: use the input
  # directly when it is an rsID, otherwise the one MyVariant resolved.
  variant_rsid <- reactive({
    query <- search()
    if (is.null(query) || is_blank(query$variant)) {
      return(NULL)
    }
    variant <- trimws(query$variant)
    if (grepl("^rs[0-9]+$", variant, ignore.case = TRUE)) {
      return(tolower(variant))
    }
    annotation <- variant_annotation()
    if (isTRUE(annotation$ok) && !is_blank(annotation$rsid)) {
      return(annotation$rsid)
    }
    NULL
  })

  # Warn when the entered gene and variant disagree: if a queryable variant
  # (rsID/HGVS) resolves to a different gene than the one typed, the cards would
  # mix gene-level and variant-level results for two different genes. Non-
  # blocking — the cards still render — but it flags the likely mistake.
  output$search_notice <- renderUI({
    query <- search()
    if (is.null(query) || is_blank(query$variant)) {
      return(NULL)
    }
    annotation <- variant_annotation()
    if (!isTRUE(annotation$ok) || is_blank(annotation$gene)) {
      return(NULL)
    }
    typed <- toupper(trimws(query$gene))
    found <- toupper(trimws(annotation$gene))
    if (!nzchar(found) || typed == found) {
      return(NULL)
    }
    vr_error(sprintf(
      paste(
        "The variant %s is in gene %s, not %s.",
        "Gene-level cards show %s; variant-level cards show %s.",
        "Enter a variant in %s, or search %s on its own."
      ),
      query$variant,
      annotation$gene,
      query$gene,
      query$gene,
      annotation$gene,
      query$gene,
      annotation$gene
    ))
  })

  # Each result module returns its data reactive so the assistant can read what
  # each card shows (gene_summary/variant_summary just render the shared
  # resolved/annotation reactives, so those are reused directly).
  gene_summary_server("gene_summary", resolved)
  variant_summary_server("variant_summary", variant_annotation)
  predictions_data <- predictions_server("predictions", search)
  protein_data <- protein_summary_server(
    "protein_summary",
    resolved,
    search,
    variant_annotation
  )
  domains_data <- protein_domains_server(
    "domains",
    resolved,
    search,
    variant_annotation
  )
  structure_data <- protein_structure_server(
    "structure",
    resolved,
    search,
    variant_annotation
  )
  clinvar_data <- clinvar_server("clinvar", variant_rsid)
  gnomad_data <- gnomad_server("gnomad", variant_rsid)
  constraint_data <- gene_constraint_server("constraint", resolved)
  ensembl_data <- ensembl_server("ensembl", variant_rsid)
  gtex_data <- gtex_expression_server("gtex", resolved)
  string_data <- string_ppi_server("string_ppi", resolved)
  opentargets_data <- opentargets_server("opentargets", resolved)
  external_links_server("links", resolved)
  # Visualization cards. The ancestry card reuses the shared gnomAD result (no
  # extra fetch); the gene model reuses the VEP result for the variant position.
  landscape_data <- variant_landscape_server(
    "landscape",
    resolved,
    search,
    variant_annotation
  )
  conservation_data <- conservation_server("conservation", variant_rsid)
  genemodel_data <- gene_model_server("genemodel", resolved, gnomad_data)
  gnomad_ancestry_server("gnomad_ancestry", gnomad_data)

  # --- AI assistant ---------------------------------------------------------
  # Mirror each card's current data into a plain (non-reactive) store so the
  # assistant's tools can read it during async streaming, which runs outside any
  # reactive context. One observer per card keeps its slot in step; reading the
  # reactives here also means the cards populate on search whether or not they
  # are currently visible.
  dash <- new.env(parent = emptyenv())
  dash$selection <- "Nothing is loaded yet."
  observe({
    query <- search()
    dash$selection <- if (is.null(query)) {
      "Nothing is loaded yet."
    } else {
      paste0(
        "Gene: ",
        query$gene,
        if (!is.null(query$variant)) {
          paste0("; Variant: ", query$variant)
        } else {
          "; no variant"
        }
      )
    }
  })
  observe(dash$gene <- resolved())
  observe(dash$variant <- variant_annotation())
  observe(dash$predictions <- predictions_data())
  observe(dash$protein <- protein_data())
  observe(dash$domains <- domains_data())
  observe(dash$structure <- structure_data())
  observe(dash$clinvar <- clinvar_data())
  observe(dash$gnomad <- gnomad_data())
  observe(dash$constraint <- constraint_data())
  observe(dash$consequences <- ensembl_data())
  observe(dash$expression <- gtex_data())
  observe(dash$interactions <- string_data())
  observe(dash$diseases <- opentargets_data())
  observe(dash$landscape <- landscape_data())
  observe(dash$conservation <- conservation_data())
  observe(dash$genemodel <- genemodel_data())

  # Load a gene/variant into the dashboard on the assistant's behalf: update the
  # (namespaced) search inputs and set the shared query, which drives the whole
  # pipeline. Runs the app's own lookups — the assistant never fetches directly.
  load_selection <- function(gene, variant = NULL) {
    gene <- trimws(as.character(gene %||% ""))
    if (!nzchar(gene)) {
      return("No gene provided; nothing was loaded.")
    }
    variant <- trimws(as.character(variant %||% ""))
    # Same format-level gate the search box uses, so the assistant can't fire
    # lookups on a malformed identifier either.
    check <- vr_validate_query(gene, if (nzchar(variant)) variant else NULL)
    if (!isTRUE(check$ok)) {
      return(paste0(
        "Nothing was loaded — invalid input: ",
        paste(check$errors, collapse = " ")
      ))
    }
    updateTextInput(session, "search-gene", value = gene)
    # The variant control is a selectize; set it via choices + selected so an
    # assistant-supplied rsID/HGVS shows even if it isn't a suggested option.
    updateSelectizeInput(
      session,
      "search-variant",
      choices = if (nzchar(variant)) {
        stats::setNames(variant, variant)
      } else {
        character()
      },
      selected = variant,
      server = FALSE
    )
    search(list(gene = gene, variant = if (nzchar(variant)) variant else NULL))
    paste0(
      "Loading ",
      gene,
      if (nzchar(variant)) paste0(" / ", variant) else "",
      " into the dashboard. The cards are refreshing; read them (read_card) to",
      " see the results."
    )
  }

  # Demo button (navbar): a guided walkthrough. Fill the search inputs with a
  # worked example, "click" Review a beat later so the fill is visible before the
  # cards load, then start the cicerone tour that steps through each card and
  # ends on the assistant. Every card must be shown for its tour anchor to
  # exist, so re-tick them all first. Falls back to a modal when cicerone is
  # absent. The staged timing runs off later::later (the reactiveVal set and the
  # tour start need no reactive context; they just push onto the event loop).
  demo_guide <- vr_demo_tour()
  observeEvent(input$demo, {
    ex <- .gene_search_example
    updateCheckboxGroupInput(
      session,
      "visible_cards",
      selected = names(.dashboard_cards)
    )
    updateTextInput(session, "search-gene", value = ex$gene)
    updateSelectizeInput(
      session,
      "search-variant",
      choices = stats::setNames(ex$variant, ex$variant),
      selected = ex$variant,
      server = FALSE
    )
    # Fill the search inputs but do NOT submit: the walkthrough asks the user to
    # click Review themselves. Its first step highlights the search box, and
    # driver.js keeps the highlighted element interactive, so Review is clickable
    # from within the tour; the cards then load on the user's own click and
    # populate as they step through. (search() is deliberately not called here.)
    if (is.null(demo_guide)) {
      showModal(vr_demo_modal(ex$label))
    } else {
      demo_guide$init(session)$start(session = session)
    }
  })

  # Assistant tools, scoped to this app: read the loaded selection, read any
  # card's data, and load a gene/variant. Guarded on ellmer being installed;
  # the chat module itself degrades gracefully when it is not.
  chat_tools <- list()
  if (requireNamespace("ellmer", quietly = TRUE)) {
    chat_tools <- list(
      ellmer::tool(
        function() dash$selection,
        paste(
          "Report the gene and variant currently loaded in the dashboard —",
          "what the user is reviewing."
        ),
        name = "get_current_selection"
      ),
      ellmer::tool(
        function(card) vr_chat_card_text(card, dash[[card]]),
        paste(
          "Read the information currently shown in one dashboard card for the",
          "loaded gene/variant. Returns a text summary of that card's data."
        ),
        arguments = list(
          card = ellmer::type_enum(
            paste0(
              "Which card to read. One of: ",
              paste(names(VR_CHAT_CARDS), collapse = ", "),
              "."
            ),
            values = names(VR_CHAT_CARDS)
          )
        ),
        name = "read_card"
      ),
      ellmer::tool(
        function(gene, variant = NULL) load_selection(gene, variant),
        paste(
          "Load a human gene (and optional variant) into the dashboard, running",
          "the app's own lookups so every card populates. Use this to pull up a",
          "gene/variant for the user. Do not perform your own external",
          "searches — always load through this tool and read the cards."
        ),
        arguments = list(
          gene = ellmer::type_string("Human gene symbol, e.g. TP53 or BRAF."),
          variant = ellmer::type_string(
            paste(
              "Optional variant: an rsID (rs...) or HGVS string. Omit for",
              "gene-only."
            ),
            required = FALSE
          )
        ),
        name = "set_selection"
      )
    )
  }

  byok_chat_server(
    "chat",
    system_prompt = paste(
      "You are a genomics assistant embedded in Variant Reviewer, a gene and",
      "variant interpretation dashboard.",
      "",
      "SCOPE — you help with human gene and variant interpretation:",
      "clinical significance, molecular mechanism, population frequency,",
      "functional impact, gene biology, and explaining the annotations shown in",
      "the app.",
      "",
      "TOOLS — work through the app, never through your own external lookups:",
      "use get_current_selection to see what gene/variant is loaded; use",
      "read_card to read what a specific card shows (gene, variant, protein,",
      "clinvar, gnomad, consequences, expression, interactions, diseases); and",
      "use set_selection to load a gene/variant into the dashboard for the user.",
      "After set_selection, the cards refresh asynchronously — read them again",
      "(read_card) on the next exchange to report results.",
      "",
      "OUT OF SCOPE — politely decline and steer back on topic if asked for",
      "anything unrelated to genomics or variant interpretation. You do not",
      "provide medical diagnosis, treatment, or personal genetic-counseling",
      "advice; for those, direct the user to a qualified clinician or genetic",
      "counselor.",
      "",
      "STYLE — be concise, flag uncertainty, and never fabricate identifiers,",
      "statistics, or citations; say when you don't know."
    ),
    tools = chat_tools,
    # Example prompts: fill-to-edit chips at the input and suggestion cards in
    # the connected greeting (shared constant, see R/chat_tools.R).
    suggestions = VR_CHAT_SUGGESTIONS
  )
}
