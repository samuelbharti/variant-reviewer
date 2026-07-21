# Reactive-logic tests for modules using shiny::testServer() (no browser).

test_that("gene_search_server is NULL before submit, emits the query after", {
  testServer(gene_search_server, {
    # Before any submit the value is NULL (not an error), so result cards can
    # show their placeholder messages.
    expect_null(session$returned())

    session$setInputs(gene = "TP53", variant = "R175H", submit = 1)
    query <- session$returned()
    expect_equal(query$gene, "TP53")
    expect_equal(query$variant, "R175H")
  })
})

test_that("gene_search_server treats a blank gene as no query", {
  testServer(gene_search_server, {
    session$setInputs(gene = "   ", variant = "", submit = 1)
    expect_null(session$returned())
  })
})

test_that("gene_search_server example click fills inputs but does not submit", {
  testServer(gene_search_server, {
    # Clicking the example only populates the inputs; the user still clicks
    # Review, so no query is emitted yet.
    session$setInputs(example = 1)
    expect_null(session$returned())
  })
})

test_that("gene_search_server prefetches variant suggestions for the gene", {
  # Stub the network lookup so the debounced prefetch runs offline.
  orig <- myvariant_gene_variants
  myvariant_gene_variants <<- function(symbol, ...) {
    list(
      ok = TRUE,
      variants = data.frame(
        rsid = "rs1",
        label = "V600E",
        significance = "Pathogenic",
        cadd = 30,
        stringsAsFactors = FALSE
      )
    )
  }
  on.exit(myvariant_gene_variants <<- orig, add = TRUE)

  testServer(gene_search_server, {
    session$setInputs(gene = "BRAF")
    session$elapse(700) # advance past the 600ms debounce
    hint <- paste(as.character(output$variant_hint), collapse = " ")
    expect_match(hint, "known pathogenic")
    expect_match(hint, "BRAF")
  })
})

test_that("an outside request runs through the same submit as a Review click", {
  # The assistant hands the module a request; the module fills its inputs and
  # submits, so the query comes out exactly as if the user had clicked Review.
  requested <- reactiveVal(NULL)
  testServer(gene_search_server, args = list(requested = requested), {
    expect_null(session$returned())

    requested(list(gene = "BRAF", variant = "rs113488022", nonce = 1L))
    session$flushReact()
    query <- session$returned()
    expect_equal(query$gene, "BRAF")
    expect_equal(query$variant, "rs113488022")
  })
})

test_that("an outside request is validated like a typed one", {
  requested <- reactiveVal(NULL)
  testServer(gene_search_server, args = list(requested = requested), {
    # A malformed gene is rejected at the same gate, so no query is emitted and
    # no downstream lookups fire.
    requested(list(gene = "not a gene!", variant = NULL, nonce = 1L))
    session$flushReact()
    expect_null(session$returned())
    expect_false(is.null(validation()))
  })
})

test_that("gene_search_server returns NULL variant when omitted", {
  testServer(gene_search_server, {
    session$setInputs(gene = "BRCA1", variant = "", submit = 1)
    query <- session$returned()
    expect_equal(query$gene, "BRCA1")
    expect_null(query$variant)
  })
})

test_that("gene_summary_server renders an error state without crashing", {
  resolved <- reactive(list(ok = FALSE, error = "No gene found."))
  testServer(gene_summary_server, args = list(resolved = resolved), {
    # renderUI evaluates lazily; force it and confirm it builds HTML.
    expect_no_error(output$content)
  })
})
