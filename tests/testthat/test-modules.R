# Reactive-logic tests for modules using shiny::testServer() (no browser).

test_that("gene_search_server emits the query on submit", {
  testServer(gene_search_server, {
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
