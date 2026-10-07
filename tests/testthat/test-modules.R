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
        id = "chr7:g.140753336A>T",
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
  retried <- 0L
  testServer(
    gene_summary_server,
    args = list(resolved = resolved, retry_resolved = function() {
      retried <<- retried + 1L
    }),
    {
      # renderUI evaluates lazily; force it and confirm it builds HTML.
      expect_no_error(output$content)

      # The header's refresh button (see vr_card_header()) bumps the retry
      # function this card was given, since it has no fetch of its own.
      session$setInputs(refresh = 1)
      expect_equal(retried, 1L)
    }
  )
})

test_that("allele-level cards ask for an allele when an rsID has several", {
  # Nothing is fetched: each card answers from the unpicked allele alone.
  rsid <- reactive("rs113488022")
  allele <- reactive(list(ambiguous = TRUE))
  testServer(clinvar_server, args = list(rsid = rsid, allele = allele), {
    expect_match(session$returned()$error, "Pick one")
  })
  testServer(gnomad_server, args = list(rsid = rsid, allele = allele), {
    expect_match(session$returned$data()$error, "Pick one")
  })
  testServer(ensembl_server, args = list(rsid = rsid, allele = allele), {
    expect_match(session$returned()$error, "Pick one")
  })
})

test_that("a picked allele without a position is never looked up by rsID", {
  # By rsID, gnomAD and VEP can answer for another allele of the rsID (the CTT
  # duplication at rs113993960 got F508del's frequency), so they refuse.
  rsid <- reactive("rs113993960")
  allele <- reactive(list(
    ok = TRUE,
    id = "chr7:g.117559592_117559594dup",
    vcf_id = NA_character_
  ))
  testServer(gnomad_server, args = list(rsid = rsid, allele = allele), {
    expect_match(session$returned$data()$error, "no genomic position")
  })
  testServer(ensembl_server, args = list(rsid = rsid, allele = allele), {
    expect_match(session$returned()$error, "no genomic position")
  })
})

test_that("gnomad_allele_frequency() keeps an rsID answer only for the same change", {
  # Stub the network lookup: the right-shifted id is unknown to gnomAD, and
  # the rsID answers with a variant id of its own choosing.
  orig <- gnomad_frequency
  answer <- "7-117559590-ATCT-A"
  gnomad_frequency <<- function(
    rsid,
    dataset = GNOMAD_DATASET,
    variant_id = NULL
  ) {
    if (!is.null(variant_id)) {
      return(list(ok = FALSE, missing = TRUE, error = "gnomAD has no record."))
    }
    list(ok = TRUE, variant_id = answer)
  }
  on.exit(gnomad_frequency <<- orig, add = TRUE)

  # F508del written at the right end of its repeat: same change, kept.
  res <- gnomad_allele_frequency("rs113993960", "7-117559591-TCTT-T")
  expect_true(res$ok)
  expect_equal(res$variant_id, "7-117559590-ATCT-A")

  # Another allele of the rsID: not kept, the "no record" stands.
  answer <- "7-117559594-T-TCTT"
  res <- gnomad_allele_frequency("rs113993960", "7-117559591-TCTT-T")
  expect_false(res$ok)
  expect_true(res$missing)
})
