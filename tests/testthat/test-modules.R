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
        position = 600L,
        stringsAsFactors = FALSE
      ),
      total = 1L
    )
  }
  on.exit(myvariant_gene_variants <<- orig, add = TRUE)

  testServer(gene_search_server, {
    session$setInputs(gene = "BRAF")
    session$elapse(700) # advance past the 600ms debounce
    hint <- paste(as.character(output$variant_hint), collapse = " ")
    expect_match(hint, "pathogenic or likely pathogenic ClinVar submission")
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

test_that("gnomad_allele_frequency() moves an unknown indel left before giving up", {
  # Stub the network: gnomAD knows only the ids in `known`, an rsID lookup
  # answers with `by_rsid`, and the reference is available or not.
  orig_freq <- gnomad_frequency
  orig_align <- ensembl_left_align
  known <- "13-32339421-CA-C"
  by_rsid <- NULL
  aligned <- "13-32339421-CA-C"
  gnomad_frequency <<- function(
    rsid,
    dataset = GNOMAD_DATASET,
    variant_id = NULL
  ) {
    if (is.null(variant_id)) {
      return(by_rsid)
    }
    if (variant_id %in% known) {
      return(list(ok = TRUE, variant_id = variant_id))
    }
    list(ok = FALSE, missing = TRUE, error = "gnomAD has no record.")
  }
  ensembl_left_align <<- function(vcf_id, window = 200L) aligned
  on.exit(
    {
      gnomad_frequency <<- orig_freq
      ensembl_left_align <<- orig_align
    },
    add = TRUE
  )

  # BRCA2 c.5073del as MyVariant writes it: found at its leftmost position.
  res <- gnomad_allele_frequency("rs80359479", "13-32339427-AA-A")
  expect_true(res$ok)
  expect_equal(res$variant_id, "13-32339421-CA-C")

  # A substitution is not moved: its "no record" stands.
  res <- gnomad_allele_frequency("rs1", "7-140753336-A-G")
  expect_true(res$missing)

  # No reference: the rsID answer is kept only for the same change.
  aligned <- NULL
  known <- "7-117559590-ATCT-A"
  by_rsid <- list(ok = TRUE, variant_id = "7-117559590-ATCT-A")
  res <- gnomad_allele_frequency("rs113993960", "7-117559591-TCTT-T")
  expect_equal(res$variant_id, "7-117559590-ATCT-A")

  # No reference and another allele from the rsID: the card says the
  # absence was not checked, not that gnomAD has no record.
  by_rsid <- list(ok = TRUE, variant_id = "7-117559594-T-TCTT")
  res <- gnomad_allele_frequency("rs113993960", "7-117559591-TCTT-T")
  expect_false(res$ok)
  expect_null(res$missing)
  expect_match(res$error, "does not show the variant is missing")
})

test_that("an rsID's alleles are offered in the Variant box", {
  # Spy on the selectize updates the module makes. The module is defined in
  # the global environment, so a function of the same name there is found
  # before shiny's.
  sent <- NULL
  assign(
    "updateSelectizeInput",
    function(session, inputId, choices = NULL, selected = NULL, ...) {
      sent <<- list(choices = choices, selected = selected)
    },
    envir = globalenv()
  )
  on.exit(rm("updateSelectizeInput", envir = globalenv()), add = TRUE)

  requested <- reactiveVal(NULL)
  alleles <- reactiveVal(NULL)
  testServer(
    gene_search_server,
    args = list(requested = requested, allele_choices = alleles),
    {
      # The assistant searches an rsID; the browser still shows the old value.
      session$setInputs(variant = "chr7:g.140753336A>T")
      requested(list(gene = "", variant = "rs80338939", nonce = 1L))
      session$flushReact()
      alleles(c(
        "p.Gly12fs (chr13:g.20189547del)" = "chr13:g.20189547del",
        "chr13:g.20189546_20189547dup" = "chr13:g.20189546_20189547dup"
      ))
      session$flushReact()
      expect_true("chr13:g.20189547del" %in% sent$choices)
      expect_true("p.Gly12fs (chr13:g.20189547del)" %in% names(sent$choices))
      # The searched rsID stays selected, not the stale browser value.
      expect_equal(sent$selected, "rs80338939")

      # The next search has no alleles: the old ones leave the list.
      requested(list(gene = "", variant = "chr12:g.25245350C>T", nonce = 2L))
      session$flushReact()
      alleles(NULL)
      session$flushReact()
      expect_false("chr13:g.20189547del" %in% sent$choices)
      expect_equal(sent$selected, "chr12:g.25245350C>T")
    }
  )
})
