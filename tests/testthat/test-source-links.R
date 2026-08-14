# Pure URL builders for the per-card "view on source" links (R/source_links.R).

test_that("source URL builders return the expected links", {
  expect_equal(src_ncbi_gene("673"), "https://www.ncbi.nlm.nih.gov/gene/673")
  expect_equal(
    src_dbsnp("rs113488022"),
    "https://www.ncbi.nlm.nih.gov/snp/rs113488022"
  )
  expect_equal(src_gtex("BRAF"), "https://gtexportal.org/home/gene/BRAF")
  expect_match(src_string("BRAF"), "identifiers=BRAF&species=9606")
  expect_equal(
    src_opentargets_gene("ENSG00000157764"),
    "https://platform.opentargets.org/target/ENSG00000157764"
  )
  expect_match(src_gnomad_gene("ENSG00000157764"), "gene/ENSG00000157764")
  expect_match(src_gnomad_variant("7-140753336-A-T"), "variant/7-140753336-A-T")
  expect_equal(
    src_uniprot("P15056"),
    "https://www.uniprot.org/uniprotkb/P15056/entry"
  )
  expect_match(
    src_uniprot("P15056", "family_and_domains"),
    "#family_and_domains$"
  )
  expect_equal(
    src_alphafold("P15056"),
    "https://alphafold.ebi.ac.uk/entry/P15056"
  )
  expect_match(
    src_ensembl_variant("rs113488022"),
    "Variation/Explore\\?v=rs113488022"
  )
  expect_match(src_clinvar_variation("40389"), "clinvar/variation/40389/$")
  expect_equal(
    src_monarch_gene("11998"),
    "https://monarchinitiative.org/HGNC:11998"
  )
  expect_equal(src_monarch_gene("HGNC:11998"), src_monarch_gene("11998"))
  expect_match(
    src_opentargets_drugs("ENSG00000157764"),
    "target/ENSG00000157764/known_drugs$"
  )
  expect_match(
    src_opentargets_pgx("ENSG00000157764"),
    "target/ENSG00000157764/pharmacogenomics$"
  )
  expect_match(
    src_europepmc_search("\"BRAF\""),
    "^https://europepmc.org/search\\?query=%22BRAF%22$"
  )
})

test_that("source URL builders return NULL for blank identifiers", {
  expect_null(src_ncbi_gene(NA))
  expect_null(src_dbsnp(""))
  expect_null(src_uniprot(NULL))
  expect_null(src_alphafold(NA_character_))
})

test_that("vr_source_link renders an anchor, or NULL without an href", {
  expect_null(vr_source_link(NULL))
  html <- as.character(vr_source_link("https://example.org", "Source"))
  expect_match(html, "href=\"https://example.org\"")
  expect_match(html, "Source")
})

test_that("vr_card_header includes a refresh button at ns('refresh')", {
  ns <- NS("mycard")
  html <- as.character(vr_card_header("My card", ns))
  expect_match(html, 'id="mycard-refresh"')
})

test_that("vr_retry_counter's dep/bump wire a reactive to re-run on demand", {
  shiny::testServer(
    function(id) {
      moduleServer(id, function(input, output, session) {
        retry <- vr_retry_counter()
        n_runs <- 0L
        data <- reactive({
          retry$dep()
          n_runs <<- n_runs + 1L
          n_runs
        })
        list(data = data, bump = retry$bump)
      })
    },
    {
      expect_equal(session$returned$data(), 1L)
      expect_equal(session$returned$data(), 1L) # cached: no re-run yet
      session$returned$bump()
      session$flushReact()
      expect_equal(session$returned$data(), 2L)
    }
  )
})

test_that("vr_card_refresh_observer's bump fires on the header's refresh click", {
  bumped <- 0L
  testServer(
    function(id) {
      moduleServer(id, function(input, output, session) {
        vr_card_refresh_observer(input, function() bumped <<- bumped + 1L)
      })
    },
    {
      session$setInputs(refresh = 1)
      expect_equal(bumped, 1L)
    }
  )
})
