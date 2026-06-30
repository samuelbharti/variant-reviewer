# Parser tests for the API clients. These exercise the pure parse functions
# against recorded JSON fixtures, so they run offline with no network.

read_fixture <- function(name) {
  jsonlite::fromJSON(
    test_path("fixtures", name),
    simplifyVector = FALSE
  )
}

test_that("mygene_parse_hit() normalizes a MyGene hit", {
  hits <- read_fixture("mygene_tp53.json")$hits
  res <- mygene_parse_hit(hits[[1]], fallback_symbol = "TP53")

  expect_true(res$ok)
  expect_equal(res$symbol, "TP53")
  expect_equal(res$entrez, "7157")
  expect_equal(res$ensembl_gene, "ENSG00000141510")
  expect_equal(res$uniprot, "P04637")
  expect_match(res$summary, "tumor suppressor")
})

test_that("mygene query helpers detect id types and clean input", {
  expect_equal(
    mygene_query_term("ENSG00000141510"),
    "ensembl.gene:ENSG00000141510"
  )
  expect_equal(mygene_query_term("7157"), "entrezgene:7157")
  expect_equal(mygene_query_term("TP53"), "TP53")
  expect_equal(mygene_clean_symbol("  TP53; "), "TP53")
  expect_null(mygene_clean_symbol("   "))
})

test_that("myvariant_is_queryable() accepts rsIDs/HGVS, rejects bare changes", {
  expect_true(myvariant_is_queryable("rs113488022"))
  expect_true(myvariant_is_queryable("chr7:g.140453136A>G"))
  expect_true(myvariant_is_queryable("NM_004333.4:c.1799T>A"))
  expect_false(myvariant_is_queryable("R175H")) # bare protein change
  expect_false(myvariant_is_queryable("175"))
  expect_false(myvariant_is_queryable(""))
})

test_that("myvariant_parse_hit() extracts key annotations", {
  hits <- read_fixture("myvariant_braf.json")$hits
  res <- myvariant_parse_hit(hits[[1]], term = "rs113488022")

  expect_true(res$ok)
  expect_equal(res$rsid, "rs113488022")
  expect_equal(res$gene, "BRAF")
  expect_match(res$hgvsp, "^p\\.")
  expect_true(is.numeric(res$cadd_phred) || is.na(res$cadd_phred))
})

test_that("gtex_parse_rows() builds a tissue/median data.frame", {
  rows <- read_fixture("gtex_tp53.json")$data
  df <- gtex_parse_rows(rows)

  expect_s3_class(df, "data.frame")
  expect_named(df, c("tissue", "median_tpm"))
  expect_true(all(df$median_tpm > 0))
  expect_false(any(grepl("_", df$tissue))) # underscores prettified to spaces
})

test_that("string_parse_rows() builds a sorted partner table", {
  rows <- read_fixture("string_tp53.json")
  df <- string_parse_rows(rows)

  expect_named(
    df,
    c(
      "partner",
      "score",
      "experimental",
      "database",
      "coexpression",
      "textmining"
    )
  )
  expect_equal(df$score, sort(df$score, decreasing = TRUE)) # sorted by score
  expect_true(all(df$score >= 0 & df$score <= 1))
})

test_that("protvar parsers extract function text and variants", {
  fn <- read_fixture("protvar_function_p04637_175.json")
  expect_match(protvar_function_text(fn), "transcription factor")

  pop <- read_fixture("protvar_population_p04637_175.json")
  variants <- protvar_variants_df(pop)
  expect_s3_class(variants, "data.frame")
  expect_named(variants, c("change", "sources"))
  expect_true(nrow(variants) >= 1)
})

test_that("protvar_parse_position() pulls the residue number", {
  expect_equal(protvar_parse_position("R175H"), 175L)
  expect_equal(protvar_parse_position("p.Arg175His"), 175L)
  expect_equal(protvar_parse_position("175"), 175L)
  expect_null(protvar_parse_position("no digits"))
  expect_null(protvar_parse_position(NULL))
})

test_that("opentargets_parse_rows() builds a disease/score data.frame", {
  rows <- read_fixture(
    "opentargets_tp53.json"
  )$data$target$associatedDiseases$rows
  df <- opentargets_parse_rows(rows)

  expect_named(df, c("disease", "disease_id", "score"))
  expect_match(df$disease_id[1], "^MONDO_")
  expect_true(all(df$score >= 0 & df$score <= 1))
  expect_equal(df$score, sort(df$score, decreasing = TRUE)) # API returns sorted
})

test_that("clinvar_parse_record() extracts classification and conditions", {
  record <- read_fixture("clinvar_40389.json")
  res <- clinvar_parse_record(record, uid = "40389")

  expect_true(res$ok)
  expect_equal(res$uid, "40389")
  expect_equal(res$significance, "Pathogenic")
  expect_match(res$review_status, "expert panel")
  expect_match(res$conditions, "RASopathy")
  expect_match(res$accession, "^VCV")
})

test_that("gnomad_freq_part() normalizes a frequency block and handles NULL", {
  part <- gnomad_freq_part(list(af = 1.37e-6, ac = 2, an = 1460618))
  expect_equal(part$ac, 2)
  expect_equal(part$an, 1460618)
  expect_true(part$af > 0)
  expect_null(gnomad_freq_part(NULL))
})

test_that("gnomad_fmt_af() keeps tiny frequencies readable", {
  expect_equal(gnomad_fmt_af(1.37e-6), "1.37e-06")
  expect_equal(gnomad_fmt_af(NA), "—")
})

test_that("external_links_build() only includes links with ids present", {
  full <- external_links_build(list(
    symbol = "TP53",
    ensembl_gene = "ENSG00000141510",
    uniprot = "P04637"
  ))
  expect_true(all(
    c("GeneCards", "Ensembl", "UniProt", "Open Targets", "gnomAD") %in%
      names(full)
  ))
  expect_match(full$UniProt, "P04637")

  partial <- external_links_build(list(
    symbol = "TP53",
    ensembl_gene = NA_character_,
    uniprot = NA_character_
  ))
  expect_true("GeneCards" %in% names(partial))
  expect_false("Ensembl" %in% names(partial)) # no ensembl id -> no link
})
