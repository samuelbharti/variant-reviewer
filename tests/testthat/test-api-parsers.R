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
  expect_equal(res$hgnc, "11998")
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

test_that("myvariant_parse_hit() describes the one allele it was given", {
  hits <- read_fixture("myvariant_braf_v600e_hg38.json")$hits
  res <- myvariant_parse_hit(hits[[1]], term = "chr7:g.140753336A>T")

  expect_true(res$ok)
  expect_equal(res$id, "chr7:g.140753336A>T")
  expect_equal(res$rsid, "rs113488022")
  expect_equal(res$gene, "BRAF")
  # snpEff's RefSeq protein change, not dbNSFP's first isoform (p.Val640Glu).
  expect_equal(res$hgvsp, "p.Val600Glu")
  expect_true("p.Val600Glu" %in% res$hgvsp_all)
  expect_true("c.1799T>A" %in% res$hgvsc_all)
  # hg38 records have no CADD block of their own; dbNSFP's score is used.
  expect_true(is.numeric(res$cadd_phred) && !is.na(res$cadd_phred))
  expect_equal(res$clinvar_id, "13961")
  expect_equal(res$vcf_id, "7-140753336-A-T")
})

test_that("myvariant_query_term() quotes HGVS so MyVariant matches it", {
  expect_equal(myvariant_query_term("rs113488022"), "rs113488022")
  expect_equal(
    myvariant_query_term(" chr7:g.140753336A>T "),
    "\"chr7:g.140753336A>T\""
  )
  # A quote typed inside is escaped, so it cannot end the phrase early.
  expect_equal(myvariant_query_term("chr7:g.1\"A>T"), "\"chr7:g.1\\\"A>T\"")
})

test_that("myvariant_rsid() falls back to the rsID in the ClinVar block", {
  expect_equal(
    myvariant_rsid(list(clinvar = list(rsid = "rs113993960"))),
    "rs113993960"
  )
  expect_equal(
    myvariant_rsid(list(
      dbsnp = list(rsid = "rs1"),
      clinvar = list(rsid = "rs2")
    )),
    "rs1"
  )
  expect_true(is.na(myvariant_rsid(list())))
})

test_that("myvariant_hgvsc_all() collects cDNA changes in ClinVar's form", {
  hit <- list(
    snpeff = list(ann = list(hgvs_c = "c.1521_1523delCTT")),
    dbnsfp = list(hgvsc = list("c.620T>A", "c.1799T>A"))
  )
  expect_setequal(
    myvariant_hgvsc_all(hit),
    c("c.1521_1523del", "c.620T>A", "c.1799T>A")
  )
  expect_length(myvariant_hgvsc_all(list()), 0)
})

test_that("myvariant_same_vcf_id() matches one indel written two ways", {
  expect_true(myvariant_same_vcf_id("7-117559590-ATCT-A", "7-117559591-TCTT-T"))
  expect_false(myvariant_same_vcf_id("7-140753336-A-T", "7-140753336-A-C"))
  expect_false(myvariant_same_vcf_id("7-140753336-A-T", NA_character_))
  expect_false(myvariant_same_vcf_id("7-140753336-A-T", "rs113488022"))
})

test_that("myvariant_pick_allele() lists the alleles of a multi-allele rsID", {
  hits <- read_fixture("myvariant_rs113488022_hg38.json")$hits
  res <- myvariant_pick_allele(hits, "rs113488022", "No annotation found for")

  expect_false(res$ok)
  expect_true(res$ambiguous)
  expect_equal(res$gene, "BRAF")
  expect_equal(nrow(res$alleles), 3L)
  expect_setequal(
    res$alleles$hgvsp,
    c("p.Val600Ala", "p.Val600Glu", "p.Val600Gly")
  )
  expect_match(res$error, "p.Val600Glu (chr7:g.140753336A>T)", fixed = TRUE)
})

test_that("myvariant_pick_allele() takes a single hit and reports none", {
  one <- myvariant_pick_allele(
    list(list(`_id` = "chr7:g.140753336A>T")),
    "rs1",
    "No annotation found for"
  )
  expect_true(one$ok)
  expect_equal(one$hit$`_id`, "chr7:g.140753336A>T")

  none <- myvariant_pick_allele(list(), "rs1", "No annotation found for")
  expect_false(none$ok)
  expect_equal(none$error, "No annotation found for 'rs1'.")
  expect_null(none$ambiguous)
})

test_that("myvariant_hgvsp() prefers snpEff over dbNSFP's isoform order", {
  dbnsfp <- list(hgvsp = list("p.Val640Glu", "p.Val600Glu"))
  # snpEff's annotation as one object, and as a list of them.
  expect_equal(
    myvariant_hgvsp(list(
      dbnsfp = dbnsfp,
      snpeff = list(ann = list(hgvs_p = "p.Val600Glu"))
    )),
    "p.Val600Glu"
  )
  expect_equal(
    myvariant_hgvsp(list(
      snpeff = list(
        ann = list(list(hgvs_p = "p.Gly12Asp"), list(hgvs_p = "p.Gly12Asp"))
      )
    )),
    "p.Gly12Asp"
  )
  # No snpEff protein change: the first dbNSFP one.
  expect_equal(
    myvariant_hgvsp(list(
      dbnsfp = dbnsfp,
      snpeff = list(ann = list(effect = "intron_variant"))
    )),
    "p.Val640Glu"
  )
  expect_true(is.na(myvariant_hgvsp(list())))
})

test_that("myvariant_vcf_id() needs every VCF field", {
  hit <- list(
    chrom = "7",
    vcf = list(position = "140753336", ref = "A", alt = "T")
  )
  expect_equal(myvariant_vcf_id(hit), "7-140753336-A-T")
  hit$vcf$alt <- NULL
  expect_true(is.na(myvariant_vcf_id(hit)))
})

test_that("myvariant_distinct_alleles() merges one indel written two ways", {
  hit <- function(id, pos, ref, alt, clinvar = NULL) {
    list(
      `_id` = id,
      chrom = "7",
      vcf = list(position = pos, ref = ref, alt = alt),
      clinvar = if (!is.null(clinvar)) list(variant_id = clinvar)
    )
  }
  # CFTR F508del, shifted right and left inside its repeat, plus the CTT
  # duplication at the same rsID, which MyVariant has without VCF fields.
  right <- hit("chr7:g.117559592_117559594del", "117559591", "TCTT", "T")
  left <- hit(
    "chr7:g.117559591_117559593del",
    "117559590",
    "ATCT",
    "A",
    clinvar = 7105
  )
  dup <- list(`_id` = "chr7:g.117559592_117559594dup")

  kept <- myvariant_distinct_alleles(list(right, left, dup))
  ids <- vapply(kept, function(h) h$`_id`, character(1))
  # The two deletions are one allele, and the one ClinVar knows is kept. The
  # record without VCF fields cannot be compared, so it stays.
  expect_equal(
    ids,
    c("chr7:g.117559591_117559593del", "chr7:g.117559592_117559594dup")
  )

  # Different bases at one position are different alleles.
  snv <- function(alt) list(chrom = "7", pos = 140753336L, ref = "A", alt = alt)
  expect_false(.mv_same_change(snv("T"), snv("C")))
  expect_true(.mv_same_change(snv("T"), snv("T")))
})

test_that("vr_variant_allele() passes one allele and flags an unpicked rsID", {
  expect_equal(vr_variant_allele(list(ok = TRUE, id = "x"))$id, "x")
  expect_true(
    vr_variant_allele(list(ok = FALSE, ambiguous = TRUE, error = "e"))$ambiguous
  )
  expect_null(vr_variant_allele(list(ok = FALSE, error = "down")))
  expect_null(vr_variant_allele(NULL))
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
  # The inline UniProt citations are stripped on the way out.
  expect_no_match(protvar_function_text(fn), "PubMed:")

  pop <- read_fixture("protvar_population_p04637_175.json")
  variants <- protvar_variants_df(pop)
  expect_s3_class(variants, "data.frame")
  expect_named(variants, c("change", "sources"))
  expect_true(nrow(variants) >= 1)
})

test_that("protvar_strip_citations() drops evidence, keeps the prose", {
  # The real TP53 shape: a citation group closing the sentence.
  expect_identical(
    protvar_strip_citations(
      "Induces cell cycle arrest (PubMed:11025664, PubMed:12524540)."
    ),
    "Induces cell cycle arrest."
  )
  # Mixed parenthetical: the note survives, the citation goes.
  expect_identical(
    protvar_strip_citations("Binds DNA (By similarity, PubMed:9840937)."),
    "Binds DNA (By similarity)."
  )
  # Several groups across a longer passage.
  stripped <- protvar_strip_citations(paste(
    "Acts as a tumor suppressor (PubMed:11025664, PubMed:12524540).",
    "Regulates the circadian clock (PubMed:24051492)"
  ))
  expect_no_match(stripped, "PubMed")
  expect_match(stripped, "tumor suppressor\\.")
  expect_match(stripped, "circadian clock$")

  # Nothing to strip, and blank input, both pass through untouched.
  expect_identical(protvar_strip_citations("Plain text."), "Plain text.")
  expect_true(is_blank(protvar_strip_citations(NA_character_)))
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

test_that("opentargets_parse_drugs() builds a drug/phase/disease data.frame", {
  rows <- read_fixture(
    "opentargets_drugs_braf.json"
  )$data$target$drugAndClinicalCandidates$rows
  df <- opentargets_parse_drugs(rows)

  expect_named(df, c("drug", "drug_id", "drug_type", "max_phase", "disease"))
  expect_true(all(nzchar(df$drug)))
  # Clinical stage is prettified from the API's SCREAMING_SNAKE form.
  expect_true(any(grepl("^Phase ", df$max_phase)))
  expect_false(any(grepl("_", df$max_phase))) # no PHASE_2 left
})

test_that("opentargets_pretty_phase() humanizes the stage enum", {
  expect_equal(opentargets_pretty_phase("PHASE_2"), "Phase 2")
  expect_equal(opentargets_pretty_phase("PRE_CLINICAL"), "Pre clinical")
  expect_true(is.na(opentargets_pretty_phase("")))
})

test_that("opentargets_parse_pgx() builds a variant/drug/effect data.frame", {
  rows <- read_fixture(
    "opentargets_pgx_cyp2c19.json"
  )$data$target$pharmacogenomics
  df <- opentargets_parse_pgx(rows)

  expect_named(df, c("rsid", "drug", "phenotype", "genotype", "evidence"))
  expect_equal(nrow(df), length(rows))
  # At least one row names a drug and carries an effect description.
  expect_true(any(!is.na(df$drug)))
  expect_true(any(nzchar(df$phenotype)))
})

test_that("europepmc_parse_results() builds a citation data.frame", {
  results <- read_fixture("europepmc_braf_v600e.json")$resultList$result
  df <- europepmc_parse_results(results)

  expect_named(
    df,
    c("title", "authors", "journal", "year", "id", "source", "doi", "cited_by")
  )
  expect_true(all(nzchar(df$title)))
  expect_true(is.integer(df$cited_by))
  # Europe PMC escapes inline markup in titles (e.g. "&lt;i&gt;"); it should
  # come back decoded so it renders as real tags, not literal "<i>" text.
  expect_true(any(grepl("<i>BRAF V600E</i>", df$title, fixed = TRUE)))
  expect_false(any(grepl("&lt;", df$title, fixed = TRUE)))
})

test_that("europepmc_decode_title() decodes HTML entities", {
  expect_equal(europepmc_decode_title("&lt;i&gt;BRAF&lt;/i&gt;"), "<i>BRAF</i>")
  expect_equal(europepmc_decode_title("A &amp; B"), "A & B")
  expect_equal(europepmc_decode_title("5' &amp; 3&#39;"), "5' & 3'")
  expect_true(is.na(europepmc_decode_title(NA_character_)))
})

test_that("europepmc_query() quotes the gene, ANDs a refinement, and sorts by date", {
  expect_equal(europepmc_query("BRAF"), "\"BRAF\" sort_date:y")
  expect_equal(
    europepmc_query("BRAF", "rs113488022"),
    "\"BRAF\" AND \"rs113488022\" sort_date:y"
  )
})

test_that("monarch_parse_items() builds an HPO id/phenotype data.frame", {
  items <- read_fixture("monarch_phenotypes_tp53.json")$items
  df <- monarch_parse_items(items)

  expect_s3_class(df, "data.frame")
  expect_named(df, c("hpo_id", "phenotype"))
  expect_true(all(grepl("^HP:", df$hpo_id)))
  expect_true(all(nzchar(df$phenotype)))
})

test_that("monarch_hgnc_id() normalizes to the HGNC CURIE", {
  expect_equal(monarch_hgnc_id("11998"), "HGNC:11998")
  expect_equal(monarch_hgnc_id("hgnc:11998"), "HGNC:11998")
  expect_equal(monarch_hgnc_id(" HGNC:11998 "), "HGNC:11998")
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

test_that("clinvar_pick_uid() keeps the record for the allele's cDNA change", {
  result <- read_fixture("clinvar_rs113488022_esummary.json")$result
  ids <- c("40389", "13961")
  records <- lapply(ids, function(id) result[[id]])

  expect_equal(clinvar_pick_uid(ids, records, cdna = "c.1799T>A"), "13961")
  expect_equal(
    clinvar_pick_uid(ids, records, cdna = c("c.620T>G", "c.1799T>G")),
    "40389"
  )
  # V600A (c.1799T>C) has no ClinVar record of its own. A known cDNA change
  # is the only thing matched, so a protein change cannot pull in another
  # allele's record.
  expect_null(
    clinvar_pick_uid(ids, records, cdna = "c.1799T>C", protein = "p.Val600Glu")
  )
  # With no cDNA change known, the protein change is used.
  expect_equal(clinvar_pick_uid(ids, records, protein = "p.Val600Glu"), "13961")
  expect_null(clinvar_pick_uid(ids, records))
})

test_that("gnomad_freq_part() normalizes a frequency block and handles NULL", {
  part <- gnomad_freq_part(list(af = 1.37e-6, ac = 2, an = 1460618))
  expect_equal(part$ac, 2)
  expect_equal(part$an, 1460618)
  expect_true(part$af > 0)
  expect_null(gnomad_freq_part(NULL))
})

test_that("gnomad_error_kind() tells a missing from an unresolved variant", {
  expect_equal(
    gnomad_error_kind(list(list(message = "Variant not found"))),
    "missing"
  )
  expect_equal(
    gnomad_error_kind(list(list(
      message = "Multiple variants found, query using variant ID to select one."
    ))),
    "multiple"
  )
  expect_equal(gnomad_error_kind(list(list(message = "Syntax error"))), "other")
})

test_that("gnomad_fmt_af() keeps tiny frequencies readable", {
  expect_equal(gnomad_fmt_af(1.37e-6), "1.37e-06")
  expect_equal(gnomad_fmt_af(NA), "—")
})

test_that("ensembl_parse_vep() extracts consequence summary and table", {
  record <- read_fixture("ensembl_vep_rs113488022.json")
  res <- ensembl_parse_vep(record)

  expect_true(res$ok)
  expect_equal(res$most_severe, "missense_variant")
  expect_s3_class(res$data, "data.frame")
  expect_named(
    res$data,
    c("gene", "transcript", "consequence", "impact", "sift", "polyphen")
  )
  expect_true(all(res$data$gene == "BRAF"))
})

test_that("ensembl_vep_region() turns a VCF allele into VEP's region form", {
  expect_equal(
    ensembl_vep_region("7-140753336-A-T"),
    "7:140753336-140753336:1/T"
  )
  # Deletion and insertion: VCF's shared leading base is dropped.
  expect_equal(
    ensembl_vep_region("7-117559590-ATCT-A"),
    "7:117559591-117559593:1/-"
  )
  expect_equal(
    ensembl_vep_region("7-117559594-T-TCTT"),
    "7:117559595-117559594:1/CTT"
  )
  expect_null(ensembl_vep_region(NA_character_))
  expect_null(ensembl_vep_region("rs113488022"))
})

test_that("ensembl_consequences_df() keeps only protein-coding rows", {
  coding <- list(
    list(
      biotype = "protein_coding",
      gene_symbol = "BRAF",
      transcript_id = "T1",
      consequence_terms = list("missense_variant"),
      impact = "MODERATE"
    ),
    list(
      biotype = "retained_intron",
      gene_symbol = "BRAF",
      transcript_id = "T2",
      consequence_terms = list("intron_variant"),
      impact = "MODIFIER"
    )
  )
  df <- ensembl_consequences_df(coding)
  expect_equal(nrow(df), 1)
  expect_equal(df$transcript, "T1")

  expect_null(ensembl_consequences_df(NULL))
  expect_null(ensembl_consequences_df(list(list(biotype = "lncRNA"))))
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

test_that("gnomad_parse_constraint() extracts constraint metrics", {
  data <- read_fixture("gnomad_constraint_braf.json")
  res <- gnomad_parse_constraint(data, "BRAF")
  expect_true(res$ok)
  expect_gt(res$pli, 0.99)
  expect_equal(round(res$loeuf, 2), 0.23)
  expect_gt(res$mis_z, 5)
})

test_that("gnomad_parse_constraint() reports missing data", {
  res <- gnomad_parse_constraint(list(data = list(gene = NULL)), "XYZ")
  expect_false(res$ok)
  expect_match(res$error, "no constraint")
})

test_that("myvariant_parse_predictions() summarizes in-silico scores", {
  hit <- read_fixture("myvariant_predictions_braf.json")$hits[[1]]
  res <- myvariant_parse_predictions(hit)
  expect_true(res$ok)
  by_name <- function(nm) {
    Filter(function(p) p$name == nm, res$predictions)[[1]]
  }
  expect_equal(by_name("REVEL")$score, 0.672)
  expect_match(by_name("REVEL")$call, "damaging")
  # AlphaMissense score is a per-transcript array; the parser keeps the max.
  expect_equal(by_name("AlphaMissense")$score, 0.9878)
  expect_match(by_name("AlphaMissense")$call, "pathogenic")
})

test_that("myvariant_predictions() rejects non-queryable input", {
  res <- myvariant_predictions("R175H")
  expect_false(res$ok)
  expect_match(res$error, "rsID")
})

test_that("myvariant_parse_gene_variants() builds a ranked variant table", {
  hits <- read_fixture("myvariant_gene_variants_braf.json")$hits
  res <- myvariant_parse_gene_variants(hits)
  expect_true(res$ok)
  v <- res$variants
  expect_true(all(
    c("id", "rsid", "label", "significance", "cadd") %in% names(v)
  ))
  # One row per allele; rsIDs are lower-cased.
  expect_equal(anyDuplicated(v$id), 0L)
  expect_true(all(grepl("^rs[0-9]+$", v$rsid)))
  # Amino-acid labels take the one-letter ref+pos+alt form (e.g. L485S).
  expect_true(any(grepl("^[A-Z][0-9]+[A-Z*]$", v$label)))
  # Every suggestion is pathogenic/likely-pathogenic, most severe ranked first.
  expect_true(all(v$significance %in% c("Pathogenic", "Likely pathogenic")))
  expect_false(is.unsorted(match(
    v$significance,
    c("Pathogenic", "Likely pathogenic")
  )))
})

test_that("myvariant_parse_gene_variants() keeps each allele of an rsID", {
  allele <- function(id, alt, sig, rsid = "rs1") {
    list(
      `_id` = id,
      dbsnp = list(rsid = rsid),
      dbnsfp = list(aa = list(ref = "V", alt = alt, pos = list(600))),
      clinvar = list(rcv = list(clinical_significance = sig))
    )
  }
  hits <- list(
    allele("chr7:g.140753336A>T", "E", "Pathogenic"),
    # Same rsID, another allele: its own row.
    allele("chr7:g.140753336A>C", "G", "Likely pathogenic", rsid = "RS1"),
    # The same allele twice (kept once), and a hit with no rsID (dropped).
    allele("chr7:g.140753336A>T", "E", "Pathogenic"),
    list(
      `_id` = "chr1:g.1A>T",
      dbnsfp = list(aa = list(ref = "A", alt = "T", pos = list(1)))
    )
  )
  res <- myvariant_parse_gene_variants(hits)
  expect_true(res$ok)
  expect_equal(nrow(res$variants), 2L)
  expect_equal(
    res$variants$id,
    c("chr7:g.140753336A>T", "chr7:g.140753336A>C")
  )
  expect_equal(res$variants$rsid, c("rs1", "rs1"))
  expect_equal(res$variants$label, c("V600E", "V600G"))
})

test_that("myvariant_variant_choices() maps display labels to allele ids", {
  parsed <- myvariant_parse_gene_variants(
    read_fixture("myvariant_gene_variants_braf.json")$hits
  )
  choices <- myvariant_variant_choices(parsed, max_n = 3)
  expect_length(choices, 3)
  expect_true(all(grepl("^chr[0-9XY]+:g\\.", unname(choices))))
  expect_match(names(choices)[1], "\\(")
  # No suggestions for a failed parse.
  expect_length(myvariant_variant_choices(list(ok = FALSE)), 0)
})

test_that("proteins_parse_features() tidies features and finds those at a residue", {
  data <- read_fixture("proteins_features_p15056.json")
  res <- proteins_parse_features(data, "P15056")
  expect_true(res$ok)
  expect_gt(nrow(res$features), 0)
  expect_true(all(
    c("type", "label", "description", "begin", "end") %in% names(res$features)
  ))
  # BRAF V600 sits inside the protein kinase domain (457-717).
  at <- proteins_features_at(res$features, 600)
  expect_true(any(at$description == "Protein kinase"))
  # ...and a position past the last annotated feature matches nothing.
  expect_equal(nrow(proteins_features_at(res$features, 5000)), 0)
})

test_that("proteins_features() rejects a missing accession", {
  res <- proteins_features("")
  expect_false(res$ok)
  expect_match(res$error, "accession")
})

test_that("alphafold_parse_model() extracts the model URL", {
  data <- read_fixture("alphafold_p15056.json")
  res <- alphafold_parse_model(data, "P15056")
  expect_true(res$ok)
  expect_match(res$pdb_url, "^https://.*AF-P15056.*\\.pdb$")
})

test_that("alphafold_parse_model() reports a missing model", {
  res <- alphafold_parse_model(list(), "XYZ")
  expect_false(res$ok)
  expect_match(res$error, "No AlphaFold model")
})

test_that("vr_http_error_message hides technical detail behind plain language", {
  # Transport errors (timeout/DNS/connection) never leak curl internals.
  timeout <- simpleError(
    "Failed to perform HTTP request. Timeout was reached [rest.ensembl.org]"
  )
  msg <- vr_http_error_message("Ensembl VEP", condition = timeout)
  expect_match(msg, "took too long", fixed = TRUE)
  expect_no_match(msg, "curl|Timeout was reached|http request")

  dns <- simpleError("Could not resolve host: mygene.info")
  expect_match(
    vr_http_error_message("MyGene", condition = dns),
    "check your internet connection"
  )

  # HTTP statuses map to their own short messages.
  expect_match(
    vr_http_error_message("ClinVar", status = 404L),
    "No ClinVar data"
  )
  expect_match(vr_http_error_message("gnomAD", status = 429L), "busy right now")
  expect_match(
    vr_http_error_message("STRING", status = 503L),
    "temporarily unavailable"
  )

  # An unclassified transport error still degrades to a safe generic message.
  expect_match(
    vr_http_error_message("GTEx", condition = simpleError("weird boom")),
    "temporarily unavailable"
  )
})

test_that("gnomad_parse_populations sums exome+genome and drops sex splits", {
  exome <- list(
    list(id = "nfe", ac = 10, an = 1000),
    list(id = "afr", ac = 1, an = 500),
    list(id = "nfe_XX", ac = 5, an = 500), # sex split -> dropped
    list(id = "XY", ac = 9, an = 900) # overall sex group -> dropped
  )
  genome <- list(
    list(id = "nfe", ac = 2, an = 200),
    list(id = "eas", ac = 0, an = 300)
  )
  df <- gnomad_parse_populations(exome, genome)
  # Only real ancestry groups survive (no XX/XY, no …_XX/…_XY, no sub-pops).
  expect_setequal(df$pop, c("nfe", "afr", "eas"))
  expect_equal(df$ac[df$pop == "nfe"], 12) # summed across sample sets
  expect_equal(df$an[df$pop == "nfe"], 1200)
  expect_equal(df$label[df$pop == "nfe"], "European (non-Finnish)")
  expect_equal(df$af[df$pop == "afr"], 1 / 500)
  expect_equal(df$af, sort(df$af, decreasing = TRUE)) # sorted by frequency
  expect_null(gnomad_parse_populations(list(), list()))
})

test_that("gnomad_sig_category and gnomad_hgvsp_residue classify inputs", {
  expect_equal(gnomad_sig_category("Pathogenic"), "Pathogenic / likely")
  expect_equal(gnomad_sig_category("Likely pathogenic"), "Pathogenic / likely")
  expect_equal(gnomad_sig_category("Likely benign"), "Benign / likely")
  expect_equal(gnomad_sig_category("Uncertain significance"), "Uncertain")
  expect_equal(
    gnomad_sig_category("Conflicting interpretations of pathogenicity"),
    "Conflicting"
  )
  expect_equal(gnomad_hgvsp_residue("p.Val600Glu"), 600L)
  expect_equal(gnomad_hgvsp_residue("p.V600E"), 600L)
  expect_true(is.na(gnomad_hgvsp_residue(NULL)))
  expect_true(is.na(gnomad_hgvsp_residue("p.=")))
})

test_that("gnomad_parse_clinvar_variants keeps only residue-bearing variants", {
  records <- list(
    list(
      pos = 140753336,
      hgvsp = "p.Val600Glu",
      major_consequence = "missense_variant",
      clinical_significance = "Pathogenic"
    ),
    list(
      pos = 1,
      hgvsp = NULL, # UTR variant, no residue -> dropped
      major_consequence = "3_prime_UTR_variant",
      clinical_significance = "Benign"
    )
  )
  res <- gnomad_parse_clinvar_variants(records, "BRAF")
  expect_true(res$ok)
  expect_equal(nrow(res$variants), 1)
  expect_equal(res$variants$residue, 600L)
  expect_equal(res$variants$category, "Pathogenic / likely")

  expect_false(gnomad_parse_clinvar_variants(list(), "BRAF")$ok)
})

test_that("myvariant_parse_conservation extracts the four metrics", {
  hit <- list(
    dbnsfp = list(
      phylop = list(
        `100way_vertebrate` = list(score = 9.236, rankscore = 0.944)
      ),
      phastcons = list(
        `100way_vertebrate` = list(score = 1.0, rankscore = 0.716)
      ),
      `gerp++` = list(rs = 5.65, rs_rankscore = 0.868),
      siphy_29way = list(logodds_score = 15.93, logodds_rankscore = 0.794)
    )
  )
  res <- myvariant_parse_conservation(hit)
  expect_true(res$ok)
  expect_equal(nrow(res$metrics), 4)
  expect_equal(res$metrics$score[res$metrics$metric == "GERP++ RS"], 5.65)
  expect_true(all(res$metrics$rankscore >= 0 & res$metrics$rankscore <= 1))

  expect_false(myvariant_parse_conservation(list(dbnsfp = list()))$ok)
})

test_that("ensembl_parse_gene_model picks the canonical transcript's exons", {
  record <- list(
    seq_region_name = "7",
    start = 100,
    end = 400,
    Transcript = list(
      list(
        id = "ENST_OTHER",
        is_canonical = 0,
        strand = -1,
        Exon = list(
          list(start = 100, end = 150)
        )
      ),
      list(
        id = "ENST_CANON",
        is_canonical = 1,
        strand = -1,
        Exon = list(
          list(start = 300, end = 400),
          list(start = 100, end = 200)
        )
      )
    )
  )
  res <- ensembl_parse_gene_model(record)
  expect_true(res$ok)
  expect_equal(res$transcript, "ENST_CANON")
  expect_equal(nrow(res$exons), 2)
  expect_equal(res$exons$start, c(100, 300)) # rows sorted by genomic start
  # Minus strand: numbered 5'->3', so the highest-coordinate exon is exon 1.
  expect_equal(res$exons$number, c(2L, 1L))
  expect_equal(res$strand, -1)

  # Plus strand numbers in genomic order instead.
  plus <- ensembl_parse_gene_model(list(
    seq_region_name = "1",
    Transcript = list(list(
      id = "T",
      is_canonical = 1,
      strand = 1,
      Exon = list(
        list(start = 100, end = 200),
        list(start = 300, end = 400)
      )
    ))
  ))
  expect_equal(plus$exons$number, c(1L, 2L))

  expect_false(ensembl_parse_gene_model(list(Transcript = list()))$ok)
})
