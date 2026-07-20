# Pure formatters that turn each card's data into text for the assistant's
# read_card tool (R/chat_tools.R). No network, no Shiny.

test_that("formatters handle the not-loaded (NULL) and error states", {
  expect_match(vr_chat_gene(NULL), "No gene loaded")
  expect_match(vr_chat_clinvar(NULL), "needs an rsID")
  expect_match(
    vr_chat_gene(list(ok = FALSE, error = "boom")),
    "No data: boom"
  )
})

test_that("gene and variant summaries include the key fields", {
  gene <- list(
    ok = TRUE,
    symbol = "TP53",
    name = "tumor protein p53",
    type_of_gene = "protein-coding",
    entrez = "7157",
    ensembl_gene = "ENSG00000141510",
    uniprot = "P04637",
    summary = "Acts as a tumor suppressor."
  )
  txt <- vr_chat_gene(gene)
  expect_match(txt, "TP53")
  expect_match(txt, "P04637")
  expect_match(txt, "tumor suppressor")

  variant <- list(
    ok = TRUE,
    id = "chr7:g.140453136A>T",
    rsid = "rs113488022",
    gene = "BRAF",
    hgvsp = "p.Val600Glu",
    cadd_phred = 32,
    clinvar_significance = "Pathogenic"
  )
  expect_match(vr_chat_variant(variant), "rs113488022")
  expect_match(vr_chat_variant(variant), "Pathogenic")
})

test_that("gnomAD summary reports exome/genome allele frequencies", {
  res <- list(
    ok = TRUE,
    variant_id = "7-140753336-A-T",
    dataset = "gnomad_r4",
    exome = list(af = 1.37e-06, ac = 2, an = 1461000),
    genome = NULL
  )
  txt <- vr_chat_gnomad(res)
  expect_match(txt, "Exome AF")
  expect_match(txt, "7-140753336-A-T", fixed = TRUE)
})

test_that("table-backed cards summarize their top rows", {
  expr <- list(
    ok = TRUE,
    data = data.frame(
      tissue = c("Testis", "Brain", "Liver"),
      median_tpm = c(120.5, 10.2, 45.9)
    )
  )
  txt <- vr_chat_expression(expr)
  # Highest first: Testis should lead.
  expect_match(txt, "Highest: Testis")
  expect_match(txt, "3 tissues")

  diseases <- list(
    ok = TRUE,
    count = 42,
    data = data.frame(
      disease = c("melanoma", "colorectal cancer"),
      score = c(0.812, 0.744)
    )
  )
  txt <- vr_chat_diseases(diseases)
  expect_match(txt, "melanoma")
  expect_match(txt, "of 42")
})

test_that("constraint and predictions summaries format for the assistant", {
  con <- list(
    ok = TRUE,
    pli = 0.999,
    loeuf = 0.23,
    oe_lof = 0.15,
    oe_mis = 0.58,
    mis_z = 5.5,
    syn_z = 1.0,
    lof_z = 7.3
  )
  txt <- vr_chat_constraint(con)
  expect_match(txt, "pLI")
  expect_match(txt, "LOEUF")

  pred <- list(
    ok = TRUE,
    predictions = list(
      list(name = "REVEL", score = 0.672, call = "damaging-leaning"),
      list(name = "SIFT", score = NA_real_, call = "deleterious")
    )
  )
  txt <- vr_chat_predictions(pred)
  expect_match(txt, "REVEL")
  expect_match(txt, "SIFT deleterious")

  expect_match(vr_chat_constraint(NULL), "No gene constraint")
  expect_match(vr_chat_predictions(NULL), "No in-silico")
})

test_that("protein domains summary calls out the feature at the residue", {
  df <- data.frame(
    type = c("DOMAIN", "DOMAIN"),
    label = c("Domain", "Domain"),
    description = c("RBD", "Protein kinase"),
    begin = c(155L, 457L),
    end = c(227L, 717L),
    stringsAsFactors = FALSE
  )
  res <- list(ok = TRUE, accession = "P15056", features = df, position = 600)
  txt <- vr_chat_domains(res)
  expect_match(txt, "Residue 600 is in: Protein kinase")
  expect_match(txt, "2 feature")
  expect_match(vr_chat_domains(NULL), "No protein features")
})

test_that("3D structure summary reports the model and highlighted residue", {
  res <- list(ok = TRUE, accession = "P15056", pdb_url = "x", position = 600)
  expect_match(vr_chat_structure(res), "P15056")
  expect_match(vr_chat_structure(res), "residue 600 highlighted")
  no_pos <- list(
    ok = TRUE,
    accession = "P15056",
    pdb_url = "x",
    position = NULL
  )
  expect_match(vr_chat_structure(no_pos), "is shown")
  expect_match(vr_chat_structure(NULL), "No 3D structure")
})

test_that("dispatch routes by card id and flags unknown cards", {
  expect_identical(
    vr_chat_card_text("gene", NULL),
    vr_chat_gene(NULL)
  )
  expect_match(vr_chat_card_text("bogus", NULL), "Unknown card")
})
