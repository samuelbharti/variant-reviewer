# Format-level validation of the search inputs (R/validate.R), which gates every
# API call. Uses biobouncer's offline "pattern" mode.

test_that("vr_validate_gene accepts well-formed symbols and rejects malformed", {
  expect_true(vr_validate_gene("BRAF")$ok)
  expect_true(vr_validate_gene("tp53")$ok)
  expect_true(vr_validate_gene("  BRAF  ")$ok)

  blank <- vr_validate_gene("")
  expect_false(blank$ok)
  expect_match(blank$error, "Enter a gene symbol")

  bad <- vr_validate_gene("TP 53")
  expect_false(bad$ok)
  expect_match(bad$error, "not a valid gene symbol")
})

test_that("vr_validate_variant checks rsIDs but passes other forms through", {
  # Blank is allowed (variant is optional).
  expect_true(vr_validate_variant("")$ok)
  expect_true(vr_validate_variant(NULL)$ok)

  # Well-formed rsIDs pass (case-insensitively).
  expect_true(vr_validate_variant("rs113488022")$ok)
  expect_true(vr_validate_variant("RS113488022")$ok)

  # Malformed rsIDs are rejected.
  bad <- vr_validate_variant("rs12ab")
  expect_false(bad$ok)
  expect_match(bad$error, "not a valid dbSNP rsID")

  # Non-rsID forms are passed through for the annotation API to resolve.
  expect_true(vr_validate_variant("R175H")$ok)
  expect_true(vr_validate_variant("NM_004333.4:c.1799T>A")$ok)
})

test_that("vr_validate_query collects errors from both fields", {
  ok <- vr_validate_query("BRAF", "rs113488022")
  expect_true(ok$ok)

  ok_gene_only <- vr_validate_query("TP53")
  expect_true(ok_gene_only$ok)

  bad <- vr_validate_query("TP 53", "rs12ab")
  expect_false(bad$ok)
  expect_length(bad$errors, 2)
})

test_that("vr_validate_query allows a variant alone but needs at least one field", {
  expect_true(vr_validate_query(NULL, "rs113488022")$ok)
  expect_true(vr_validate_query("", "rs113488022")$ok)

  both_blank <- vr_validate_query("", "")
  expect_false(both_blank$ok)
  expect_match(both_blank$errors, "gene symbol or a variant")
})
