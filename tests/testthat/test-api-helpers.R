# Unit tests for the pure HTTP/utility helpers in R/api_http.R.

test_that("is_blank() detects empty-ish values", {
  expect_true(is_blank(NULL))
  expect_true(is_blank(NA))
  expect_true(is_blank(""))
  expect_true(is_blank("   "))
  expect_true(is_blank(character(0)))
  expect_false(is_blank("TP53"))
  expect_false(is_blank(0))
})

test_that("pluck_at() walks nested lists with a default", {
  x <- list(ensembl = list(gene = "ENSG00000141510"))
  expect_equal(pluck_at(x, "ensembl", "gene"), "ENSG00000141510")
  expect_null(pluck_at(x, "ensembl", "missing"))
  expect_equal(pluck_at(x, "missing", default = "fallback"), "fallback")
})

test_that("%||% coalesces NULL", {
  expect_equal("a" %||% "b", "a")
  expect_equal(NULL %||% "b", "b")
})

test_that("vr_num() formats numbers and handles missing", {
  expect_equal(vr_num(21.234, 1), "21.2")
  expect_equal(vr_num(NA), "—")
  expect_equal(vr_num(NULL), "—")
})
