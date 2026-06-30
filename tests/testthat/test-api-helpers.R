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

test_that("vr_cache_key() is stable and distinguishes inputs", {
  expect_identical(
    vr_cache_key("GET", "a", list(x = 1)),
    vr_cache_key("GET", "a", list(x = 1))
  )
  expect_false(identical(vr_cache_key("GET", "a"), vr_cache_key("GET", "b")))
})

test_that("vr_cached() caches successes but not failures", {
  vr_cache$reset()

  ok_calls <- 0
  ok_fetch <- function() {
    ok_calls <<- ok_calls + 1
    list(ok = TRUE, data = ok_calls)
  }
  key_ok <- vr_cache_key("test", "ok")
  first <- vr_cached(key_ok, ok_fetch)
  second <- vr_cached(key_ok, ok_fetch)
  expect_equal(first$data, 1)
  expect_equal(second$data, 1) # served from cache, not re-fetched
  expect_equal(ok_calls, 1)

  fail_calls <- 0
  fail_fetch <- function() {
    fail_calls <<- fail_calls + 1
    list(ok = FALSE, error = "boom")
  }
  key_fail <- vr_cache_key("test", "fail")
  vr_cached(key_fail, fail_fetch)
  vr_cached(key_fail, fail_fetch)
  expect_equal(fail_calls, 2) # failures are not cached

  vr_cache$reset()
})
