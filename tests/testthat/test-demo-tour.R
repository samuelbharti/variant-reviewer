# Demo walkthrough helpers (R/demo_tour.R).

test_that("vr_demo_tour builds a cicerone guide with a step per anchor", {
  skip_if_not_installed("cicerone")
  guide <- vr_demo_tour()
  expect_s3_class(guide, "Cicerone")
  # One step per entry in .vr_demo_steps; the first anchors the search box and
  # the last the chat panel.
  expect_identical(.vr_demo_steps[[1]]$el, "tour_search")
  expect_identical(.vr_demo_steps[[length(.vr_demo_steps)]]$el, "tour_chat")
  # Every card in the dashboard grid has a matching tour step.
  card_anchors <- paste0("tour_", names(.dashboard_cards))
  step_anchors <- vapply(.vr_demo_steps, function(s) s$el, character(1))
  expect_true(all(card_anchors %in% step_anchors))
})

test_that("vr_demo_modal names the example and the assistant steps", {
  html <- as.character(vr_demo_modal("BRAF V600E"))
  expect_match(html, "BRAF V600E", fixed = TRUE)
  # The ampersand is HTML-escaped once rendered.
  expect_match(html, "Model &amp; key", fixed = TRUE)
})
