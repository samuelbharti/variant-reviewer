# Record the demo example's API responses into demo_cache.rds.
#
# Run this whenever the demo example, an API endpoint, or a query parameter
# changes, so the snapshot the app seeds at startup stays in step:
#
#   Rscript dev/record_demo_cache.R
#
# It needs a live network. It drives the real app server (via shiny::testServer)
# for BRAF V600E with every card visible, so the requests it makes -- and
# therefore the cache keys it records -- are exactly the ones the running app
# makes. Nothing here reconstructs calls by hand, so the snapshot cannot drift
# out of step with what the modules actually fetch. The result is written to
# demo_cache.rds at the project root, which R/demo_cache.R seeds from.

# Project root = parent of this script's dev/ directory. Works from any wd.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg))
} else {
  getwd()
}
app_dir <- normalizePath(file.path(script_dir, ".."))
setwd(app_dir)

out_file <- file.path(app_dir, "demo_cache.rds")
# Passed to the testServer expr through the environment, so it does not depend
# on the working directory testServer runs in.
Sys.setenv(VR_DEMO_SNAPSHOT_OUT = out_file)

library(shiny)

message("Recording demo cache for BRAF V600E from live APIs...")

shiny::testServer(app_dir, {
  # Start from an empty cache so the snapshot holds only the demo's requests.
  vr_cache_clear()

  # Drive the search box exactly as the Demo button does, but with every card
  # visible so the snapshot covers every source (including the default-off gene
  # model), not just the demo's default-on subset.
  session$setInputs(
    visible_cards = names(.dashboard_cards),
    `search-gene` = .gene_search_example$gene,
    `search-variant` = .gene_search_example$variant,
    `search-submit` = 1
  )

  # Fetches are synchronous, so each flush blocks until the sources it touches
  # have returned. A few passes let the reactive graph settle (later sources
  # depend on earlier ones, e.g. the gene model reads the gnomAD result).
  for (i in 1:6) {
    session$flushReact()
  }

  # Dump every cache entry to the snapshot. Keys are the request hashes; values
  # are the normalized result lists. Drop any that pruned out between listing
  # and reading (belt and braces; nothing expires in a short script run).
  keys <- vr_cache$keys()
  entries <- lapply(keys, function(k) vr_cache$get(k, missing = NULL))
  keep <- !vapply(entries, is.null, logical(1))
  snapshot <- stats::setNames(entries[keep], keys[keep])

  saveRDS(snapshot, Sys.getenv("VR_DEMO_SNAPSHOT_OUT"))
  message(sprintf("Recorded %d cache entries.", length(snapshot)))
})

size_kb <- round(file.info(out_file)$size / 1024, 1)
message(sprintf("Wrote %s (%s KB).", out_file, size_kb))
