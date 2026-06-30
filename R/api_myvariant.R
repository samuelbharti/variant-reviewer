# MyVariant.info client: annotate a variant given an rsID or HGVS string.
# Docs: https://docs.myvariant.info/en/latest/

MYVARIANT_BASE <- "https://myvariant.info/v1"

# MyVariant resolves rsIDs and HGVS strings, but free-text matches a bare
# protein change (e.g. "R175H") to an arbitrary variant. Only query for inputs
# that look like an rsID or an HGVS string (which contains a ":").
myvariant_is_queryable <- function(variant) {
  if (is_blank(variant)) {
    return(FALSE)
  }
  term <- trimws(as.character(variant))
  grepl("^rs[0-9]+$", term, ignore.case = TRUE) ||
    grepl(":", term, fixed = TRUE)
}

# Returns:
#   list(ok = TRUE, id, rsid, gene, hgvsp, cadd_phred, clinvar_significance)
#   list(ok = FALSE, error = "...")
myvariant_annotate <- function(variant) {
  if (is_blank(variant)) {
    return(list(ok = FALSE, error = "No variant supplied."))
  }
  if (!myvariant_is_queryable(variant)) {
    return(list(
      ok = FALSE,
      error = "Enter an rsID (rs...) or HGVS (e.g. chr7:g.140453136A>G) for variant-level annotation."
    ))
  }
  term <- trimws(as.character(variant))

  res <- vr_api_get(
    MYVARIANT_BASE,
    path = "query",
    query = list(
      q = term,
      size = 1,
      fields = paste(
        "dbsnp.rsid",
        "dbnsfp.genename",
        "dbnsfp.hgvsp",
        "cadd.phred",
        "clinvar.rcv.clinical_significance",
        sep = ","
      )
    ),
    source = "MyVariant"
  )
  if (!res$ok) {
    return(list(ok = FALSE, error = res$error))
  }

  hits <- res$data$hits
  if (is.null(hits) || length(hits) == 0) {
    return(list(
      ok = FALSE,
      error = paste0("No annotation found for '", term, "'.")
    ))
  }
  myvariant_parse_hit(hits[[1]], term)
}

# Pure parser: turn a single MyVariant hit into the normalized result.
myvariant_parse_hit <- function(hit, term = NA_character_) {
  list(
    ok = TRUE,
    id = pluck_at(hit, "_id", default = term),
    rsid = mygene_first(pluck_at(hit, "dbsnp", "rsid")),
    gene = mygene_first(pluck_at(hit, "dbnsfp", "genename")),
    hgvsp = mygene_first(pluck_at(hit, "dbnsfp", "hgvsp")),
    cadd_phred = pluck_at(hit, "cadd", "phred", default = NA),
    clinvar_significance = myvariant_clinvar_sig(hit)
  )
}

# clinvar.rcv may be a single object or a list of RCV records; collapse the
# distinct clinical significance values into one readable string.
myvariant_clinvar_sig <- function(hit) {
  rcv <- pluck_at(hit, "clinvar", "rcv")
  if (is.null(rcv)) {
    return(NA_character_)
  }
  sigs <- if (!is.null(rcv$clinical_significance)) {
    rcv$clinical_significance
  } else {
    lapply(rcv, function(r) r$clinical_significance)
  }
  sigs <- unique(unlist(sigs, use.names = FALSE))
  sigs <- sigs[!is.na(sigs) & nzchar(sigs)]
  if (length(sigs) == 0) NA_character_ else paste(sigs, collapse = "; ")
}
