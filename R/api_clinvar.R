# ClinVar client via NCBI E-utilities (JSON): clinical significance, review
# status, and associated conditions for a variant.
# Docs: https://www.ncbi.nlm.nih.gov/books/NBK25500/

EUTILS_BASE <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Look up the ClinVar classification for a search term (an rsID works best).
#
# `allele` is the myvariant_annotate() result for the allele being reviewed, or
# NULL. ClinVar keeps one record per allele, so an rsID can return several
# (rs113488022 returns V600G's and V600E's). With an allele, the record is its
# own: the ClinVar variation id MyVariant has for it, else the record whose
# title names its protein change. Without one, a single record is used and
# several are an error, since picking one would be a guess.
# Returns:
#   list(ok = TRUE, uid, accession, title, significance, review_status,
#        last_evaluated, conditions)
#   list(ok = FALSE, error = "...")
clinvar_classification <- function(term, allele = NULL) {
  if (!is_blank(allele$clinvar_id)) {
    return(clinvar_fetch_record(as.character(allele$clinvar_id)))
  }
  if (is_blank(term)) {
    return(list(
      ok = FALSE,
      error = "No variant identifier for ClinVar lookup."
    ))
  }

  search <- vr_api_get(
    EUTILS_BASE,
    path = "esearch.fcgi",
    query = list(db = "clinvar", term = term, retmode = "json"),
    source = "ClinVar"
  )
  if (!search$ok) {
    return(list(ok = FALSE, error = search$error))
  }
  ids <- as.character(unlist(
    pluck_at(search$data, "esearchresult", "idlist"),
    use.names = FALSE
  ))
  if (length(ids) == 0) {
    return(list(
      ok = FALSE,
      error = paste0("No ClinVar record for ", term, ".")
    ))
  }
  if (is.null(allele) && length(ids) == 1) {
    return(clinvar_fetch_record(ids))
  }
  if (is.null(allele)) {
    return(list(
      ok = FALSE,
      error = paste0(
        "ClinVar has ",
        length(ids),
        " records for ",
        term,
        ", one per allele. Pick one allele in the Variant box to see its classification."
      )
    ))
  }

  # One esummary call for every record, then keep the allele's own.
  summary <- clinvar_esummary(ids)
  if (!summary$ok) {
    return(list(ok = FALSE, error = summary$error))
  }
  records <- lapply(ids, function(id) pluck_at(summary$data, "result", id))
  uid <- clinvar_pick_uid(ids, records, allele$hgvsp_all)
  if (is.null(uid)) {
    return(list(
      ok = FALSE,
      error = paste0(
        "No ClinVar record for ",
        if (is_blank(allele$hgvsp)) allele$id else allele$hgvsp,
        " (",
        term,
        ")."
      )
    ))
  }
  clinvar_parse_record(records[[match(uid, ids)]], uid)
}

# esummary for one or more ClinVar uids, in one request.
clinvar_esummary <- function(ids) {
  vr_api_get(
    EUTILS_BASE,
    path = "esummary.fcgi",
    query = list(
      db = "clinvar",
      id = paste(ids, collapse = ","),
      retmode = "json"
    ),
    source = "ClinVar"
  )
}

# Fetch and parse the record for one uid.
clinvar_fetch_record <- function(uid) {
  summary <- clinvar_esummary(uid)
  if (!summary$ok) {
    return(list(ok = FALSE, error = summary$error))
  }
  record <- pluck_at(summary$data, "result", uid)
  if (is.null(record)) {
    return(list(ok = FALSE, error = "ClinVar summary was unavailable."))
  }
  clinvar_parse_record(record, uid)
}

# Pure helper: the uid of the record whose title names one of `changes` (e.g.
# "(p.Val600Glu)" in "NM_004333.6(BRAF):c.1799T>A (p.Val600Glu)"), or NULL.
clinvar_pick_uid <- function(ids, records, changes) {
  if (length(changes) == 0) {
    return(NULL)
  }
  titles <- vapply(
    records,
    function(r) as.character(pluck_at(r, "title", default = "")),
    character(1)
  )
  needles <- paste0("(", changes, ")")
  named <- vapply(
    titles,
    function(t) any(vapply(needles, grepl, logical(1), x = t, fixed = TRUE)),
    logical(1)
  )
  if (!any(named)) NULL else as.character(ids[[which(named)[[1]]]])
}

# Pure parser: an esummary ClinVar record -> normalized classification list.
clinvar_parse_record <- function(record, uid = NA_character_) {
  germline <- pluck_at(record, "germline_classification")
  list(
    ok = TRUE,
    uid = uid,
    accession = pluck_at(record, "accession", default = NA_character_),
    title = pluck_at(record, "title", default = NA_character_),
    significance = pluck_at(germline, "description", default = NA_character_),
    review_status = pluck_at(
      germline,
      "review_status",
      default = NA_character_
    ),
    last_evaluated = pluck_at(
      germline,
      "last_evaluated",
      default = NA_character_
    ),
    conditions = clinvar_conditions(germline)
  )
}

# Collapse the trait set into a readable, comma-separated condition string.
clinvar_conditions <- function(germline) {
  traits <- pluck_at(germline, "trait_set")
  if (is.null(traits) || length(traits) == 0) {
    return(NA_character_)
  }
  names <- vapply(
    traits,
    function(t) {
      as.character(pluck_at(t, "trait_name", default = NA_character_))
    },
    character(1)
  )
  names <- unique(names[!is.na(names) & nzchar(names)])
  if (length(names) == 0) NA_character_ else paste(names, collapse = "; ")
}
