# Open Targets Platform client: disease associations for a gene (target).
# GraphQL API docs: https://platform-docs.opentargets.org/data-access/graphql-api

OPENTARGETS_URL <- "https://api.platform.opentargets.org/api/v4/graphql"

OPENTARGETS_QUERY <- paste(
  "query($id: String!, $size: Int!) {",
  "  target(ensemblId: $id) {",
  "    approvedSymbol",
  "    associatedDiseases(page: {index: 0, size: $size}) {",
  "      count",
  "      rows { score disease { id name } }",
  "    }",
  "  }",
  "}",
  sep = "\n"
)

# Disease associations for an Ensembl gene id (target), ordered by association
# score (Open Targets returns them already sorted, highest first).
# Returns:
#   list(ok = TRUE, count, data = data.frame(disease, disease_id, score))
#   list(ok = FALSE, error = "...")
opentargets_diseases <- function(ensembl_id, size = 15) {
  if (is_blank(ensembl_id)) {
    return(list(
      ok = FALSE,
      error = "No Ensembl gene ID available for this gene."
    ))
  }

  res <- vr_api_post_json(
    OPENTARGETS_URL,
    body = list(
      query = OPENTARGETS_QUERY,
      variables = list(id = ensembl_id, size = size)
    ),
    source = "Open Targets"
  )
  if (!res$ok) {
    return(list(ok = FALSE, error = res$error))
  }
  # GraphQL reports query errors in a 200 response body.
  if (!is.null(res$data$errors)) {
    return(list(ok = FALSE, error = "Open Targets returned a query error."))
  }

  target <- pluck_at(res$data, "data", "target")
  if (is.null(target)) {
    return(list(
      ok = FALSE,
      error = paste0("Open Targets has no record for '", ensembl_id, "'.")
    ))
  }
  associated <- pluck_at(target, "associatedDiseases")
  rows <- pluck_at(associated, "rows")
  if (is.null(rows) || length(rows) == 0) {
    return(list(ok = FALSE, error = "No disease associations found."))
  }

  list(
    ok = TRUE,
    count = pluck_at(associated, "count", default = NA),
    data = opentargets_parse_rows(rows)
  )
}

# Pure parser: associated-disease rows -> data.frame(disease, disease_id, score).
opentargets_parse_rows <- function(rows) {
  data.frame(
    disease = vapply(
      rows,
      function(r) {
        as.character(pluck_at(r, "disease", "name", default = NA_character_))
      },
      character(1)
    ),
    disease_id = vapply(
      rows,
      function(r) {
        as.character(pluck_at(r, "disease", "id", default = NA_character_))
      },
      character(1)
    ),
    score = vapply(
      rows,
      function(r) {
        as.numeric(pluck_at(r, "score", default = NA))
      },
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
}
