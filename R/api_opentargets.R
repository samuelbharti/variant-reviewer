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

# Known drugs and clinical candidates for a target. The field returns one row
# per drug with its highest clinical stage and the diseases it has been tried
# against; the element type carries `diseaseFromSource` as a fallback label when
# the mapped disease name is blank.
OPENTARGETS_DRUGS_QUERY <- paste(
  "query($id: String!) {",
  "  target(ensemblId: $id) {",
  "    drugAndClinicalCandidates {",
  "      count",
  "      rows {",
  "        maxClinicalStage",
  "        drug { id name drugType }",
  "        diseases { diseaseFromSource disease { id name } }",
  "      }",
  "    }",
  "  }",
  "}",
  sep = "\n"
)

# Known drugs / clinical candidates for an Ensembl gene id (target).
# Returns:
#   list(ok = TRUE, count,
#        data = data.frame(drug, drug_id, drug_type, max_phase, disease))
#   list(ok = FALSE, error = "...")
opentargets_drugs <- function(ensembl_id) {
  if (is_blank(ensembl_id)) {
    return(list(
      ok = FALSE,
      error = "No Ensembl gene ID available for this gene."
    ))
  }

  res <- vr_api_post_json(
    OPENTARGETS_URL,
    body = list(
      query = OPENTARGETS_DRUGS_QUERY,
      variables = list(id = ensembl_id)
    ),
    source = "Open Targets"
  )
  if (!res$ok) {
    return(list(ok = FALSE, error = res$error))
  }
  if (!is.null(res$data$errors)) {
    return(list(ok = FALSE, error = "Open Targets returned a query error."))
  }

  known <- pluck_at(res$data, "data", "target", "drugAndClinicalCandidates")
  rows <- pluck_at(known, "rows")
  if (is.null(rows) || length(rows) == 0) {
    return(list(
      ok = FALSE,
      error = "No known drugs or clinical candidates found."
    ))
  }

  list(
    ok = TRUE,
    count = pluck_at(known, "count", default = NA),
    data = opentargets_parse_drugs(rows)
  )
}

# Pure parser: known-drug rows ->
# data.frame(drug, drug_id, drug_type, max_phase, disease).
opentargets_parse_drugs <- function(rows) {
  first_disease <- function(r) {
    diseases <- pluck_at(r, "diseases")
    if (is.null(diseases) || length(diseases) == 0) {
      return(NA_character_)
    }
    d <- diseases[[1]]
    as.character(pluck_at(
      d,
      "disease",
      "name",
      default = pluck_at(d, "diseaseFromSource", default = NA_character_)
    ))
  }
  data.frame(
    drug = vapply(
      rows,
      function(r) {
        as.character(pluck_at(r, "drug", "name", default = NA_character_))
      },
      character(1)
    ),
    drug_id = vapply(
      rows,
      function(r) {
        as.character(pluck_at(r, "drug", "id", default = NA_character_))
      },
      character(1)
    ),
    drug_type = vapply(
      rows,
      function(r) {
        as.character(pluck_at(r, "drug", "drugType", default = NA_character_))
      },
      character(1)
    ),
    max_phase = vapply(
      rows,
      function(r) {
        opentargets_pretty_phase(
          pluck_at(r, "maxClinicalStage", default = NA_character_)
        )
      },
      character(1)
    ),
    disease = vapply(rows, first_disease, character(1)),
    stringsAsFactors = FALSE
  )
}

# "PHASE_2" -> "Phase 2", "PRE_CLINICAL" -> "Pre clinical".
opentargets_pretty_phase <- function(x) {
  if (is_blank(x)) {
    return(NA_character_)
  }
  s <- gsub("_", " ", tolower(as.character(x)))
  paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
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
