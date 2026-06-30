# Ensembl VEP client: predicted consequences of a variant.
# REST API docs: https://rest.ensembl.org/documentation/info/vep_id_get

ENSEMBL_BASE <- "https://rest.ensembl.org"

# Run VEP for a variant id (rsID).
# Returns:
#   list(ok = TRUE, most_severe, assembly,
#        data = data.frame(gene, transcript, consequence, impact, sift, polyphen))
#   list(ok = FALSE, error = "...")
ensembl_vep <- function(rsid) {
  if (is_blank(rsid)) {
    return(list(ok = FALSE, error = "No rsID available for VEP lookup."))
  }

  res <- vr_api_get(
    ENSEMBL_BASE,
    path = paste0("vep/human/id/", rsid),
    query = list(`content-type` = "application/json"),
    source = "Ensembl VEP"
  )
  if (!res$ok) {
    return(list(ok = FALSE, error = res$error))
  }
  records <- res$data
  if (is.null(records) || length(records) == 0) {
    return(list(
      ok = FALSE,
      error = paste0("Ensembl VEP has no record for ", rsid, ".")
    ))
  }

  ensembl_parse_vep(records[[1]])
}

# Pure parser: a VEP record -> normalized result.
ensembl_parse_vep <- function(record) {
  list(
    ok = TRUE,
    most_severe = pluck_at(
      record,
      "most_severe_consequence",
      default = NA_character_
    ),
    assembly = pluck_at(record, "assembly_name", default = NA_character_),
    data = ensembl_consequences_df(pluck_at(record, "transcript_consequences"))
  )
}

# Build a data.frame of protein-coding transcript consequences (the rows worth
# showing), or NULL when there are none.
ensembl_consequences_df <- function(consequences) {
  if (is.null(consequences) || length(consequences) == 0) {
    return(NULL)
  }
  is_coding <- vapply(
    consequences,
    function(x) {
      identical(pluck_at(x, "biotype"), "protein_coding")
    },
    logical(1)
  )
  consequences <- consequences[is_coding]
  if (length(consequences) == 0) {
    return(NULL)
  }

  chr_field <- function(key) {
    vapply(
      consequences,
      function(x) {
        as.character(pluck_at(x, key, default = NA_character_))
      },
      character(1)
    )
  }
  data.frame(
    gene = chr_field("gene_symbol"),
    transcript = chr_field("transcript_id"),
    consequence = vapply(
      consequences,
      function(x) {
        terms <- pluck_at(x, "consequence_terms")
        if (is.null(terms)) {
          NA_character_
        } else {
          paste(unlist(terms), collapse = ", ")
        }
      },
      character(1)
    ),
    impact = chr_field("impact"),
    sift = chr_field("sift_prediction"),
    polyphen = chr_field("polyphen_prediction"),
    stringsAsFactors = FALSE
  )
}
