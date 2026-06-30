# gnomAD client: population allele frequencies for a variant, looked up by
# rsID (avoids hg19/hg38 coordinate mismatches). GraphQL API:
# https://gnomad.broadinstitute.org/api

GNOMAD_URL <- "https://gnomad.broadinstitute.org/api"
GNOMAD_DATASET <- "gnomad_r4"

# Allele frequencies for an rsID.
# Returns:
#   list(ok = TRUE, variant_id, dataset, exome = list(af, ac, an)|NULL,
#        genome = list(af, ac, an)|NULL)
#   list(ok = FALSE, error = "...")
gnomad_frequency <- function(rsid, dataset = GNOMAD_DATASET) {
  if (is_blank(rsid)) {
    return(list(ok = FALSE, error = "No rsID available for gnomAD lookup."))
  }

  query <- sprintf(
    paste(
      "query($rsid: String!) {",
      "  variant(rsid: $rsid, dataset: %s) {",
      "    variant_id",
      "    exome { af ac an }",
      "    genome { af ac an }",
      "  }",
      "}",
      sep = "\n"
    ),
    dataset
  )

  res <- vr_api_post_json(
    GNOMAD_URL,
    body = list(query = query, variables = list(rsid = rsid)),
    source = "gnomAD"
  )
  if (!res$ok) {
    return(list(ok = FALSE, error = res$error))
  }
  if (!is.null(res$data$errors)) {
    return(list(ok = FALSE, error = "gnomAD returned a query error."))
  }

  variant <- pluck_at(res$data, "data", "variant")
  if (is.null(variant)) {
    return(list(
      ok = FALSE,
      error = paste0("gnomAD has no record for ", rsid, ".")
    ))
  }

  list(
    ok = TRUE,
    variant_id = pluck_at(variant, "variant_id", default = NA_character_),
    dataset = dataset,
    exome = gnomad_freq_part(pluck_at(variant, "exome")),
    genome = gnomad_freq_part(pluck_at(variant, "genome"))
  )
}

# Normalize one frequency block (exome or genome) to a numeric list, or NULL
# when gnomAD has no data for that sample set.
gnomad_freq_part <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  list(
    af = as.numeric(pluck_at(x, "af", default = NA)),
    ac = as.numeric(pluck_at(x, "ac", default = NA)),
    an = as.numeric(pluck_at(x, "an", default = NA))
  )
}
