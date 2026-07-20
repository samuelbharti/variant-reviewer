# Syntactic validation of the search inputs (gene symbol, variant identifier)
# before any external API is queried, using the {biobouncer} package. This is a
# fast, offline, format-level gate (biobouncer's "pattern" mode): it rejects
# malformed identifiers up front so a network lookup is never fired on
# obviously-bad input. It does not confirm the identifier exists — the resolving
# APIs (MyGene, MyVariant) do that.

# Guard on biobouncer being available so the app degrades gracefully (skips the
# gate rather than erroring) if the package is somehow absent.
.vr_has_biobouncer <- function() {
  requireNamespace("biobouncer", quietly = TRUE)
}

# Validate a gene symbol against the HGNC symbol grammar.
# Returns list(ok = TRUE) or list(ok = FALSE, error = <message>).
vr_validate_gene <- function(gene) {
  gene <- trimws(as.character(gene %||% ""))
  if (!nzchar(gene)) {
    return(list(ok = FALSE, error = "Enter a gene symbol."))
  }
  if (!.vr_has_biobouncer()) {
    return(list(ok = TRUE))
  }
  valid <- isTRUE(biobouncer::is_valid_id(gene, "hgnc", how = "pattern"))
  if (!valid) {
    return(list(
      ok = FALSE,
      error = sprintf(
        "“%s” is not a valid gene symbol — enter an HGNC symbol like TP53 or BRAF.",
        gene
      )
    ))
  }
  list(ok = TRUE)
}

# Validate a variant identifier. Only rsIDs are format-checked (against the
# dbSNP grammar); other accepted forms (HGVS, protein shorthand like R175H) are
# passed through for the annotation API to resolve. A blank variant is allowed
# (the variant is optional). Returns list(ok, error) as above.
vr_validate_variant <- function(variant) {
  variant <- trimws(as.character(variant %||% ""))
  if (!nzchar(variant)) {
    return(list(ok = TRUE))
  }
  looks_like_rsid <- grepl("^rs", variant, ignore.case = TRUE)
  if (.vr_has_biobouncer() && looks_like_rsid) {
    valid <- isTRUE(
      biobouncer::is_valid_id(tolower(variant), "dbsnp", how = "pattern")
    )
    if (!valid) {
      return(list(
        ok = FALSE,
        error = sprintf(
          "“%s” is not a valid dbSNP rsID — use a form like rs113488022.",
          variant
        )
      ))
    }
  }
  list(ok = TRUE)
}

# Validate a full gene/variant query. Returns list(ok = TRUE) when both pass,
# otherwise list(ok = FALSE, errors = <character vector of messages>).
vr_validate_query <- function(gene, variant = NULL) {
  errors <- character()
  g <- vr_validate_gene(gene)
  if (!isTRUE(g$ok)) {
    errors <- c(errors, g$error)
  }
  v <- vr_validate_variant(variant)
  if (!isTRUE(v$ok)) {
    errors <- c(errors, v$error)
  }
  if (length(errors) > 0) {
    return(list(ok = FALSE, errors = errors))
  }
  list(ok = TRUE)
}
