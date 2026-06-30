# ProtVar (EBI) client: protein-level functional context and known variants at
# a residue. Docs: https://www.ebi.ac.uk/ProtVar/
#
# ProtVar exposes per-position endpoints keyed by UniProt accession + protein
# position: /function/{acc}/{pos} and /population/{acc}/{pos}. We get the
# accession from MyGene and the position by parsing the variant string (or a
# MyVariant hgvsp such as "p.Arg175His").

PROTVAR_BASE <- "https://www.ebi.ac.uk/ProtVar/api"

# Extract a protein position (first integer) from a protein-change string such
# as "R175H", "p.Arg175His", "Arg175His", or "175". Returns an integer or NULL.
protvar_parse_position <- function(variant) {
  if (is_blank(variant)) {
    return(NULL)
  }
  m <- regmatches(variant, regexpr("\\d+", variant))
  if (length(m) == 0 || m == "") {
    return(NULL)
  }
  as.integer(m)
}

# Combine ProtVar function context + variants at a residue into one result.
# Returns:
#   list(ok = TRUE, accession, position, function_text,
#        variants = data.frame(change, sources))
#   list(ok = FALSE, error = "...")
protvar_annotate <- function(accession, position) {
  if (is_blank(accession)) {
    return(list(ok = FALSE, error = "No UniProt accession available."))
  }
  if (is_blank(position)) {
    return(list(ok = FALSE, error = "No protein position available."))
  }
  path_pos <- paste0(accession, "/", position)

  fn <- vr_api_get(
    PROTVAR_BASE,
    path = paste0("function/", path_pos),
    source = "ProtVar"
  )
  if (!fn$ok) {
    return(list(ok = FALSE, error = fn$error))
  }

  pop <- vr_api_get(
    PROTVAR_BASE,
    path = paste0("population/", path_pos),
    source = "ProtVar"
  )

  list(
    ok = TRUE,
    accession = accession,
    position = position,
    function_text = protvar_function_text(fn$data),
    variants = if (pop$ok) protvar_variants_df(pop$data) else NULL
  )
}

# Pull the first FUNCTION comment text from a /function response.
protvar_function_text <- function(data) {
  comments <- pluck_at(data, "comments")
  if (is.null(comments)) {
    return(NA_character_)
  }
  for (cm in comments) {
    if (identical(pluck_at(cm, "type"), "FUNCTION")) {
      txt <- pluck_at(cm, "text")
      if (!is.null(txt) && length(txt) > 0) {
        return(as.character(pluck_at(
          txt[[1]],
          "value",
          default = NA_character_
        )))
      }
    }
  }
  NA_character_
}

# Summarize known variants at the residue into a data.frame(change, sources).
protvar_variants_df <- function(data) {
  variants <- pluck_at(data, "variants")
  if (is.null(variants) || length(variants) == 0) {
    return(NULL)
  }
  change <- vapply(
    variants,
    function(v) {
      as.character(pluck_at(v, "alternativeSequence", default = NA_character_))
    },
    character(1)
  )
  sources <- vapply(
    variants,
    function(v) {
      xrefs <- pluck_at(v, "xrefs")
      if (is.null(xrefs)) {
        return(NA_character_)
      }
      names <- unique(vapply(
        xrefs,
        function(x) {
          as.character(pluck_at(x, "name", default = NA_character_))
        },
        character(1)
      ))
      names <- names[!is.na(names)]
      if (length(names) == 0) NA_character_ else paste(names, collapse = ", ")
    },
    character(1)
  )

  df <- data.frame(change = change, sources = sources, stringsAsFactors = FALSE)
  df[!is.na(df$change), , drop = FALSE]
}
