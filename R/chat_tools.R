# Formatters that turn each result card's data (the list its module returns)
# into a compact text summary for the assistant's read_card tool. Kept pure and
# Shiny-free so they are unit-testable. Each takes the card's value — a list with
# ok = TRUE/FALSE, or NULL before a search — and returns a single string.

# Cards the assistant can read, id -> human label. The ids are the read_card
# enum values and the keys the dashboard snapshot is stored under (server.R).
VR_CHAT_CARDS <- c(
  gene = "Gene summary",
  variant = "Variant annotation",
  predictions = "In-silico predictions",
  protein = "Protein context (ProtVar)",
  domains = "Protein domains & features",
  structure = "3D structure (AlphaFold)",
  clinvar = "Clinical significance (ClinVar)",
  gnomad = "Population frequency (gnomAD)",
  constraint = "Gene constraint (gnomAD)",
  consequences = "Variant consequences (Ensembl VEP)",
  expression = "Tissue expression (GTEx)",
  interactions = "Protein interactions (STRING)",
  diseases = "Disease associations (Open Targets)"
)

# First value, or a default when blank.
.vr_or <- function(x, default = "n/a") {
  if (is_blank(x)) default else as.character(x)[1]
}

# Shared guard: NULL means nothing loaded; !ok means the lookup failed. Returns
# a string to short-circuit on, or NULL to continue formatting.
.vr_card_guard <- function(res, empty) {
  if (is.null(res)) {
    return(empty)
  }
  if (!isTRUE(res$ok)) {
    return(paste0("No data: ", .vr_or(res$error, "unavailable"), "."))
  }
  NULL
}

vr_chat_gene <- function(res) {
  g <- .vr_card_guard(res, "No gene loaded yet.")
  if (!is.null(g)) {
    return(g)
  }
  paste0(
    "Gene ",
    .vr_or(res$symbol),
    " — ",
    .vr_or(res$name),
    " (type: ",
    .vr_or(res$type_of_gene),
    "). ",
    "IDs: Entrez ",
    .vr_or(res$entrez),
    ", Ensembl ",
    .vr_or(res$ensembl_gene),
    ", UniProt ",
    .vr_or(res$uniprot),
    ".",
    if (!is_blank(res$summary)) paste0(" Summary: ", res$summary) else ""
  )
}

vr_chat_variant <- function(res) {
  g <- .vr_card_guard(res, "No variant loaded yet.")
  if (!is.null(g)) {
    return(g)
  }
  paste0(
    "Variant ",
    .vr_or(res$id),
    " (dbSNP ",
    .vr_or(res$rsid),
    ") in gene ",
    .vr_or(res$gene),
    ". Protein change ",
    .vr_or(res$hgvsp),
    ", CADD phred ",
    .vr_or(res$cadd_phred),
    ", ClinVar: ",
    .vr_or(res$clinvar_significance),
    "."
  )
}

vr_chat_protein <- function(res) {
  g <- .vr_card_guard(
    res,
    "No protein context yet (needs a gene with a UniProt id and a variant)."
  )
  if (!is.null(g)) {
    return(g)
  }
  n_var <- if (is.data.frame(res$variants)) nrow(res$variants) else 0L
  paste0(
    "UniProt ",
    .vr_or(res$accession),
    " position ",
    .vr_or(res$position),
    ". ",
    if (!is_blank(res$function_text)) {
      paste0("Function: ", res$function_text, " ")
    } else {
      ""
    },
    n_var,
    " catalogued variant(s) at this residue",
    if (n_var > 0) {
      paste0(": ", paste(res$variants$change, collapse = ", "), ".")
    } else {
      "."
    }
  )
}

vr_chat_clinvar <- function(res) {
  g <- .vr_card_guard(res, "No ClinVar data yet (needs an rsID/HGVS variant).")
  if (!is.null(g)) {
    return(g)
  }
  paste0(
    "ClinVar ",
    .vr_or(res$accession),
    ": significance ",
    .vr_or(res$significance),
    " (review: ",
    .vr_or(res$review_status),
    "). Condition(s): ",
    .vr_or(res$conditions),
    ". Last evaluated: ",
    .vr_or(res$last_evaluated),
    "."
  )
}

vr_chat_gnomad <- function(res) {
  g <- .vr_card_guard(res, "No gnomAD data yet (needs an rsID/HGVS variant).")
  if (!is.null(g)) {
    return(g)
  }
  part <- function(label, p) {
    if (is.null(p) || is_blank(p$af)) {
      return(NULL)
    }
    paste0(
      label,
      " AF ",
      formatC(as.numeric(p$af), format = "g", digits = 3),
      " (",
      .vr_or(p$ac),
      "/",
      .vr_or(p$an),
      " alleles)"
    )
  }
  parts <- c(part("Exome", res$exome), part("Genome", res$genome))
  if (length(parts) == 0) {
    parts <- "no allele-frequency data"
  }
  paste0(
    "Variant ",
    .vr_or(res$variant_id),
    ": ",
    paste(parts, collapse = "; "),
    "."
  )
}

vr_chat_predictions <- function(res) {
  g <- .vr_card_guard(res, "No in-silico predictions yet (needs a variant).")
  if (!is.null(g)) {
    return(g)
  }
  parts <- vapply(
    res$predictions,
    function(p) {
      score <- if (is.na(p$score)) {
        ""
      } else {
        formatC(p$score, format = "g", digits = 3)
      }
      call <- if (is.na(p$call)) "" else paste0(" ", p$call)
      trimws(gsub("\\s+", " ", paste0(p$name, " ", score, call)))
    },
    character(1)
  )
  paste0("In-silico predictions: ", paste(parts, collapse = "; "), ".")
}

vr_chat_constraint <- function(res) {
  g <- .vr_card_guard(res, "No gene constraint yet (needs a gene).")
  if (!is.null(g)) {
    return(g)
  }
  n <- function(x, d = 2) {
    if (is_blank(x)) "n/a" else formatC(x, format = "f", digits = d)
  }
  paste0(
    "Gene constraint: pLI ",
    n(res$pli, 3),
    ", LOEUF ",
    n(res$loeuf),
    ", missense Z ",
    n(res$mis_z),
    ", synonymous Z ",
    n(res$syn_z),
    ", observed/expected LoF ",
    n(res$oe_lof),
    ", observed/expected missense ",
    n(res$oe_mis),
    "."
  )
}

# Summarize the first `n` rows of a card's data frame with a per-row formatter.
.vr_rows <- function(df, n, fmt) {
  k <- min(n, nrow(df))
  vapply(seq_len(k), fmt, character(1))
}

vr_chat_domains <- function(res) {
  g <- .vr_card_guard(
    res,
    "No protein features yet (needs a gene with a UniProt id)."
  )
  if (!is.null(g)) {
    return(g)
  }
  df <- res$features
  if (!is.data.frame(df) || nrow(df) == 0) {
    return("No annotated domains or features for this protein.")
  }
  pos <- res$position
  at <- ""
  if (!is.null(pos) && !is.na(suppressWarnings(as.integer(pos)))) {
    hits <- proteins_features_at(df, pos)
    at <- if (nrow(hits) > 0) {
      labels <- vapply(
        seq_len(nrow(hits)),
        function(i) proteins_feature_text(hits, i),
        character(1)
      )
      paste0("Residue ", pos, " is in: ", paste(labels, collapse = ", "), ". ")
    } else {
      paste0("Residue ", pos, " is not within an annotated domain. ")
    }
  }
  top <- .vr_rows(df, 8, function(i) proteins_feature_text(df, i))
  paste0(at, nrow(df), " feature(s): ", paste(top, collapse = "; "), ".")
}

vr_chat_structure <- function(res) {
  g <- .vr_card_guard(
    res,
    "No 3D structure yet (needs a gene with a UniProt id)."
  )
  if (!is.null(g)) {
    return(g)
  }
  pos <- res$position
  if (!is.null(pos) && !is.na(suppressWarnings(as.integer(pos)))) {
    paste0(
      "An AlphaFold predicted structure for ",
      .vr_or(res$accession),
      " is shown, with variant residue ",
      pos,
      " highlighted."
    )
  } else {
    paste0(
      "An AlphaFold predicted structure for ",
      .vr_or(res$accession),
      " is shown."
    )
  }
}

vr_chat_consequences <- function(res) {
  g <- .vr_card_guard(
    res,
    "No VEP consequences yet (needs an rsID/HGVS variant)."
  )
  if (!is.null(g)) {
    return(g)
  }
  df <- res$data
  if (!is.data.frame(df) || nrow(df) == 0) {
    return("No protein-coding transcript consequences.")
  }
  rows <- .vr_rows(df, 5, function(i) {
    paste0(
      df$consequence[i],
      " (",
      df$transcript[i],
      ", impact ",
      df$impact[i],
      ")"
    )
  })
  paste0(
    nrow(df),
    " consequence(s). Top: ",
    paste(rows, collapse = "; "),
    "."
  )
}

vr_chat_expression <- function(res) {
  g <- .vr_card_guard(res, "No GTEx expression yet (needs a gene).")
  if (!is.null(g)) {
    return(g)
  }
  df <- res$data
  if (!is.data.frame(df) || nrow(df) == 0) {
    return("No expression data.")
  }
  ord <- order(df$median_tpm, decreasing = TRUE)
  rows <- .vr_rows(df[ord, , drop = FALSE], 5, function(i) {
    paste0(
      df[ord, ][["tissue"]][i],
      " ",
      formatC(df[ord, ][["median_tpm"]][i], format = "f", digits = 1),
      " TPM"
    )
  })
  paste0(
    "Median expression across ",
    nrow(df),
    " tissues. Highest: ",
    paste(rows, collapse = "; "),
    "."
  )
}

vr_chat_interactions <- function(res) {
  g <- .vr_card_guard(res, "No STRING interactions yet (needs a gene).")
  if (!is.null(g)) {
    return(g)
  }
  df <- res$data
  if (!is.data.frame(df) || nrow(df) == 0) {
    return("No interaction partners.")
  }
  rows <- .vr_rows(df, 8, function(i) {
    paste0(
      df$partner[i],
      " (",
      formatC(df$score[i], format = "f", digits = 3),
      ")"
    )
  })
  paste0(
    nrow(df),
    " partner(s). Top by combined score: ",
    paste(rows, collapse = ", "),
    "."
  )
}

vr_chat_diseases <- function(res) {
  g <- .vr_card_guard(res, "No disease associations yet (needs a gene).")
  if (!is.null(g)) {
    return(g)
  }
  df <- res$data
  if (!is.data.frame(df) || nrow(df) == 0) {
    return("No disease associations.")
  }
  rows <- .vr_rows(df, 8, function(i) {
    paste0(
      df$disease[i],
      " (",
      formatC(df$score[i], format = "f", digits = 3),
      ")"
    )
  })
  extra <- if (!is_blank(res$count)) {
    paste0(" (top ", nrow(df), " of ", res$count, ")")
  } else {
    ""
  }
  paste0(
    "Disease associations",
    extra,
    ": ",
    paste(rows, collapse = "; "),
    "."
  )
}

# Dispatch: format one card's data by id. Unknown ids report themselves rather
# than error, so a bad tool argument degrades gracefully.
vr_chat_card_text <- function(card, data) {
  fn <- switch(
    card,
    gene = vr_chat_gene,
    variant = vr_chat_variant,
    predictions = vr_chat_predictions,
    protein = vr_chat_protein,
    domains = vr_chat_domains,
    structure = vr_chat_structure,
    clinvar = vr_chat_clinvar,
    gnomad = vr_chat_gnomad,
    constraint = vr_chat_constraint,
    consequences = vr_chat_consequences,
    expression = vr_chat_expression,
    interactions = vr_chat_interactions,
    diseases = vr_chat_diseases,
    NULL
  )
  if (is.null(fn)) {
    return(paste0(
      "Unknown card '",
      card,
      "'. Available: ",
      paste(names(VR_CHAT_CARDS), collapse = ", "),
      "."
    ))
  }
  fn(data)
}
