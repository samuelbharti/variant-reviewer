# "View on source" links shown in each card header, so every card points to the
# authoritative external page for what it displays. The URL builders are pure and
# return NULL when the required identifier is missing.

# Small right-aligned external link for a card header (or NULL if no href).
vr_source_link <- function(href, label = "Source") {
  if (is_blank(href)) {
    return(NULL)
  }
  tags$a(
    href = href,
    target = "_blank",
    rel = "noopener noreferrer",
    class = "small text-decoration-none text-nowrap",
    label,
    " ",
    icon("up-right-from-square")
  )
}

# Card header with the title on the left and the module's source-link slot
# (`uiOutput(ns("source"))`) on the right. `ns` is the module's NS().
vr_card_header <- function(title, ns) {
  card_header(
    class = "d-flex justify-content-between align-items-center gap-2",
    tags$span(title),
    uiOutput(ns("source"), inline = TRUE)
  )
}

# --- Per-source URL builders --------------------------------------------------

src_ncbi_gene <- function(entrez) {
  if (is_blank(entrez)) {
    NULL
  } else {
    paste0("https://www.ncbi.nlm.nih.gov/gene/", entrez)
  }
}

src_dbsnp <- function(rsid) {
  if (is_blank(rsid)) {
    NULL
  } else {
    paste0("https://www.ncbi.nlm.nih.gov/snp/", rsid)
  }
}

src_gtex <- function(symbol) {
  if (is_blank(symbol)) {
    NULL
  } else {
    paste0("https://gtexportal.org/home/gene/", symbol)
  }
}

src_string <- function(symbol) {
  if (is_blank(symbol)) {
    return(NULL)
  }
  paste0(
    "https://string-db.org/cgi/network?identifiers=",
    utils::URLencode(symbol, reserved = TRUE),
    "&species=9606"
  )
}

src_opentargets_gene <- function(ensembl) {
  if (is_blank(ensembl)) {
    NULL
  } else {
    paste0("https://platform.opentargets.org/target/", ensembl)
  }
}

src_gnomad_gene <- function(ensembl) {
  if (is_blank(ensembl)) {
    NULL
  } else {
    paste0(
      "https://gnomad.broadinstitute.org/gene/",
      ensembl,
      "?dataset=gnomad_r4"
    )
  }
}

src_gnomad_variant <- function(variant_id, dataset = "gnomad_r4") {
  if (is_blank(variant_id)) {
    NULL
  } else {
    paste0(
      "https://gnomad.broadinstitute.org/variant/",
      variant_id,
      "?dataset=",
      dataset
    )
  }
}

src_uniprot <- function(uniprot, section = NULL) {
  if (is_blank(uniprot)) {
    return(NULL)
  }
  paste0(
    "https://www.uniprot.org/uniprotkb/",
    uniprot,
    "/entry",
    if (!is.null(section)) paste0("#", section) else ""
  )
}

src_alphafold <- function(uniprot) {
  if (is_blank(uniprot)) {
    NULL
  } else {
    paste0("https://alphafold.ebi.ac.uk/entry/", uniprot)
  }
}

src_ensembl_variant <- function(rsid) {
  if (is_blank(rsid)) {
    NULL
  } else {
    paste0("https://www.ensembl.org/Homo_sapiens/Variation/Explore?v=", rsid)
  }
}

src_ensembl_gene <- function(ensembl_gene) {
  if (is_blank(ensembl_gene)) {
    NULL
  } else {
    paste0("https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=", ensembl_gene)
  }
}

src_clinvar_variation <- function(uid) {
  if (is_blank(uid)) {
    NULL
  } else {
    paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/", uid, "/")
  }
}

src_monarch_gene <- function(hgnc) {
  if (is_blank(hgnc)) {
    return(NULL)
  }
  id <- if (grepl("^HGNC:", hgnc, ignore.case = TRUE)) {
    toupper(hgnc)
  } else {
    paste0("HGNC:", hgnc)
  }
  paste0("https://monarchinitiative.org/", id)
}
