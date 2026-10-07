# MyVariant.info client: annotate a variant given an rsID or HGVS string.
# Docs: https://docs.myvariant.info/en/latest/

MYVARIANT_BASE <- "https://myvariant.info/v1"

# The dashboard annotates GRCh38, so variant lookups ask MyVariant for hg38
# records. Left out, MyVariant answers in hg19.
MYVARIANT_ASSEMBLY <- "hg38"

# One rsID can cover several alleles (rs113488022 is V600A, V600E and V600G),
# so a variant lookup asks for up to this many hits and then picks one.
MYVARIANT_MAX_ALLELES <- 10

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

# The `q` value for a variant. MyVariant reads an unquoted HGVS string such as
# chr7:g.140753336A>T as a field:value query and finds nothing, so HGVS is
# quoted.
myvariant_query_term <- function(variant) {
  term <- trimws(as.character(variant))
  if (grepl(":", term, fixed = TRUE)) paste0("\"", term, "\"") else term
}

# Fetch the hit for the one allele `variant` names. An HGVS string names one
# allele. An rsID can name several, and then no hit is picked: the error lists
# the alleles so the user can choose, since any choice made here would be a
# guess (taking the first hit is what showed V600A for BRAF V600E).
# `fields` are the fields the caller needs; the ones used to label alleles are
# added here.
# Returns:
#   list(ok = TRUE, hit)
#   list(ok = FALSE, error, ambiguous = TRUE, gene, alleles) for several alleles
#   list(ok = FALSE, error = "...") otherwise
myvariant_fetch_allele <- function(variant, fields, not_found) {
  term <- trimws(as.character(variant))
  allele_fields <- c(
    "chrom",
    "vcf",
    "clinvar.variant_id",
    "dbnsfp.genename",
    "dbnsfp.hgvsp",
    "snpeff.ann.hgvs_p"
  )
  res <- vr_api_get(
    MYVARIANT_BASE,
    path = "query",
    query = list(
      q = myvariant_query_term(term),
      size = MYVARIANT_MAX_ALLELES,
      assembly = MYVARIANT_ASSEMBLY,
      fields = paste(unique(c(fields, allele_fields)), collapse = ",")
    ),
    source = "MyVariant"
  )
  if (!res$ok) {
    return(list(ok = FALSE, error = res$error))
  }
  myvariant_pick_allele(res$data$hits, term, not_found)
}

# Pure helper for myvariant_fetch_allele(): the one hit, or an error naming
# every allele when there are several.
myvariant_pick_allele <- function(hits, term, not_found) {
  if (is.null(hits) || length(hits) == 0) {
    return(list(ok = FALSE, error = paste0(not_found, " '", term, "'.")))
  }
  hits <- myvariant_distinct_alleles(hits)
  if (length(hits) == 1) {
    return(list(ok = TRUE, hit = hits[[1]]))
  }
  ids <- vapply(
    hits,
    function(h) as.character(pluck_at(h, "_id", default = NA_character_)),
    character(1)
  )
  changes <- vapply(hits, myvariant_hgvsp, character(1))
  alleles <- data.frame(id = ids, hgvsp = changes, stringsAsFactors = FALSE)
  alleles <- alleles[order(alleles$id), , drop = FALSE]
  rownames(alleles) <- NULL
  listed <- ifelse(
    is.na(alleles$hgvsp),
    alleles$id,
    paste0(alleles$hgvsp, " (", alleles$id, ")")
  )
  genes <- unique(stats::na.omit(vapply(
    hits,
    function(h) mygene_first(pluck_at(h, "dbnsfp", "genename")),
    character(1)
  )))
  list(
    ok = FALSE,
    ambiguous = TRUE,
    # Every allele of an rsID sits in the same gene, so the gene-level cards
    # can still load while the user picks one.
    gene = if (length(genes) == 1) genes[[1]] else NULL,
    alleles = alleles,
    error = paste0(
      term,
      " covers ",
      nrow(alleles),
      " alleles: ",
      paste(listed, collapse = ", "),
      ". Pick one from the Variant list, or enter its HGVS."
    )
  )
}

# MyVariant can hold one indel twice, shifted left and right inside a repeat:
# CFTR F508del is both 7-117559590-ATCT-A and 7-117559591-TCTT-T. Those are one
# allele, so they are merged, keeping the hit with a ClinVar record (the one
# ClinVar and gnomAD also use). Hits without VCF fields are kept as they are.
myvariant_distinct_alleles <- function(hits) {
  kept <- list()
  for (hit in hits) {
    vcf <- myvariant_vcf(hit)
    same <- which(vapply(
      kept,
      function(k) .mv_same_change(myvariant_vcf(k), vcf),
      logical(1)
    ))
    if (length(same) == 0) {
      kept[[length(kept) + 1]] <- hit
    } else if (
      is_blank(pluck_at(kept[[same[[1]]]], "clinvar", "variant_id")) &&
        !is_blank(pluck_at(hit, "clinvar", "variant_id"))
    ) {
      kept[[same[[1]]]] <- hit
    }
  }
  kept
}

# TRUE when two VCF records (see myvariant_vcf()) make the same edit to the
# reference. Their ref bases are laid over one stretch of the reference (they
# must agree where they overlap and leave no gap), each edit is applied, and
# the two results are compared.
.mv_same_change <- function(a, b) {
  if (is.null(a) || is.null(b) || !identical(a$chrom, b$chrom)) {
    return(FALSE)
  }
  start <- min(a$pos, b$pos)
  end <- max(a$pos + nchar(a$ref), b$pos + nchar(b$ref)) - 1
  ref <- rep(NA_character_, end - start + 1)
  for (r in list(a, b)) {
    at <- r$pos - start + seq_len(nchar(r$ref))
    bases <- strsplit(r$ref, "", fixed = TRUE)[[1]]
    if (any(!is.na(ref[at]) & ref[at] != bases)) {
      return(FALSE)
    }
    ref[at] <- bases
  }
  if (anyNA(ref)) {
    return(FALSE)
  }
  edit <- function(r) {
    before <- seq_len(r$pos - start)
    after <- seq_along(ref) > r$pos - start + nchar(r$ref)
    paste0(
      paste(ref[before], collapse = ""),
      r$alt,
      paste(ref[after], collapse = "")
    )
  }
  identical(edit(a), edit(b))
}

# Returns:
#   list(ok = TRUE, id, rsid, gene, hgvsp, hgvsp_all, cadd_phred,
#        clinvar_significance, clinvar_id, vcf_id)
#   list(ok = FALSE, error = "...") (plus ambiguous, gene and alleles when an
#   rsID covers several alleles; see myvariant_fetch_allele())
myvariant_annotate <- function(variant) {
  if (is_blank(variant)) {
    return(list(ok = FALSE, error = "No variant supplied."))
  }
  if (!myvariant_is_queryable(variant)) {
    return(list(
      ok = FALSE,
      error = "Enter an rsID (rs...) or HGVS (e.g. chr7:g.140753336A>T) for variant-level annotation."
    ))
  }
  term <- trimws(as.character(variant))

  found <- myvariant_fetch_allele(
    term,
    fields = c(
      "dbsnp.rsid",
      "cadd.phred",
      "dbnsfp.cadd.phred",
      "clinvar.rcv.clinical_significance"
    ),
    not_found = "No annotation found for"
  )
  if (!found$ok) {
    return(found)
  }
  myvariant_parse_hit(found$hit, term)
}

# The allele the variant cards describe, from a myvariant_annotate() result: the
# result itself when it names one allele, list(ambiguous = TRUE) when the input
# is an rsID covering several, and NULL when there is nothing to go on.
vr_variant_allele <- function(annotation) {
  if (isTRUE(annotation$ok)) {
    return(annotation)
  }
  if (isTRUE(annotation$ambiguous)) {
    return(list(ambiguous = TRUE))
  }
  NULL
}

# What an allele-level card shows while an rsID's allele is not picked yet.
vr_allele_needed <- function(what) {
  list(
    ok = FALSE,
    error = paste0(
      "This rsID covers more than one allele. Pick one in the Variant box to see its ",
      what,
      "."
    )
  )
}

# Pure parser: turn a single MyVariant hit into the normalized result. Besides
# what the Variant card shows, it carries what the other variant cards need to
# find this same allele: its ClinVar variation id, its VCF id (for gnomAD and
# VEP), and every spelling of its protein change.
myvariant_parse_hit <- function(hit, term = NA_character_) {
  list(
    ok = TRUE,
    id = pluck_at(hit, "_id", default = term),
    rsid = mygene_first(pluck_at(hit, "dbsnp", "rsid")),
    gene = mygene_first(pluck_at(hit, "dbnsfp", "genename")),
    hgvsp = myvariant_hgvsp(hit),
    hgvsp_all = myvariant_hgvsp_all(hit),
    cadd_phred = .mv_cadd(hit),
    clinvar_significance = myvariant_clinvar_sig(hit),
    clinvar_id = mygene_first(pluck_at(hit, "clinvar", "variant_id")),
    vcf_id = myvariant_vcf_id(hit)
  )
}

# snpEff's protein changes for a hit. `snpeff.ann` is one object, or a list of
# them when the variant hits several transcripts.
.mv_snpeff_hgvsp <- function(hit) {
  ann <- pluck_at(hit, "snpeff", "ann")
  if (is.null(ann)) {
    return(character())
  }
  if (!is.null(names(ann))) {
    ann <- list(ann)
  }
  changes <- unlist(
    lapply(ann, function(a) pluck_at(a, "hgvs_p")),
    use.names = FALSE
  )
  changes[!is.na(changes) & nzchar(changes)]
}

# The protein change to show for a hit (e.g. p.Val600Glu). snpEff's comes
# first: it is on the RefSeq transcript. dbNSFP lists one per Ensembl
# transcript in no set order, so its first entry can be a minor isoform
# (p.Val640Glu for BRAF V600E).
myvariant_hgvsp <- function(hit) {
  snpeff <- .mv_snpeff_hgvsp(hit)
  if (length(snpeff) > 0) {
    return(snpeff[[1]])
  }
  mygene_first(pluck_at(hit, "dbnsfp", "hgvsp"))
}

# Every three-letter protein change for a hit, across transcripts. ClinVar
# names a record by one of them, e.g. "NM_004333.6(BRAF):c.1799T>A
# (p.Val600Glu)", so this is what a ClinVar title is matched against.
myvariant_hgvsp_all <- function(hit) {
  changes <- c(
    .mv_snpeff_hgvsp(hit),
    unlist(pluck_at(hit, "dbnsfp", "hgvsp"), use.names = FALSE)
  )
  unique(changes[grepl("^p\\.[A-Z][a-z]{2}[0-9]", changes)])
}

# The hit's GRCh38 VCF record as list(chrom, pos, ref, alt), or NULL when a
# field is missing.
myvariant_vcf <- function(hit) {
  chrom <- as.character(pluck_at(hit, "chrom", default = NA))
  pos <- suppressWarnings(as.integer(pluck_at(hit, "vcf", "position")))
  ref <- as.character(pluck_at(hit, "vcf", "ref", default = NA))
  alt <- as.character(pluck_at(hit, "vcf", "alt", default = NA))
  parts <- c(chrom, ref, alt)
  if (length(pos) != 1 || is.na(pos) || anyNA(parts) || !all(nzchar(parts))) {
    return(NULL)
  }
  list(chrom = chrom, pos = pos, ref = ref, alt = alt)
}

# The allele as chrom-pos-ref-alt ("7-140753336-A-T"), the form gnomAD uses for
# its variant ids, or NA when the VCF fields are missing.
myvariant_vcf_id <- function(hit) {
  vcf <- myvariant_vcf(hit)
  if (is.null(vcf)) {
    return(NA_character_)
  }
  paste(vcf$chrom, vcf$pos, vcf$ref, vcf$alt, sep = "-")
}

# CADD phred for a hit. MyVariant's own CADD block exists only for hg19
# records, so hg38 records fall back to the CADD score dbNSFP carries.
.mv_cadd <- function(hit) {
  cadd <- .mv_max_num(pluck_at(hit, "cadd", "phred"))
  if (is.na(cadd)) {
    cadd <- .mv_max_num(pluck_at(hit, "dbnsfp", "cadd", "phred"))
  }
  cadd
}

# In-silico pathogenicity predictions for a variant, from dbNSFP (+ CADD) via
# MyVariant. Fetched through myvariant_fetch_allele() like
# myvariant_annotate(), so it resolves the same allele. Returns:
#   list(ok = TRUE, predictions = list(list(name, score, call), ...))
#   list(ok = FALSE, error = "...")
myvariant_predictions <- function(variant) {
  if (is_blank(variant)) {
    return(list(ok = FALSE, error = "No variant supplied."))
  }
  if (!myvariant_is_queryable(variant)) {
    return(list(
      ok = FALSE,
      error = "Enter an rsID (rs...) or HGVS for in-silico predictions."
    ))
  }
  found <- myvariant_fetch_allele(
    variant,
    fields = c(
      "cadd.phred",
      "dbnsfp.cadd.phred",
      "dbnsfp.revel",
      "dbnsfp.alphamissense",
      "dbnsfp.sift",
      "dbnsfp.polyphen2",
      "dbnsfp.metalr",
      "dbnsfp.metasvm"
    ),
    not_found = "No predictions found for"
  )
  if (!found$ok) {
    return(list(ok = FALSE, error = found$error))
  }
  myvariant_parse_predictions(found$hit)
}

# dbNSFP prediction-code dictionaries (per predictor). Codes come as a scalar or
# a per-transcript array; the parser collapses them to one representative call.
.mv_pred_maps <- list(
  alphamissense = c(
    P = "likely pathogenic",
    B = "likely benign",
    A = "ambiguous"
  ),
  polyphen2 = c(D = "probably damaging", P = "possibly damaging", B = "benign"),
  sift = c(D = "deleterious", T = "tolerated"),
  meta = c(D = "damaging", T = "tolerated")
)

# Max numeric across a scalar/array (most-damaging), NA if none.
.mv_max_num <- function(x) {
  if (is.null(x)) {
    return(NA_real_)
  }
  v <- suppressWarnings(as.numeric(unlist(x, use.names = FALSE)))
  v <- v[!is.na(v)]
  if (length(v) == 0) NA_real_ else max(v)
}

# First non-empty prediction code mapped through `map`, NA if none.
.mv_call <- function(x, map) {
  codes <- unlist(x, use.names = FALSE)
  codes <- codes[nzchar(codes)]
  if (length(codes) == 0) {
    return(NA_character_)
  }
  code <- codes[[1]]
  if (code %in% names(map)) unname(map[[code]]) else code
}

# Pure parser: turn a MyVariant hit into an ordered list of predictor readouts.
myvariant_parse_predictions <- function(hit) {
  d <- pluck_at(hit, "dbnsfp")
  revel <- .mv_max_num(pluck_at(d, "revel", "score"))
  cadd <- .mv_cadd(hit)
  entries <- list(
    list(
      name = "REVEL",
      score = revel,
      call = if (is.na(revel)) {
        NA_character_
      } else if (revel >= 0.5) {
        "damaging-leaning"
      } else {
        "benign-leaning"
      }
    ),
    list(
      name = "AlphaMissense",
      score = .mv_max_num(pluck_at(d, "alphamissense", "score")),
      call = .mv_call(
        pluck_at(d, "alphamissense", "pred"),
        .mv_pred_maps$alphamissense
      )
    ),
    list(
      name = "CADD (phred)",
      score = cadd,
      call = if (!is.na(cadd) && cadd >= 20) {
        "top ~1% deleterious"
      } else {
        NA_character_
      }
    ),
    list(
      name = "PolyPhen-2",
      score = .mv_max_num(pluck_at(d, "polyphen2", "hdiv", "score")),
      call = .mv_call(
        pluck_at(d, "polyphen2", "hdiv", "pred"),
        .mv_pred_maps$polyphen2
      )
    ),
    list(
      name = "SIFT",
      score = .mv_max_num(pluck_at(d, "sift", "score")),
      call = .mv_call(pluck_at(d, "sift", "pred"), .mv_pred_maps$sift)
    ),
    list(
      name = "MetaLR",
      score = .mv_max_num(pluck_at(d, "metalr", "score")),
      call = .mv_call(pluck_at(d, "metalr", "pred"), .mv_pred_maps$meta)
    ),
    list(
      name = "MetaSVM",
      score = .mv_max_num(pluck_at(d, "metasvm", "score")),
      call = .mv_call(pluck_at(d, "metasvm", "pred"), .mv_pred_maps$meta)
    )
  )
  entries <- Filter(function(e) !is.na(e$score) || !is.na(e$call), entries)
  if (length(entries) == 0) {
    return(list(
      ok = FALSE,
      error = "No in-silico predictions available for this variant."
    ))
  }
  list(ok = TRUE, predictions = entries)
}

# Evolutionary conservation scores for a variant's position, from dbNSFP via
# MyVariant. Higher scores/ranks mean a more conserved (less tolerant) position,
# a supporting line of evidence in variant interpretation. The whole dbnsfp
# block is fetched because the "gerp++" key cannot be requested through the
# field selector. Returns:
#   list(ok = TRUE, metrics = data.frame(metric, score, rankscore))
#   list(ok = FALSE, error = "...")
myvariant_conservation <- function(variant) {
  if (is_blank(variant)) {
    return(list(ok = FALSE, error = "No variant supplied."))
  }
  if (!myvariant_is_queryable(variant)) {
    return(list(
      ok = FALSE,
      error = "Enter an rsID (rs...) or HGVS for conservation scores."
    ))
  }
  term <- trimws(as.character(variant))
  res <- vr_api_get(
    MYVARIANT_BASE,
    path = "query",
    # Conservation is a property of the position, so whichever allele of an
    # rsID comes first gives the same scores.
    query = list(
      q = term,
      size = 1,
      assembly = MYVARIANT_ASSEMBLY,
      fields = "dbnsfp"
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
      error = paste0("No conservation scores found for '", term, "'.")
    ))
  }
  myvariant_parse_conservation(hits[[1]])
}

# One conservation-metric row (numeric score + 0-1 rankscore, NA when absent).
.mv_cons_row <- function(metric, score, rankscore) {
  data.frame(
    metric = metric,
    score = suppressWarnings(as.numeric(score %||% NA)),
    rankscore = suppressWarnings(as.numeric(rankscore %||% NA)),
    stringsAsFactors = FALSE
  )
}

# Pure parser: pull the four common conservation metrics out of a dbNSFP hit.
myvariant_parse_conservation <- function(hit) {
  d <- pluck_at(hit, "dbnsfp")
  rows <- list(
    .mv_cons_row(
      "phyloP (100-way vertebrate)",
      pluck_at(d, "phylop", "100way_vertebrate", "score"),
      pluck_at(d, "phylop", "100way_vertebrate", "rankscore")
    ),
    .mv_cons_row(
      "phastCons (100-way vertebrate)",
      pluck_at(d, "phastcons", "100way_vertebrate", "score"),
      pluck_at(d, "phastcons", "100way_vertebrate", "rankscore")
    ),
    .mv_cons_row(
      "GERP++ RS",
      pluck_at(d, "gerp++", "rs"),
      pluck_at(d, "gerp++", "rs_rankscore")
    ),
    .mv_cons_row(
      "SiPhy (29-way)",
      pluck_at(d, "siphy_29way", "logodds_score"),
      pluck_at(d, "siphy_29way", "logodds_rankscore")
    )
  )
  df <- do.call(rbind, rows)
  df <- df[!(is.na(df$score) & is.na(df$rankscore)), , drop = FALSE]
  if (nrow(df) == 0) {
    return(list(
      ok = FALSE,
      error = "No conservation scores available for this variant."
    ))
  }
  rownames(df) <- NULL
  list(ok = TRUE, metrics = df)
}

# Notable variants for a gene: ClinVar pathogenic / likely-pathogenic variants
# that carry an rsID, used to populate the search box's variant suggestions.
# One row per allele, keyed by its hg38 HGVS id, so the two pathogenic alleles
# of rs113488022 (V600E and V600G) are two separate choices.
# Returns:
#   list(ok = TRUE, variants = data.frame(id, rsid, label, significance, cadd))
#   list(ok = FALSE, error = "...")
myvariant_gene_variants <- function(symbol, size = 200) {
  if (is_blank(symbol)) {
    return(list(ok = FALSE, error = "No gene supplied."))
  }
  sym <- trimws(as.character(symbol))
  res <- vr_api_get(
    MYVARIANT_BASE,
    path = "query",
    query = list(
      q = paste0(
        "clinvar.gene.symbol:",
        sym,
        " AND clinvar.rcv.clinical_significance:",
        "(\"Pathogenic\" OR \"Likely pathogenic\")",
        " AND _exists_:dbsnp.rsid"
      ),
      size = size,
      assembly = MYVARIANT_ASSEMBLY,
      fields = paste(
        "dbsnp.rsid",
        "dbnsfp.aa.ref",
        "dbnsfp.aa.alt",
        "dbnsfp.aa.pos",
        "clinvar.rcv.clinical_significance",
        "cadd.phred",
        "dbnsfp.cadd.phred",
        sep = ","
      )
    ),
    source = "MyVariant"
  )
  if (!res$ok) {
    return(list(ok = FALSE, error = res$error))
  }
  myvariant_parse_gene_variants(res$data$hits)
}

# One-letter amino-acid change (e.g. "V600E") from a dbnsfp.aa block, or NA.
# `pos` arrives as a per-transcript array; the first position is representative.
.mv_aa_label <- function(aa) {
  ref <- mygene_first(pluck_at(aa, "ref"))
  alt <- mygene_first(pluck_at(aa, "alt"))
  pos <- suppressWarnings(as.integer(unlist(
    pluck_at(aa, "pos"),
    use.names = FALSE
  )))
  pos <- pos[!is.na(pos)]
  if (is_blank(ref) || is_blank(alt) || length(pos) == 0) {
    return(NA_character_)
  }
  paste0(ref, pos[[1]], if (identical(alt, "X")) "*" else alt)
}

# Primary clinical significance (label + severity rank) from the "; "-joined
# significance string, so suggestions can lead with the most severe call.
.mv_sig_primary <- function(sig) {
  terms <- tolower(trimws(strsplit(sig %||% "", ";", fixed = TRUE)[[1]]))
  if ("pathogenic" %in% terms) {
    return(list(label = "Pathogenic", rank = 1L))
  }
  if ("likely pathogenic" %in% terms) {
    return(list(label = "Likely pathogenic", rank = 2L))
  }
  first <- terms[nzchar(terms)]
  list(
    label = if (length(first) == 0) {
      "ClinVar"
    } else {
      tools::toTitleCase(first[[1]])
    },
    rank = 3L
  )
}

# Pure parser: turn gene-scoped hits into a ranked, de-duplicated variant table
# (most severe first, then highest CADD). One row per allele.
myvariant_parse_gene_variants <- function(hits) {
  empty <- list(ok = FALSE, error = "No notable variants found for this gene.")
  if (is.null(hits) || length(hits) == 0) {
    return(empty)
  }
  rows <- lapply(hits, function(h) {
    id <- pluck_at(h, "_id")
    rsid <- mygene_first(pluck_at(h, "dbsnp", "rsid"))
    if (is_blank(id) || is_blank(rsid)) {
      return(NULL)
    }
    prim <- .mv_sig_primary(myvariant_clinvar_sig(h))
    data.frame(
      id = as.character(id),
      rsid = tolower(rsid),
      label = .mv_aa_label(pluck_at(h, "dbnsfp", "aa")),
      significance = prim$label,
      rank = prim$rank,
      cadd = .mv_cadd(h),
      stringsAsFactors = FALSE
    )
  })
  rows <- do.call(rbind, rows)
  if (is.null(rows) || nrow(rows) == 0) {
    return(empty)
  }
  rows <- rows[order(rows$rank, -ifelse(is.na(rows$cadd), -Inf, rows$cadd)), ]
  rows <- rows[!duplicated(rows$id), ]
  rows$rank <- NULL
  rownames(rows) <- NULL
  list(ok = TRUE, variants = rows)
}

# Named character vector for a selectizeInput: value = the allele's hg38 HGVS
# id, name = display label like "V600E, rs113488022 (Pathogenic)". The value is
# the HGVS id, not the rsID, because one rsID can be several alleles and the
# choice has to say which one. Falls back to the rsID in the label when there
# is no amino-acid change (e.g. splice/frameshift variants).
myvariant_variant_choices <- function(parsed, max_n = 100) {
  if (is.null(parsed) || !isTRUE(parsed$ok)) {
    return(character())
  }
  v <- parsed$variants
  if (nrow(v) > max_n) {
    v <- v[seq_len(max_n), ]
  }
  disp <- ifelse(
    is.na(v$label),
    sprintf("%s (%s)", v$rsid, v$significance),
    sprintf("%s, %s (%s)", v$label, v$rsid, v$significance)
  )
  stats::setNames(v$id, disp)
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
