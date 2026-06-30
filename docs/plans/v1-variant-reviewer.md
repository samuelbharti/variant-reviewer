# Variant Reviewer — v1 Design

## Context

Variant Reviewer is a lightweight gene/variant interpretation companion. For a
given gene (and optionally a variant) it answers: what does this gene do, where
is it expressed, what does it interact with, and what does the variant mean at
the protein level. v1 is a single dashboard page driven by public, no-auth APIs.

## Scope (v1)

- **Gene-first** search box (gene symbol required; optional variant field).
- **Single scrolling dashboard** of bslib cards.
- Summary cards: **gene** (MyGene), **variant** (MyVariant), **protein** (ProtVar).
- **GTEx** tissue expression (ggplot).
- **STRING** PPI table (reactable).
- **External links** card.
- Stack: httr2 (API calls), reactable (tables), jsonlite, ggplot2,
  shinycssloaders; renv for reproducibility.
- Later (seams left open): ClinVar, gnomAD, Ensembl, Open Targets.

## Architecture

Data flow: the search module emits a reactive `{gene, variant}` → `server.R`
resolves the gene once via MyGene (symbol → Ensembl / Entrez / UniProt) and
fans the resolved reactive out to each result module → each module calls a
pure-R API client and renders its card.

Separation of concerns:

- Pure-R API clients + parsers live in `R/` (no Shiny, unit-testable).
- Shiny modules in `modules/` only orchestrate reactivity and rendering.

### Files

API clients (`R/`, auto-sourced by `R/load_components.R`):

- `api_http.R` — `vr_api_get()`: one httr2 wrapper for URL/query, timeout,
  retries, and normalized `list(ok, status, data, error)` results. Also the
  shared `is_blank()` / `pluck_at()` helpers.
- `api_mygene.R` — `mygene_resolve()` + `mygene_parse_hit()`.
- `api_myvariant.R` — `myvariant_annotate()` + `myvariant_parse_hit()`.
- `api_protvar.R` — `protvar_annotate()` (per-position `/function` +
  `/population`), `protvar_parse_position()`, parsers for function text/variants.
- `api_gtex.R` — `gtex_median_expression()` (resolves a versioned GENCODE id via
  GTEx's own gene lookup, then queries with `datasetId=gtex_v8`) +
  `gtex_parse_rows()`.
- `api_string.R` — `string_interaction_partners()` + `string_parse_rows()`.
- `ui_helpers.R` — `vr_field()`, `vr_empty()`, `vr_error()`, `vr_num()`, `%||%`.

Modules (`modules/`, `*_ui()`/`*_server()`): `gene_search`, `gene_summary`,
`variant_summary`, `protein_summary`, `gtex_expression`, `string_ppi`,
`external_links`. Each result card uses `shinycssloaders::withSpinner()` and
shows a friendly empty/error placeholder.

App shell: `userInterface/dashboard_ui.R` composes the cards; `ui.R` is a
`page_fluid` with `bs_theme(brand = TRUE)`; `server.R` wires search → resolved →
modules.

### API quirks worth remembering

- GTEx needs **both** `datasetId=gtex_v8` and a **versioned** GENCODE id; the
  versioned id comes from GTEx's `/reference/gene` lookup, not MyGene.
- ProtVar uses per-position paths (`/function/{acc}/{pos}`,
  `/population/{acc}/{pos}`), not a single mappings endpoint.
- STRING returns `text/json`, so JSON parsing must skip the content-type check.
- `brand.yml` is a phantom renv dependency (referenced only via
  `bs_theme(brand = TRUE)`); it is recorded into `renv.lock` explicitly.

## Testing

- Offline parser tests against recorded fixtures in `tests/testthat/fixtures/`.
- `shiny::testServer()` reactive tests for the search module.
- `shinytest2` end-to-end: app-launch smoke test (always) + live-API gene
  search (skipped on CI).
- `lintr` + `air` clean.

## Open seams

ClinVar / gnomAD / Ensembl / Open Targets each become an additional
`R/api_*.R` + `modules/*_mod.R` pair added to the dashboard grid — no
architectural change required.
