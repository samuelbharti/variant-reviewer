# Variant Reviewer

A lightweight gene and variant interpretation companion, built with Shiny. Enter
a gene symbol (and, optionally, a variant) to pull together — on a single
dashboard — what the gene does, where it is expressed, what it interacts with,
and protein-level context for the variant.

The app is a thin, reactive front end over several public bioinformatics APIs:

| Card / section          | Source                            |
| ----------------------- | --------------------------------- |
| Gene summary            | MyGene                            |
| Variant annotation      | MyVariant                         |
| In-silico predictions   | dbNSFP (MyVariant)                |
| Protein context         | ProtVar (EBI)                     |
| Variant landscape       | ClinVar variants (gnomAD) + UniProt domains |
| Conservation            | dbNSFP (MyVariant)                |
| Ancestry frequency      | gnomAD                            |
| Protein domains         | UniProt (EBI Proteins)            |
| 3D structure            | AlphaFold DB                      |
| Clinical significance   | ClinVar (NCBI E-utilities)        |
| Population frequency    | gnomAD                            |
| Gene constraint         | gnomAD                            |
| Gene model              | Ensembl                           |
| Variant consequences    | Ensembl VEP                       |
| Tissue expression       | GTEx                              |
| Protein interactions    | STRING                            |
| Disease associations    | Open Targets                      |
| External resource links | derived                           |
| AI assistant (chat)     | BYOK: Gemini / OpenAI / Anthropic |

All data sources are public and require no API key. The optional AI assistant is
"bring your own key" — see [AI assistant](#ai-assistant) below.

## Scope

**The app** reviews one human gene and, optionally, one variant at a time,
aggregating public annotations read-only onto a single dashboard and linking out
to the primary sources. It is **not** for batch/VCF-scale analysis or variant
calling, non-human species, or clinical diagnosis, treatment, or
genetic-counseling advice.

**The assistant** helps interpret the gene or variant on screen and answers
general genomics questions within that scope; it declines unrelated requests and
is not a source of medical, diagnostic, or treatment advice.

## Requirements

- R (>= 4.3)
- Packages managed with `renv` (see `renv.lock`)

## Installation

```r
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
renv::restore()
```

Or install the core packages manually:

```r
install.packages(c(
  "shiny", "bslib", "brand.yml", "ggplot2",
  "httr2", "reactable", "jsonlite", "shinycssloaders",
  "r3dmol", # 3D structure card (AlphaFold viewer)
  # AI assistant (optional; the app degrades gracefully without them)
  "ellmer", "shinychat"
))

# biobouncer (input validation) is published on r-universe, not CRAN:
install.packages(
  "biobouncer",
  repos = c("https://samuelbharti.r-universe.dev", getOption("repos"))
)
```

## How to run

```r
shiny::runApp()
```

Or open the project in RStudio and click **Run App**.

## Build and run with Docker

```bash
docker build -t variant-reviewer .
docker run --rm -p 3838:3838 variant-reviewer
```

Then open [http://localhost:3838](http://localhost:3838).

## Project structure

```txt
.
├── _brand.yml              # Brand colors, fonts, logo (theming)
├── global.R                # Libraries and component loading
├── ui.R                    # Navbar layout: Home dashboard + About page
├── server.R                # Wires search -> resolved gene -> result modules
├── R/                      # Pure-R API clients + helpers (no Shiny)
│   ├── api_http.R          # Shared httr2 GET wrapper (timeouts, retries, errors)
│   ├── api_mygene.R        # Gene resolution
│   ├── api_myvariant.R     # Variant annotation
│   ├── api_protvar.R       # Protein-level context
│   ├── api_gtex.R          # Tissue expression
│   ├── api_string.R        # Interaction partners
│   └── ui_helpers.R        # Small presentation helpers
├── modules/                # Shiny modules (one card each)
├── userInterface/          # Page-level layout (dashboard_ui.R, about_ui.R)
├── tests/                  # testthat unit/reactive tests + shinytest2 e2e
└── docs/                   # Project documentation, including docs/plans/
```

### Architecture

The search module emits a reactive query `{gene, variant}`. As the gene field
settles, the variant box prefetches that gene's known pathogenic/likely-
pathogenic variants (ClinVar via MyVariant) as typeahead suggestions, while
still accepting any free-typed rsID/HGVS. Before any API is
queried, the gene and variant are validated with
[biobouncer](https://github.com/samuelbharti/biobouncer)'s offline `pattern`
mode (gene against the HGNC grammar, rsIDs against dbSNP), so malformed input is
rejected up front with an inline message rather than firing failing lookups —
the same gate covers the assistant's `set_selection`. `server.R` resolves
the gene **once** via MyGene (symbol → Ensembl / Entrez / UniProt) and shares
that with every gene-level module, so each API is queried only when needed. API
clients are pure R (in `R/`, individually testable); modules only orchestrate
reactivity and rendering. Successful HTTP responses are cached in-process for 30
minutes (at the `vr_api_get()`/`vr_api_post_json()` chokepoint), so repeated
searches are instant; failures are never cached. Adding a new source (e.g. ClinVar, gnomAD,
Open Targets) is a new `R/api_*.R` + `modules/*_mod.R` pair dropped into the
dashboard grid.

## Testing

```r
shiny::runTests(".")
```

- **Parser tests** (`tests/testthat/test-api-parsers.R`) run offline against
  recorded JSON fixtures in `tests/testthat/fixtures/`.
- **Reactive tests** (`test-modules.R`) use `shiny::testServer()`.
- **End-to-end** (`test-shinytest2.R`) launches the app in a headless browser;
  the live-API search test is skipped on CI for determinism (run locally with
  `NOT_CRAN=true`).

## AI assistant

The dashboard includes an optional chat assistant ([modules/byok_chat_mod.R](modules/byok_chat_mod.R))
for discussing the gene or variant you're reviewing. It is **bring your own key
(BYOK)**: open the chat's **Model & key** drawer (the gear button in the chat
header), pick a provider (Google Gemini, OpenAI, or Anthropic), and paste your
own API key, which is held only in your session's server memory and never
written to disk. Alternatively, set a server-side key via the
matching environment variable (`GEMINI_API_KEY` / `GOOGLE_API_KEY`,
`OPENAI_API_KEY`, or `ANTHROPIC_API_KEY`) and connect with the key field blank.

The assistant is grounded in the dashboard through app-scoped tools (wired in
[server.R](server.R), formatters in [R/chat_tools.R](R/chat_tools.R)):

- `get_current_selection` — the gene/variant currently loaded.
- `read_card` — the data shown in a specific card (gene, variant, predictions,
  protein, domains, structure, clinvar, gnomad, constraint, consequences,
  expression, interactions, diseases).
- `set_selection` — load a gene (and optional variant) into the dashboard,
  running the app's own lookups.

It works through the app's data and lookups rather than searching externally. It
requires the `ellmer` and `shinychat` packages; if they are absent the card
renders a short setup panel and the rest of the app loads normally.

## Theming

Branding lives in [`_brand.yml`](_brand.yml) — colors, fonts, and logo in one
place, applied by bslib via `bs_theme(brand = TRUE)` in [ui.R](ui.R).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).
