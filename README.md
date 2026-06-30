# Variant Reviewer

A lightweight gene and variant interpretation companion, built with Shiny. Enter
a gene symbol (and, optionally, a variant) to pull together — on a single
dashboard — what the gene does, where it is expressed, what it interacts with,
and protein-level context for the variant.

The app is a thin, reactive front end over several public bioinformatics APIs:

| Card / section            | Source   |
| ------------------------- | -------- |
| Gene summary              | MyGene   |
| Variant annotation        | MyVariant |
| Protein context           | ProtVar (EBI) |
| Tissue expression         | GTEx     |
| Protein interactions      | STRING   |
| External resource links   | derived  |

All sources are public and require no API key.

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
  "httr2", "reactable", "jsonlite", "shinycssloaders"
))
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
├── ui.R                    # Single-page dashboard layout
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
├── userInterface/          # Page-level layout (dashboard_ui.R)
├── tests/                  # testthat unit/reactive tests + shinytest2 e2e
└── docs/                   # Project documentation, including docs/plans/
```

### Architecture

The search module emits a reactive query `{gene, variant}`. `server.R` resolves
the gene **once** via MyGene (symbol → Ensembl / Entrez / UniProt) and shares
that with every gene-level module, so each API is queried only when needed. API
clients are pure R (in `R/`, individually testable); modules only orchestrate
reactivity and rendering. Adding a new source (e.g. ClinVar, gnomAD,
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

## Theming

Branding lives in [`_brand.yml`](_brand.yml) — colors, fonts, and logo in one
place, applied by bslib via `bs_theme(brand = TRUE)` in [ui.R](ui.R).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).
