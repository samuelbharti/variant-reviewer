# Changelog

All notable changes to this project should be documented in this file.

## [Unreleased]

- Turned the variant box into a typeahead: once a gene is entered, it suggests
  that gene's known pathogenic / likely-pathogenic variants (from ClinVar via
  MyVariant), labelled by amino-acid change and rsID (e.g. `V600E — rs113488022
  (Pathogenic)`), filterable by typing the protein change. Any rsID/HGVS can
  still be entered freely. Suggestions are prefetched when the gene field
  settles (debounced) and gated by the same format check.
- Validate the gene symbol and variant with the
  [biobouncer](https://github.com/samuelbharti/biobouncer) package (offline
  `pattern` mode) before any API is queried, rejecting malformed input up front
  with an inline message. The same gate covers the assistant's `set_selection`.
  Documented on the About page ("Input validation") and in the README.
- Added the external-link icon to the buttons in the **External resources** card,
  matching the "view on source" links in the card headers.
- Added a "view on source" external link to every result card's header (NCBI
  Gene, dbSNP, UniProt, AlphaFold, gnomAD, ClinVar, Ensembl, GTEx, STRING, Open
  Targets), linking to the authoritative page for the loaded gene/variant.
- Tinted card headers with a soft wash of the brand primary (theme-derived via
  `color-mix`, with a hex fallback).
- Added a **3D structure** card: the AlphaFold-predicted structure for the
  gene's protein (via the AlphaFold DB API and the `r3dmol` viewer) with the
  variant residue highlighted. The coordinate file is downloaded only when the
  viewer renders; the assistant can read it via `read_card`.
- Added a **Protein domains & features** card (UniProt domains/regions/sites via
  the EBI Proteins API) that flags which feature the variant residue falls in;
  readable by the assistant (`read_card`) and toggleable in the Cards picker.
- Added an **In-silico predictions** card (REVEL, AlphaMissense, CADD,
  PolyPhen-2, SIFT, MetaLR/SVM from dbNSFP via MyVariant) and a **Gene
  constraint** card (pLI, LOEUF, Z-scores from gnomAD). Both are readable by the
  assistant (`read_card`) and toggleable in the Cards picker.
- Added a scope note under the search (human, GRCh38/hg38, Ensembl transcripts)
  and a non-blocking warning when the entered variant belongs to a different
  gene than the one typed.
- Made the chat Model & key controls fill the drawer width and widened the
  drawer by 50px.
- When an API key is found in the environment, the chat now defaults to that
  provider and its Model & key control says so — inviting the user to just pick
  a model and Connect (no pasting), instead of prompting for a key.
- Gave the assistant app-scoped tools: `get_current_selection` and `read_card`
  (read the gene/variant and any card's data), and `set_selection` (load a
  gene/variant into the dashboard). The assistant works through the app's own
  data and lookups rather than searching externally. Result modules now return
  their data reactive so the parent can surface it.
- Added a footer (app version, author, source/license links) shown on every
  page.
- Added a full-screen expand button to every result card.
- Added a "Cards" picker (popover) to show/hide individual result cards.
- Themed reactable tables to match the app theme via Bootstrap CSS variables.
- Enlarged the GTEx expression plot's axis labels and row spacing for
  legibility.
- Switched to an organic (earthy green/clay/stone) color palette across the
  bslib theme (`_brand.yml`) and plots (`vr_colors`).
- Defined and documented the scope of the app and of the assistant (About page
  "Scope" card and README); tightened the assistant's system prompt with
  explicit in/out-of-scope guardrails.
- Added a one-click "Load an example" (BRAF V600E) to the search box that
  populates the inputs for the user to review and submit.
- Laid the Home page out as results (left) beside a pinned chat column (right);
  moved the chat's provider/key/model controls into an offcanvas drawer with an
  always-visible connection badge in the chat header.
- Moved to a two-page navbar layout with a new About page (overview, scope, the
  annotations available and their sources, data-source/privacy notes) and added
  horizontal padding to the dashboard content.
- Added an optional bring-your-own-key (BYOK) AI chat assistant
  (`modules/byok_chat_mod.R`) supporting Google Gemini, OpenAI, and Anthropic.
  It is grounded in the current search via a `current_review_context` tool and
  degrades to a setup panel when `ellmer`/`shinychat` are absent.

## [2.2.0] - 2026-06-29

- Added GitHub Actions CI (lint, air format check, tests, markdownlint) and
  pre-commit hooks.
- Added air formatter and lintr configuration.
- Added a test suite (testthat unit tests, `testServer`, and a shinytest2
  smoke test).
- Added brand.yml theming applied through bslib, with optional `thematic`
  plot/table theming.
- Added a template manifest (`template.yml`) and `dev/use_template.R`
  scaffolding engine; published the repo as a GitHub template.
- Added issue and pull request templates and a Contributor Covenant code of
  conduct.
- Fixed the app entry point and ensured `R/` utilities are sourced at startup.

## [0.1.0] - 2026-05-01

- Created base Shiny template structure.

## [2.0.0] - 2026-05-02

- Bumped template version to v2.0 and updated metadata (CITATION, README).
