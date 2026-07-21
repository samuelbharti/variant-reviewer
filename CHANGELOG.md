# Changelog

All notable changes to this project should be documented in this file.

## [Unreleased]

- The assistant now connects on its own when a provider key is set in the
  environment: the matching provider starts selected, its default model is
  preselected, and the client is built on load, so the user can chat without
  opening the drawer or clicking Connect. Pasting a key or switching provider
  still connects manually.
- A search now accepts a gene, a variant, or both. A variant entered on its own
  resolves its own gene (via MyVariant), so the gene-level cards still fill.
- When a gene and a variant are both entered but name different genes, the
  search is blocked with a clear message instead of showing a mixed result. The
  old behaviour only warned.
- Widened the assistant column on wide screens for a bit more reading room.
- The **Gene model** card now starts un-ticked. rest.ensembl.org is slow and
  frequently times out, and hidden cards make no API call, so it no longer
  slows a search by default; users can enable it from the **Cards** popover. The
  demo walkthrough skips it to match.

- Fixed the 3D structure only appearing on the first search. The viewer was
  rebuilt by a `renderUI` every time, and an htmlwidget re-created that way comes
  back empty with the server's value already sent and nowhere to land. The
  element now lives in the UI permanently and only its contents change.
- Unticking a card in the **Cards** popover now skips its API call instead of
  just hiding the result. Hidden cards were already suspended by Shiny; the only
  thing still fetching them was the loop that mirrors results for the assistant.
  The progress total follows suit, so a search with cards hidden reports fewer
  sources.
- The 3D structure card now reports a coordinate download that fails, instead of
  silently rendering an empty viewer.
- Bounded the API cache and made it tunable. It now has entry and byte ceilings
  with least-recently-used eviction, so a long-running server cannot grow without
  limit, and `VR_CACHE_TTL` / `VR_CACHE_MAX_SIZE` / `VR_CACHE_MAX_N` override the
  defaults. Added `vr_cache_stats()` and `vr_cache_clear()`.
- Widened the assistant column to a fixed 460px track on wide screens, and
  render the example prompts as a bullet list instead of shinychat's card grid.
  That grid sizes its tracks from a minimum card width, which in a column this
  narrow came out wider than the column and left a horizontal scrollbar; a plain
  list has no minimum, so it wraps. The prompts stay clickable.
- Added a search progress popup (bottom right) that counts the data sources as
  they load, such as "Loading sources (7 of 16)" plus the one being fetched. The total
  reflects the search: a gene-only search skips the seven variant-level sources
  and counts 9. Fetching is synchronous, so the R process cannot flush outputs
  mid-search and the cards arrive in clumps; `shiny::Progress` is used because it
  writes straight to the websocket instead of waiting for the flush cycle. The
  sources are driven from one high-priority observer so it is the first consumer
  to touch each reactive, which is what keeps the count honest.
- Stripped the inline UniProt citations from the ProtVar **Function** text. The
  TP53 entry carried more `PubMed:` ids than prose; non-citation notes such as
  "(By similarity)" are kept.
- Simplified the assistant's **Model & key** drawer: pasting an API key now
  loads that key's models on its own, so the "List models for this key" button
  is gone. The "Forget key" button is gone too, since switching provider already
  clears the client and starts a fresh conversation. Gemini now starts on
  `gemini-flash-lite-latest`, and the provider default stays selectable even
  when the live list comes back without it.
- Limited the assistant to the search box. `set_selection` now fills the gene /
  variant inputs and runs the same submit as a **Review** click, instead of
  writing the query behind the search module's back. Reading cards is
  unchanged; the assistant cannot write card state, and its request gets the
  same validation and visible search box as a typed one.
- Fixed the assistant column changing width as the conversation grew. It sits in
  a bslib fill container (a column flexbox), where `align-self` controls the
  horizontal axis, so `align-self: start` shrink-wrapped it to its content, so
  the column resized whenever the greeting's suggestion cards gave way to a reply.
- Added four visualization cards: a **Variant landscape** protein "lollipop"
  (every ClinVar variant placed at its residue, coloured by significance, over
  the UniProt domain track, with the queried variant marked); a **Conservation**
  card (phyloP, phastCons, GERP++, SiPhy ranks from dbNSFP); an **Ancestry
  frequency** card (gnomAD allele frequency by genetic-ancestry group, reusing
  the shared gnomAD result); and a **Gene model** card (the canonical
  transcript's exons drawn 5'->3', with the variant's exon highlighted; the
  position is taken from gnomAD, and exons are numbered strand-aware). All four
  are readable by the assistant via `read_card`.
- Turned the variant box into a typeahead: once a gene is entered, it suggests
  that gene's known pathogenic / likely-pathogenic variants (from ClinVar via
  MyVariant), labelled by amino-acid change and rsID (e.g. `V600E, rs113488022
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
  provider and its Model & key control says so, inviting the user to just pick
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
