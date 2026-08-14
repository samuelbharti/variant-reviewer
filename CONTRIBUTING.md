# Contributing Guidelines

## Branching

- Create feature branches from `dev`.
- Open pull requests into `dev`, unless someone tells you to use a different
  branch.

## Local Setup

1. Restore dependencies with `renv::restore()`.
2. Run the app locally with `shiny::runApp()`.
3. (Recommended) Install the git pre-commit hooks:

   ```bash
   pip install pre-commit
   pre-commit install
   ```

## Code Style

- Keep page UI definitions in `userInterface/`.
- Keep reusable UI/server logic in `modules/`.
- Keep utility functions in `R/`. `R/load_components.R` sources these files
  automatically.
- Format R code with [air](https://posit-dev.github.io/air/): `air format .`
- Lint with `lintr::lint_dir(".")` (config in `.lintr`).

## Testing

- Tests live in `tests/testthat/`. Run them with:

  ```r
  shiny::runTests(".")
  ```

- Add tests for each kind of change:
  - Unit tests for `R/` helpers.
  - `shiny::testServer()` tests for module reactivity.
  - `shinytest2` tests for end-to-end behavior.

## Template Mechanism

This repo is a project template. Two files turn it into a real project:

- `template.yml`: the manifest. It lists variables, find-and-replace targets,
  files to strip, and files to reset. This file is the single source of truth
  for template decisions.
- `dev/use_template.R`: a generic engine that applies the manifest. A person
  can call it directly (`Rscript dev/use_template.R --project_name="..."`), or
  an initializer package can call it in code (`use_template(values = ...)`).

When you change a template decision, change `template.yml`, not the engine.
Template-only identity (citation files, the Zenodo badge, and the "How to
cite" section) lives inside `<!-- template:strip:start --> ... :end -->`
markers, so init removes it. If you mirror this template into a package's
`inst/` folder, treat that copy as a generated mirror of a tagged release. Do
not hand-edit it.

## Continuous Integration

Every push and pull request runs the `CI` workflow
(`.github/workflows/ci.yaml`):

- lint
- formatting check
- the test suite
- Markdown linting

Before you open a PR, make sure that these checks pass locally.

## Pull Request Checklist

- [ ] App runs locally (`shiny::runApp()`).
- [ ] Code is formatted (`air format .`) and lints clean (`lintr::lint_dir(".")`).
- [ ] Tests pass (`shiny::runTests(".")`).
- [ ] New or changed code follows the project structure.
- [ ] If behavior changed, README and docs are updated.
