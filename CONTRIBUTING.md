# Contributing

Thanks for looking. This is a solo project, so please open an issue before you
start on anything large. That way I can tell you early whether I want it, and
you do not waste the work. Small fixes are welcome as a pull request straight
away.

Please also read the [Code of Conduct](CODE_OF_CONDUCT.md).

## Setup

```r
renv::restore()   # install the pinned dependencies
shiny::runApp()   # run the app
```

The git hooks are optional but recommended:

```bash
pip install pre-commit
pre-commit install
```

## Where code goes

- `R/` for pure R helpers and API clients, with no Shiny.
  `R/load_components.R` sources them for you.
- `modules/` for Shiny modules, one card each.
- `userInterface/` for the page layouts.
- `tests/testthat/` for the tests.

To add a data source, add one `R/api_*.R` file and one `modules/*_mod.R` file.

## Before you open a pull request

Branch from `main` and open the pull request against `main`. Then check that
all of this passes:

```bash
air format .              # format
```

```r
lintr::lint_dir(".")      # lint, configured in .lintr
shiny::runTests(".")      # tests
shiny::runApp()           # the app still starts
```

Add a test for what you changed: a unit test for an `R/` helper, a
`shiny::testServer()` test for module reactivity, or a `shinytest2` test for
end-to-end behavior. If behavior changed, update the README too.

The `CI` workflow runs the same lint, format check, tests, and Markdown lint on
every push and pull request.
