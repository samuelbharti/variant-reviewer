# Installation

## Recommended: renv

1. If `renv` is not installed, install it.
2. Run `renv::restore()` in the project root.
3. Start the app with `shiny::runApp()`.

If you use Docker, keep `renv.lock` and the `renv/` folder in the project
root. The image then restores the project library from the lockfile.

### Quick-start helper

This template includes a helper script that sets up `renv` for a new project.
Run:

```sh
Rscript dev/init-renv.R
```

The script installs a small set of recommended packages, then creates
`renv.lock`. Review the lockfile before you commit it.

## Manual setup

Install the packages listed in the README. Then run `shiny::runApp()`.

## Docker

Build:

```bash
docker build -t my-shiny-app .
```

Run:

```bash
docker run --rm -p 3838:3838 my-shiny-app
```

The Dockerfile restores packages from `renv.lock`. It does not install
packages one by one.
