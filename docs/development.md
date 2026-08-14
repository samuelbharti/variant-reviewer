# Development Guide

## Structure

- Keep global setup in `global.R`.
- Keep tab and page layouts in `userInterface/`.
- Keep reusable module pairs in `modules/`.
- Keep utility helpers in `R/`.

## Workflow

1. Create a branch from `dev`.
2. Add UI or server changes to the matching folder.
3. Run the app locally.
4. Whenever a package version changes, restore dependencies with
   `renv::restore()`.

## Deployment

You can publish this app directly, or deploy it with Docker. You can add
CI/CD workflows later, for your own project. This guide does not include them.
