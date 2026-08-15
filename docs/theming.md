# Theming

This template themes the app from one file: [`_brand.yml`](../_brand.yml). It
uses the [brand.yml](https://posit-dev.github.io/brand-yml/) standard, applied
through [bslib](https://rstudio.github.io/bslib/).

## How it works

- `_brand.yml` defines the brand: the color palette, semantic colors (primary,
  secondary, foreground, background), typography (fonts, sizes, weights), and
  an optional logo.
- [app_ui.R](../app_ui.R) calls `bslib::bs_theme(brand = TRUE)`. This finds
  `_brand.yml` at the app root, and applies it to the whole UI.
- `brand = TRUE` makes the file a requirement: the app fails to start without
  it. To make the file optional, use `bslib::bs_theme()` instead. It applies
  `_brand.yml` when the file exists, and does nothing when the file is
  missing.

## Customizing

Edit `_brand.yml`. For example, to change the primary color and the base
font:

```yaml
color:
  palette:
    blue: "#1d4ed8"
  primary: blue

typography:
  fonts:
    - family: Roboto
      source: google
      weight: [400, 600]
  base: Roboto
```

Restart the app to see the changes. The
[brand.yml specification](https://posit-dev.github.io/brand-yml/articles/brand-yml.html)
documents the full set of fields.

## Theming plots and tables

bslib themes the HTML and CSS of the UI. R draws plots separately. Install
[`thematic`](https://rstudio.github.io/thematic/) so base R, ggplot2, and
lattice graphics use the app's colors automatically:

```r
install.packages("thematic")
```

When `thematic` is installed, [global.R](../global.R) calls
`thematic::thematic_shiny(font = "auto")` on its own, so you do not need to
add anything else. To also render custom or Google fonts (such as Inter) in
plots, install [`showtext`](https://github.com/yixuan/showtext). Without
`showtext`, thematic applies the theme colors, but falls back to the default
graphics-device font.

## Notes

- Keep `_brand.yml` in version control, so the look of the app stays
  reproducible.
