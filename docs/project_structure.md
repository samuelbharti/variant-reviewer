# Project Structure

```txt
.
├── app.R
├── app_ui.R
├── app_server.R
├── global.R
├── R/
├── modules/
├── userInterface/
├── tests/
├── data/
├── dev/
├── docs/
└── www/
```

## Notes

- `app.R` is the entry point. It sources `global.R`, then loads the UI from
  `app_ui.R` and the server function from `app_server.R`.
- `global.R` loads dependencies and sources components.
- `R/load_components.R` sources modules and page UI files automatically.
- `www/` stores static assets, such as CSS, JavaScript, and images.
