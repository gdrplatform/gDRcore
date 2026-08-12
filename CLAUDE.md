# CLAUDE.md — gDRcore

@../CLAUDE.md

## Commands

### Pre-push check (pkgdown + R CMD check, no tests)

Run this before pushing when you've added new exported functions or changed
`_pkgdown.yml` / NAMESPACE. Takes ~2 min. Catches pkgdown index errors and
R CMD check WARNINGs/NOTEs early.

```bash
podman exec rplatform_rstudio Rscript -e "
gDRstyle::checkPackage(
  'gDRcore',
  '/home/rstudio/projects/open-sourced',
  'gDRcore',
  fail_on = 'note',
  skip_tests = TRUE,
  skip_lint = TRUE,
  build_vignettes = FALSE
)
"
```

**Why `build_vignettes = FALSE`:** `gDRcore.Rmd` calls `gDRinternal::connect()` →
`googlesheets_connect()` which requires Google credentials not available in the
local container. This is pre-existing; CI on resdev has access to gstore.
`custom-fitting.Rmd` and other vignettes build fine locally.

### Running tests

```bash
podman exec rplatform_rstudio Rscript -e "
devtools::load_all('/home/rstudio/projects/open-sourced/gDRcore', quiet = TRUE)
testthat::test_file('/home/rstudio/projects/open-sourced/gDRcore/tests/testthat/test-fit_utils.R')
"
```

### Building vignettes (excluding gDRcore.Rmd)

```bash
podman exec rplatform_rstudio Rscript -e "
devtools::load_all('/home/rstudio/projects/open-sourced/gDRcore', quiet = TRUE)
# Build specific vignette only:
rmarkdown::render('/home/rstudio/projects/open-sourced/gDRcore/vignettes/custom-fitting.Rmd')
"
```

## Known Issues

- `gDRcore.Rmd` requires `GDR_GOOGLE_TOKEN_DIR` — skip locally with
  `build_vignettes = FALSE` in `checkPackage()`.
- `_pkgdown.yml` must list every exported function. After adding new exports,
  run `devtools::document()` then add the functions to `_pkgdown.yml` before
  pushing.
