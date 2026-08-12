# CLAUDE.md — gDRcore

@../CLAUDE.md

## Commands

### Pre-push check (pkgdown + R CMD check, no tests) — always run before pushing

Takes ~2 min. Catches pkgdown index errors and R CMD check WARNINGs/NOTEs.

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
`googlesheets_connect()` which requires Google credentials not available locally.
Pre-existing; CI on resdev has gstore access. All other vignettes build fine locally.

---

### When vignettes change — run BEFORE pushing

Three additional checks to catch issues that `checkPackage(build_vignettes=FALSE)`
misses. All three together take ~3 min and replace 45+ min CI wait.

**Step 1 — Rmd dependency linter** (~5 sec):
Finds unqualified function calls in R chunks that may fail at `R CMD check --run-examples`.
False positives possible for functions exported by the package itself (not yet installed).

```bash
podman exec rplatform_rstudio Rscript \
  /usr/local/lib/R/site-library/gDRstyle/scripts/lint_rmd_deps.R \
  /home/rstudio/projects/open-sourced/gDRcore/vignettes/
```

**Step 2 — BiocCheck vignette directory** (~30 sec):
Catches YAML front matter issues (e.g. too many `---` delimiters) and other
Bioconductor vignette policy checks. Run on the package directory, not the file.

```bash
podman exec rplatform_rstudio Rscript -e "
BiocCheck::BiocCheck(
  '/home/rstudio/projects/open-sourced/gDRcore',
  'no-check-R-ver' = TRUE, 'no-check-bioc-help' = TRUE,
  'no-check-unit-tests' = TRUE, 'no-check-coding-practices' = TRUE,
  'no-check-function-len' = TRUE, 'no-check-man-pages' = TRUE
)
" 2>&1 | grep -E 'ERROR|WARNING|vignette|YAML|---'
```

**Step 3 — render the changed vignette(s)** (~2 min):
Confirms the vignette actually runs end-to-end.

```bash
podman exec rplatform_rstudio Rscript -e "
devtools::load_all('/home/rstudio/projects/open-sourced/gDRcore', quiet = TRUE)
rmarkdown::render(
  '/home/rstudio/projects/open-sourced/gDRcore/vignettes/custom-fitting.Rmd',
  quiet = TRUE
)
cat('OK\n')
"
```

---

### Running tests

```bash
podman exec rplatform_rstudio Rscript -e "
devtools::load_all('/home/rstudio/projects/open-sourced/gDRcore', quiet = TRUE)
testthat::test_file('/home/rstudio/projects/open-sourced/gDRcore/tests/testthat/test-fit_utils.R')
"
```

---

### After adding new exported functions — always run before pushing

New exported functions must appear in both NAMESPACE and `_pkgdown.yml`, and all
their `@param` entries must be documented. Missing entries cause R CMD check ERRORS
and pkgdown `build_reference_index()` failures.

```bash
# 1. Regenerate NAMESPACE + Rd files
podman exec rplatform_rstudio Rscript -e "
devtools::document('/home/rstudio/projects/open-sourced/gDRcore', quiet = TRUE)
"

# 2. Verify pkgdown reference index (fast, catches missing _pkgdown.yml entries)
podman exec rplatform_rstudio Rscript -e "
pkgdown::build_reference_index('/home/rstudio/projects/open-sourced/gDRcore')
"
```

Then add the new functions to `_pkgdown.yml` reference section manually.

---

## Known issues / gotchas

- **`gDRcore.Rmd`** requires `GDR_GOOGLE_TOKEN_DIR` — skip with `build_vignettes = FALSE`.
- **BiocCheck `---` rule**: standalone `---` lines in vignette body are misread
  as YAML delimiters. Use `***` for horizontal rules instead of `---`.
- **`_pkgdown.yml`** must list every exported function — after adding new exports,
  run `devtools::document()` then add them to `_pkgdown.yml`.
- **`data.table::` prefix**: R CMD check runs `@examples` in a clean namespace.
  Always use `data.table::as.data.table()`, never bare `as.data.table()`.
- **gDRstyle lintRmdDeps** may report false positives for functions from the
  package itself if the package isn't installed in the container yet (e.g. new
  functions on a feature branch). Verify manually if unsure.
