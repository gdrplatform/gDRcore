# Persist a data.table into a BumpyMatrix assay

Performs an idempotent upsert (merge mode) keyed by `upsert_key_cols`,
or a full overwrite (replace mode).

## Usage

``` r
.persist_assay(se, new_dt, merge, assay_name, row, col, upsert_key_cols)
```

## Arguments

- se:

  SummarizedExperiment

- new_dt:

  data.table including "row" and "column" index fields

- merge:

  `"merge"` or `"replace"`

- assay_name:

  assay name to write into

- row:

  row index column name in new_dt

- col:

  column index column name in new_dt

- upsert_key_cols:

  character vector of column(s) forming the upsert key

## Value

updated SummarizedExperiment
