# Get all registered fit profiles

Returns the list of all experiment-type profiles used by
[`apply_fit`](https://gdrplatform.github.io/gDRcore/reference/apply_fit.md)
and
[`apply_fits`](https://gdrplatform.github.io/gDRcore/reference/apply_fits.md)
to resolve default slicing configuration.

## Usage

``` r
get_fit_profiles()
```

## Value

Named list of profile definitions. Each element contains `slicing_cols`,
`slicing_values`, `input_assay`, and `description`.

## Examples

``` r
get_fit_profiles()
#> $`single-agent`
#> $`single-agent`$slicing_cols
#> [1] "normalization_type"
#> 
#> $`single-agent`$slicing_values
#> [1] "GR" "RV"
#> 
#> $`single-agent`$input_assay
#> [1] "Averaged"
#> 
#> $`single-agent`$nested_cols
#> [1] "concentration"
#> 
#> $`single-agent`$nested_cols_note
#> [1] "identifier keys resolved at runtime via gDRutils::get_env_identifiers() — actual column names may differ in custom environments"
#> 
#> $`single-agent`$description
#> [1] "1D dose-response curve; fit_fn receives one data.table per normalization_type slice, rows indexed by concentration"
#> 
#> 
#> $`time-course`
#> $`time-course`$slicing_cols
#> [1] "normalization_type"
#> 
#> $`time-course`$slicing_values
#> NULL
#> 
#> $`time-course`$input_assay
#> [1] "LogFoldChange"
#> 
#> $`time-course`$nested_cols
#> [1] "concentration" "duration"     
#> 
#> $`time-course`$nested_cols_note
#> [1] "Stage 1 of fit_SE.timecourse(): compute_growth_rates() reads this assay — it is the default value of the lfc_assay argument — and applies rate_fn (by default the slope of lm(LogFoldChange ~ Duration)) per well per half-open time window [start, end). Duration is consumed here; it does not appear in Stage 2. Stage 1 aggregates and normalises across cells, so it does not go through apply_fit() and slicing_values is null: LogFoldChange carries no GR/RV split. identifier keys resolved at runtime via gDRutils::get_env_identifiers()"
#> 
#> $`time-course`$description
#> [1] "Time-course Stage 1: log2 fold-change time-series; compute_growth_rates() iterates rows indexed by concentration x duration to compute per-well growth rates, pluggable via rate_fn. Output is the 'time-course-metrics' GrowthRates assay."
#> 
#> 
#> $combination
#> $combination$slicing_cols
#> [1] "normalization_type"
#> 
#> $combination$slicing_values
#> [1] "GR" "RV"
#> 
#> $combination$input_assay
#> [1] "Averaged"
#> 
#> $combination$nested_cols
#> [1] "concentration"  "concentration2"
#> 
#> $combination$nested_cols_note
#> [1] "identifier keys resolved at runtime via gDRutils::get_env_identifiers() — actual column names may differ in custom environments"
#> 
#> $combination$description
#> [1] "2D concentration grid; fit_fn receives one data.table per normalization_type slice, rows indexed by concentration x concentration2"
#> 
#> 
#> $`time-course-metrics`
#> $`time-course-metrics`$slicing_cols
#> [1] "normalization_type"
#> 
#> $`time-course-metrics`$slicing_values
#> [1] "NGR"
#> 
#> $`time-course-metrics`$input_assay
#> [1] "GrowthRates"
#> 
#> $`time-course-metrics`$nested_cols
#> [1] "concentration"
#> 
#> $`time-course-metrics`$nested_cols_note
#> [1] "Stage 2 of fit_SE.timecourse(): apply_fit() reads GrowthRates (produced by Stage 1) and fits a logistic curve on NormalizedGrowthRate ~ log10(Concentration). normalization_type='NGR' (NormalizedGrowthRate) — distinct from single-agent 'GR' (GR value) and 'RV' (Relative Viability): time-course has no RV/GR split, only a single rate-based metric. BumpyMatrix row key = DrugName|partner drug slots|partner concentration slots|CellLineName — every slot present in the data is used (DrugName_2/Concentration_2 for doublets, DrugName_3/Concentration_3 for triplets), so SA and each combo arm (same DrugA at different partner concentrations) are fitted as separate curves. The partner concentrations are present in GrowthRates cells alongside Concentration. identifier keys resolved at runtime via gDRutils::get_env_identifiers()"
#> 
#> $`time-course-metrics`$description
#> [1] "Time-course Stage 2: dose-response logistic fit on NormalizedGrowthRate ~ log10(Concentration) per (Drug, CellLine, [DrugName_2, Concentration_2], period) combination. fit_fn receives one data.table per normalization_type='NGR' slice, rows indexed by concentration (x axis) with Concentration_2 as additional context."
#> 
#> 
```
