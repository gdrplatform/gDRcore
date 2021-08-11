#' @keywords internal
fit_combo_cotreatments <- function(measured, series_id, cotrt_id, normalization_type) {
  series_concs <- setdiff(unique(measured[, series_id]), 0)
  cotrt_concs <- setdiff(unique(measured[, cotrt_id]), 0)

  # Single agent fits.
  sa_fit <- fit_cotreatment_series(measured, series_id = series_id, cotrt_id = cotrt_id, cotrt_value = 0,
    normalization_type = normalization_type, e_0 = 1, GR_0 = 1)

  # Fit cotreatments in matrix.
  cotrt_fittings <- vector("list", length(all_conc))
  for (i in seq_along(cotrt_concs)) {
    # TODO: Switch to ec50 once c50 -> ec50
    conc <- cotrt_concs[i]
    sa <- gDRutils::logistic_4parameters(conc, sa_fit$x_inf, sa_fit$x_0, sa_fit$c50, sa_fit$h)
    cotrt_fittings[[i]] <- fit_cotreatment_series(measured, series_id = series_id, cotrt_id = cotrt_id,
      cotrt_value = conc, normalization_type = normalization_type, e_0 = sa, GR_0 = sa)
  }

  do.call("rbind", cotrt_fittings)
}


#' @keywords internal
fit_cotreatment_series <- function(measured, series_id, cotrt_id, cotrt_value, normalization_type, e_0, GR_0) {
  assert_numeric(cotrt_value)
  df_ <- measured[measured[[cotrt_id]] == cotrt_value & measured[[series_id]] > 0, , drop = FALSE]
  cotrt_fit <- gDRutils::fit_curves(
    df_ = df_,
    series_identifiers = series_id,
    e_0 = e_0,
    GR_0 = GR_0,
    force_fit = TRUE,
    cap = 0.2,
    normalization_type = normalization_type
  )

  nested_df <- DataFrame()
  nested_df[series_id] <- df_[[series_id]]
  nested_df[cotrt_id] <- cotrt_value

  out <- DataFrame(DataFrame(nested_df), cotrt_fit)
  out
}
