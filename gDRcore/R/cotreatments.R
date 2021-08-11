#' @keywords internal
fit_combo_cotreatments <- function(measured, nested_identifiers, normalization_type) {       
  id <- nested_identifiers[1] # (i.e. conc1)
  id2 <- nested_identifiers[2] # (i.e. conc2)

  all_conc1 <- setdiff(unique(measured[, id]), 0)
  all_conc2 <- setdiff(unique(measured[, id2]), 0)

  # Single agent fits.
  sa_fit1 <- fit_cotreatment_series(measured, series_id = id, cotrt_id = id2, cotrt_value = 0,
    normalization_type = normalization_type, e_0 = 1, GR_0 = 1)
  sa_fit2 <- fit_cotreatment_series(measured, series_id = id2, cotrt_id = id, cotrt_value = 0,
    normalization_type = normalization_type, e_0 = 1, GR_0 = 1)

  # Fit cotreatments in matrix.
  cotrt_fittings <- vector("list", sum(length(all_conc1), length(all_conc2)))
  n <- 1
  for (conc1 in all_conc1) {
    # TODO: Switch to ec50 once c50 -> ec50
    single_agent2 <- gDRutils::logistic_4parameters(conc1, sa_fit1$x_inf, sa_fit1$x_0, sa_fit1$c50, sa_fit1$h)
    cotrt_fittings[[n]] <- fit_cotreatment_series(measured, series_id = id2, cotrt_id = id,
      cotrt_value = conc1, normalization_type = normalization_type, e_0 = single_agent2, GR_0 = single_agent2)
    n <- n + 1
  }
  for (conc2 in all_conc2) {
    single_agent1 <- gDRutils::logistic_4parameters(conc2, sa_fit2$x_inf, sa_fit2$x_0, sa_fit2$c50, sa_fit2$h)
    cotrt_fittings[[n]] <- fit_cotreatment_series(measured, series_id = id, cotrt_id = id2,
      cotrt_value = conc2, normalization_type = normalization_type, e_0 = single_agent1, GR_0 = single_agent1)
    n <- n + 1
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

  out <- DataFrame(IRanges::DataFrameList(nested_df), cotrt_fit)
  out
}
