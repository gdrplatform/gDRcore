#' @keywords internal
fit_combo_cotreatments <- function(measured, series_id, cotrt_id, normalization_type) {
  series_concs <- setdiff(unique(measured[, series_id]), 0)
  cotrt_concs <- unique(measured[, cotrt_id])   # keep the single agent for the series_id

  if (all(measured$normalization_type != normalization_type)) {
    measured$x <- NA
    measured$x_std <- NA
    measured$normalization_type <- normalization_type
  } else {
    measured <- measured[measured$normalization_type == normalization_type, ]
  }
  
  # Single agent fit for the cotrt (to be used as the reference) --> ids are flipped
  sa_fit <- fit_cotreatment_series(measured, series_id = cotrt_id, cotrt_id = series_id, cotrt_value = 0,
    normalization_type = normalization_type, e_0 = 1, GR_0 = 1)
    
  # Fit cotreatments in matrix.
  cotrt_fittings <- vector("list", length(cotrt_concs))
  for (i in seq_along(cotrt_concs)) {
    conc <- cotrt_concs[i]
    sa <- gDRutils::predict_efficacy_from_conc(conc, sa_fit$x_inf, sa_fit$x_0, sa_fit$ec50, sa_fit$h)

    if (is.na(sa) &
        any(conc == measured[, cotrt_id] & measured[, series_id] == 0)) {
      
      # if the fit or the prediction fails, tries to get the reference value from the actual data
      sa <- measured[measured[, cotrt_id] == conc & measured[, series_id] == 0, normalization_type]
    } # else x_0 will be NA (thus a free parameter)

    cotrt_fittings[[i]] <- fit_cotreatment_series(measured, series_id = series_id, cotrt_id = cotrt_id,
      cotrt_value = conc, normalization_type = normalization_type, e_0 = sa, GR_0 = sa)    
  }
  do.call("rbind", cotrt_fittings)
}


#' @keywords internal
fit_cotreatment_series <- function(measured, series_id, cotrt_id, cotrt_value, normalization_type, e_0, GR_0) {
  checkmate::assert_numeric(cotrt_value)
  df_ <- measured[measured[[cotrt_id]] == cotrt_value & measured[[series_id]] > 0, , drop = FALSE]
  cotrt_fit <- gDRutils::fit_curves(
    df_ = df_,
    series_identifiers = series_id,
    e_0 = e_0,
    GR_0 = GR_0,
    cap = 0.1,
    normalization_type = normalization_type
  )

  cotrt_fit$cotrt_value <- cotrt_value
  cotrt_fit
}
