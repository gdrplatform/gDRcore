#' @keywords internal
fit_cotreatments <- function(measured, nested_identifiers, normalization_type) {       
  # Identify the combo concentrations.
  id <- nested_identifiers[1] # (i.e. conc1)
  id2 <- nested_identifiers[2] # (i.e. conc2)

  all_conc1 <- setdiff(unique(measured[, id]), 0)
  all_conc2 <- setdiff(unique(measured[, id2]), 0)

  # Single agent fits.
  single_agent1 <- fit_cotreatment_series(measured, series_id = id, cotrt_id = id2, cotrt_value = 0,
    normalization_type, all_conc = all_conc1)
  single_agent2 <- fit_cotreatment_series(measured, series_id = id2, cotrt_id = id, cotrt_value = 0,
    normalization_type, all_conc = all_conc2)

  # Fit cotreatments in matrix.
  cotrt_fittings <- vector("list", sum(length(all_conc1), length(all_conc2)))
  n <- 1
  for (conc1 in all_conc1) {
    cotrt_fittings[[n]] <- fit_cotreatment_series(measured, series_id = id2, cotrt_id = id,
      cotrt_value = conc1, all_conc = all_conc2, e_0 = single_agent2, GR_0 = single_agent2)
    n <- n + 1
  }
  for (conc2 in all_conc2) {
    cotrt_fittings[[n]] <- fit_cotreatment_series(measured, series_id = id, cotrt_id = id2,
      cotrt_value = conc2, all_conc = all_conc1, e_0 = single_agent1, GR_0 = single_agent1)
    n <- n + 1
  }

  cotrt_fittings
}


#' @keywords internal
fit_cotreatment_series <- function(flat_data, series_id, cotrt_id, cotrt_value, normalization_type, all_conc = NULL) {
  assert_numeric(cotrt_value)
  df_ <- flat_data[flat_data[[cotrt_id]] == cotrt_value & flat_data[[series_id]] > 0, , drop = FALSE]
  cotrt_fit <- gDRutils::fit_curves(
    df_ = df_[!is.na(df_[, normalization_type]), , drop = FALSE],
    series_identifiers = series_id,
    force_fit = TRUE,
    cap = 0.2,
    normalization_type = normalization_type
  )

  # Fill any missing concentration readouts.
  if (any(missing_concs <- setdiff(all_conc, unique(df_[[series_id]])))) {
    vals <- gDRutils::logistic_4parameters(
      missing_concs, 
      cotrt_fit$x_inf, 
      cotrt_fit$x_0, 
      cotrt_fit$ec50,
      cotrt_fit$h)

    df_inf <- data.frame(Concentration = missing_concs)
    df_inf[[norm_method]] <- vals
    missing_cols <- grepl("std_", colnames(df_))
    missing_cols <- c(missing_cols, colnames(df_) %in% c("CorrectedReadout", setdiff(c("RelativeViability", "GRvalue"), norm_method)))
    df_inf[missing_cols] <- 0 
    cotrt_fit <- rbind(cotrt_fit, df_inf)
  }

  # TODO: Include 2 Iranges columns of the concentration columns.
  cotrt_fit
}
