#' Get predicted values for a given fit and input.
#'
#' Map fittings to identifiers and compute the predicted values for corresponding fits.
#'
#' @param pred numeric vector for which you want predictions.
#' @param match_col vector to match on \code{fittings} to get the correct fit.
#' @param fittings data.frame of fit metrics.
#' @param fitting_id_col string of the column name in \code{fittings} that should be
#' used to match with \code{match_col} .
#'
#' @return numeric vector of predicted values given \code{pred} inputs and \code{fittings} values.
#'
#' @export
map_ids_to_fits <- function(pred, match_col, fittings, fitting_id_col) {
  ridx <- S4Vectors::match(round(log10(match_col), 2), round(log10(fittings[[fitting_id_col]]), 2))
  metrics <- fittings[ridx, c(fitting_id_col, "x_inf", "x_0", "ec50", "h")]
  # Extrapolate fitted values.
  out <- gDRutils::logistic_4parameters(pred,
    metrics$x_inf,
    metrics$x_0,
    metrics$ec50,
    metrics$h)
  out
}


#' Calculate the difference between values in two data.frames
#'
#' Calculate the difference between values, likely representing the same metric, from two data.frames.
#'
#' @param metric data.frame often representing readouts derived by calculating some metric.
#' Examples of this could include hsa or bliss calculations from single-agent data. 
#' @param measured data.frame often representing measured data from an experiment.
#' @param series_identifiers character vector of identifiers in \code{measured} or \code{metric}
#' which define a unique data point.
#' @param metric_col string of the column in \code{metric} to use in the excess calculation.
#' @param measured_col string of the column in \code{measured} to use in the excess calculation.
#'
#' @return DataFrame of \code{measured}, now with an additional column named
#' \code{excess} (positive values for synergy/benefit).
#' @export
calculate_excess <- function(metric, measured, series_identifiers, metric_col, measured_col) {
  # TODO: Ensure same dims metric, measured
  # TODO: Ensure there is a unique entry for series_id1, series_id2 and no repeats
  # TODO: Check that there are no NAs
  idx <- S4Vectors::match(DataFrame(measured[, series_identifiers]), DataFrame(metric[, series_identifiers]))
  
  out <- measured[, series_identifiers]
  excess <- metric[idx, metric_col] - measured[, measured_col]
  excess[apply(as.matrix(out[, series_identifiers]) == 0, 1, any)] <- 0 
  out$excess <- excess
  as.data.frame(out)
}


#' @keywords internal
#' @noRd
convertDFtoBumpyMatrixUsingIds <- function(df, row_id = "row_id", col_id = "col_id") {
  BumpyMatrix::splitAsBumpyMatrix(
    df[!colnames(df) %in% c(row_id, col_id)],
    row = df[[row_id]], col = df[[col_id]])
}


#' Create a mapping of concentrations to standardized concentrations.
#'
#' Create a mapping of concentrations to standardized concentrations.
#'
#' @param conc1 numeric vector of the concentrations for drug 1.
#' @param conc2 numeric vector of the concentrations for drug 2.
#'
#' @return data.frame of 2 columns named \code{"concs"} and \code{"rconcs"}
#' containing the original concentrations and their closest matched standardized concentrations
#' respectively.
#' and their new standardized concentrations.
#'
#' @details The concentrations are standardized in that they will contain regularly spaced dilutions
#' and close values will be rounded.
#' @export
map_conc_to_standardized_conc <- function(conc1, conc2) {
  # Remove any single-agent data.
  conc1 <- conc1[conc1 > 0]
  conc2 <- conc2[conc2 > 0]
  concs <- c(conc1, conc2)

  conc_1 <- sort(unique(conc1))
  log10_step1 <- .calculate_dilution_ratio(conc_1)
  rconc_1 <- sort(round_concentration(10 ^ seq(log10(max(conc_1)),
                                               log10(min(conc_1)),
                                               -log10_step1), 3))

  conc_2 <- sort(unique(conc2[conc1 > 0]))
  log10_step2 <- .calculate_dilution_ratio(conc_2)
  rconc_2 <- sort(round_concentration(10 ^ seq(log10(max(conc_2)),
                                               log10(min(conc_2)),
                                               -log10_step2), 3))
  rconc <- unique(round_concentration(c(rconc_1, rconc_2), 3))
  .find_closest_match <- function(x) {
    rconc[abs(rconc - x) == min(abs(rconc - x))]
  }
  mapped_rconcs <- vapply(concs, .find_closest_match, numeric(1))

  map <- data.frame(concs = concs, rconcs = mapped_rconcs)

  mismatched <- which(round_concentration(map$conc, 2) != round_concentration(map$rconc, 2))
  for (i in mismatched) {
    warning(sprintf("mapped original concentration '%s' to '%s'",
      map[i, "concs"], map[i, "rconcs"]))
  }

  map
}



#' Calculate a dilution ratio.
#'
#' Calculation a dilution ratio from a vector of concentrations.
#'
#' @param concs numeric vector of concentrations.
#'
#' @return numeric value of the dilution ratio for a given set of concentrations.
#' @keywords internal
#' @noRd
.calculate_dilution_ratio <- function(concs) {
  first_removed <- concs[-1]
  first_two_removed <- first_removed[-1]
  last_removed <- concs[-length(concs)]
  last_two_removed <- last_removed[-length(last_removed)]

  dil_ratios <- c(log10(first_removed / last_removed), log10(first_two_removed / last_two_removed))
  rounded_dil_ratios <- round_concentration(dil_ratios, 3)

  # Get most frequent dilution ratio.
  highest_freq_ratio <- names(sort(table(rounded_dil_ratios), decreasing = TRUE))[[1]]
  as.numeric(highest_freq_ratio)
}
