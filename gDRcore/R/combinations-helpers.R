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
  out <- gDRutils::predict_efficacy_from_conc(pred,
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
#' @param nested_identifiers character vector of identifiers in \code{measured} or \code{metric}
#' which define a unique data point.
#' @param metric_col string of the column in \code{metric} to use in the excess calculation.
#' @param measured_col string of the column in \code{measured} to use in the excess calculation.
#'
#' @return DataFrame of \code{measured}, now with an additional column named
#' \code{excess} (positive values for synergy/benefit).
#' @export
calculate_excess <- function(metric, measured, nested_identifiers, metric_col, measured_col) {
  # TODO: Ensure same dims metric, measured
  # TODO: Ensure there is a unique entry for series_id1, series_id2 and no repeats
  # TODO: Check that there are no NAs
  idx <- S4Vectors::match(DataFrame(measured[, nested_identifiers]), DataFrame(metric[, nested_identifiers]))
  
  out <- measured[, nested_identifiers]
  excess <- metric[idx, metric_col] - measured[, measured_col]
  excess[apply(as.matrix(out[, nested_identifiers]) == 0, 1, any)] <- NA
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


#' Standardize concentrations.
#'
#' Utilize a map to standardize concentrations.
#'
#' @param original_concs numeric vector of concentrations to replace using the \code{conc_map}.
#' @param conc_map data.frame of two columns named \code{original_conc_col} and \code{standardized_conc_col}.
#' @param original_conc_col string of the name of the column in \code{conc_map}
#' containing the original concentrations to replace.
#' @param standardized_conc_col string of the name of the column in \code{conc_map}
#' containing the standardized concentrations to use for replacement.
#'
#' @return numeric vector of standardized concentrations.
#'
#' @seealso map_conc_to_standardized_conc
#' @export
replace_conc_with_standardized_conc <- function(original_concs, conc_map, original_conc_col, standardized_conc_col) {
  out <- conc_map[match(original_concs, conc_map[[original_conc_col]]), standardized_conc_col]
  if (length(out) != length(original_concs)) {
    stop("standardized output is not the same length as the input")
  }
  out
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
#' @seealso replace_conc_w_standardized_conc
#' @details The concentrations are standardized in that they will contain regularly spaced dilutions
#' and close values will be rounded.
#' @export
map_conc_to_standardized_conc <- function(conc1, conc2) {
  # Remove single-agent.
  conc_1 <- setdiff(conc1, 0)
  conc_2 <- setdiff(conc2, 0)

  conc_1 <- sort(unique(conc_1))
  rconc1 <- .standardize_conc(conc_1)
  
  conc_2 <- sort(unique(conc_2))
  rconc2 <- .standardize_conc(conc_2)

  rconc <- c(0, unique(c(rconc1, rconc2)))
  .find_closest_match <- function(x) {
    rconc[abs(rconc - x) == min(abs(rconc - x))]
  }
  concs <- unique(c(conc1, conc2))
  mapped_rconcs <- vapply(concs, .find_closest_match, numeric(1))
  map <- unique(data.frame(concs = concs, rconcs = mapped_rconcs))

  tol <- 1
  
  # Check if standardized values are within 5% of the original values
  round_diff <- which(abs(map$concs - map$rconcs) / map$concs > 0.05)
  
  map$rconcs[round_diff] <- map$concs[round_diff]
  mismatched <- which(round_concentration(map$conc, tol) != round_concentration(map$rconc, tol))
  for (i in mismatched) {
    warning(sprintf("mapping original concentration '%s' to '%s'",
      map[i, "concs"], map[i, "rconcs"]))
  }
  map
}


#' Standardize concentration values.
#'
#' Standardize concentration values.
#'
#' @param conc numeric vector of the concentrations
#'
#' @return vector of standardized concentrations
#' @details If no \code{conc} are passed, \code{NULL} is returned.
#'
#' @export
.standardize_conc <- function(conc) {
  rconc <- if (S4Vectors::isEmpty(conc)) {
    NULL
  } else if (length(unique(round_concentration(conc, 3))) > 4) {
    # 4 is determined by the fewest number of concentrations required to be considered a "matrix".
    log10_step <- .calculate_dilution_ratio(conc)
    num_steps <- round((log10(max(conc)) - log10(min(conc)) / log10_step), 0)
    seqs <- log10(max(conc)) - (log10_step * c(0:num_steps))
    sort(round_concentration(10 ^ seqs, 3))
  } else {
    # Few enough concentrations, don't need to calculate a series.
    round_concentration(conc, 3)
  }
  rconc
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
  rounded_dil_ratios <- round_concentration(dil_ratios, 2)

  # Get most frequent dilution ratio.
  highest_freq_ratio <- names(sort(table(rounded_dil_ratios), decreasing = TRUE))[[1]]
  as.numeric(highest_freq_ratio)
}
