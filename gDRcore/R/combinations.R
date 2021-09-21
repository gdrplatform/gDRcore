#' @keywords internal
.create_combo_control <- function(nested_identifiers) {
  out <- data.frame()
  for (i in nested_identifiers) {
    out[[i]] <- 0
  }
  out$GRvalue <- 1
  out$std_GRvalue <- 0
  out$RelativeViability <- 1
  out$std_RelativeViability <- 0
  out
}

map_ids_to_fits <- function(ids, fittings, fitting_id_col) {
  ridx <- S4Vectors::match(round(log10(ids), 2), round(log10(fittings[[fitting_id_col]]), 2))
  metrics <- fittings[ridx, c(fitting_id_col, "x_inf", "x_0", "ec50", "h")]
  # Extrapolate fitted values.
  out <- gDRutils::logistic_4parameters(metrics[[fitting_id_col]],
    metrics$x_inf,
    metrics$x_0,
    metrics$ec50,
    metrics$h)
  out
}


#' @details HSA takes the minimum of the two single agents readouts.
calculate_HSA <- function(sa1, series_id1, sa2, series_id2, metric) {
  .calculate_matrix_metric(sa1, series_id1, sa2, series_id2, metric, FXN = pmin)
}


#' @details Bliss is the mulitplication of the single agent readouts.
calculate_Bliss <- function(sa1, series_id1, sa2, series_id2, metric) {
  .calculate_matrix_metric(sa1, series_id1, sa2, series_id2, metric, FXN = "*")
}


#' @keywords internal
.calculate_matrix_metric <- function(sa1, series_id1, sa2, series_id2, metric, FXN) {
  colnames(sa1)[colnames(sa1) == metric] <- "metric1"
  colnames(sa2)[colnames(sa2) == metric] <- "metric2"

  # TODO: ensure they're unique?
  u <- expand.grid(sa1[[series_id1]], sa2[[series_id2]])
  colnames(u) <- c(series_id1, series_id2)

  idx <- match(u[[series_id1]], sa1[[series_id1]])
  u <- base::merge(u, sa1[, c(series_id1, "metric1")], by = series_id1)
  u <- base::merge(u, sa2[, c(series_id2, "metric2")], by = series_id2)

  metric <- do.call(FXN, list(u$metric1, u$metric2))
  cbind(u, metric)
}


#' Calculate the difference between values in two data.frames
#'
#' Calculate the difference between values, likely representing the same metric, from two data.frames.
#'
#' @param metric data.frame often representing readouts derived by calculating some metric.
#' Examples of this could include hsa or bliss calculations from single-agent data. 
#' @param measured data.frame often representing measured data from an experiment.
#' @param series_identifiers character vector of identifiers in \code{measured} or \code{metric} which define a unique data point.
#' @param metric_col string of the column in \code{metric} to use in the excess calculation.
#' @param measured_col string of the column in \code{measured} to use in the excess calculation.
#'
#' @return DataFrame of \code{measured}, now with an additional column named \code{excess}.
#' @export
calculate_excess <- function(metric, measured, series_identifiers, metric_col, measured_col) {
  # TODO: Ensure same dims metric, measured
  # TODO: Ensure there is a unique entry for series_id1, series_id2 and no repeats
  # TODO: Check that there are no NAs
  idx <- S4Vectors::match(DataFrame(metric[, series_identifiers]), DataFrame(measured[, series_identifiers]))
  
  out <- measured[, series_identifiers]
  excess <- metric[idx, metric_col] - measured[, measured_col]
  out$excess <- excess
  out
}


calculate_Loewe <- function(mean_matrix, row_fittings, col_fittings, dilution_fittings) {
  if (min(mean_matrix, na.rm = TRUE) > 0.7) {
    iso_cutoff <- NULL
  } else {
    if (norm_method == "GRvalue") {
      iso_cutoff <- seq(max(-0.25, ceiling(20 * min(mx_mean + 0.08, na.rm = TRUE)) / 20), 0.8,  0.05)
    } else {
      iso_cutoff <- seq(max(0.2, ceiling(20 * min(mx_mean + 0.08, na.rm = TRUE)) / 20), 0.8,  0.05)
    }
    names(iso_cutoff) <- as.character(iso_cutoff)
  }

  # create the variable for the different isobolograms
  all_iso <- vector("list", length(iso_cutoff))
  # TODO: Should these xc50s become ec50s?
  # TOOD: by_col and by_row should become encoded by the IRanges object length
  ref_x50 <- c(conc_1 = min(all_fits[["by_col"]][1, "xc50"], max(drug1_axis$conc_1) * conc_margin),
              conc_2 = min(all_fits[["by_row"]][1, "xc50"], max(drug2_axis$conc_2) * conc_margin))
  names(all_iso) <- iso_cutoff
  for (isobol_value in iso_cutoff) { # run through the different isobolograms
    # cutoff point by row
    row_fittings <- row_fittings[nrow(row_fittings):1, ]
    df_iso <- cbind(row_fittings[, "conc_1", drop = FALSE], data.frame(conc_2 =
      ifelse(row_fittings$x_0 < isobol_value, 0, ifelse(row_fittings$x_inf > isobol_value,
        Inf,
        row_fittings$ec50 * ((((row_fittings$x_0 - row_fittings$x_inf) / (isobol_value - row_fittings$x_inf)) - 1) ^
            (1 / pmax(row_fittings$h, 0.01))))
          ),
          fit_type = "by_row"))

    # cutoff point by columns
    col_fittings <- all_fits[["by_col"]]
    df_iso <- rbind(df_iso, cbind(data.frame(conc_1 =
      ifelse(col_fittings$x_0 < isobol_value, 0, ifelse(col_fittings$x_inf > isobol_value,
        Inf,
        col_fittings$ec50 * ((((col_fittings$x_0 - col_fittings$x_inf) / (isobol_value - col_fittings$x_inf)) - 1) ^
            (1 / pmax(col_fittings$h, 0.01))))
          ),
      col_fittings[, "conc_2", drop = FALSE],
        fit_type = "by_col")))


    df_iso$conc_1 <- pmin(df_iso$conc_1, max(drug1_axis$conc_1) * conc_margin)
    df_iso$conc_2 <- pmin(df_iso$conc_2, max(drug2_axis$conc_2) * conc_margin)

    ref_conc_1 <- pmin(df_iso$conc_1[df_iso$conc_2 == 0 & df_iso$fit_type == "by_col"],
                    max(drug1_axis$conc_1) * conc_margin)
    ref_conc_2 <- pmin(df_iso$conc_2[df_iso$conc_1 == 0 & df_iso$fit_type == "by_row"],
                    max(drug2_axis$conc_2) * conc_margin)

    # cutoff point by diagonal (co-dilution)
    # co-dil is given as concentration of drug 1
    if (nrow(codilution_fittings) > 1) {
      codilution_fittings <- codilution_fittings[nrow(codilution_fittings):1, ]
      codilution_fittings <- codilution_fittings[codilution_fittings$fit_type %in% "DRC3pHillFitModelFixS0", ]
      conc_mix <- if (codilution_fittings$x_0 < isobol_value) {
        0
        } else {
          ifelse(codilution_fittings$x_inf > isobol_value,
        Inf,
        codilution_fittings$ec50 * ((((codilution_fittings$x_0 - codilution_fittings$x_inf) /
                                        (isobol_value - codilution_fittings$x_inf)) - 1) ^
                                      (1 / codilution_fittings$h)))
        }  
      df_iso_codil <- data.frame(conc_1 = conc_mix / (1 + 1 / codilution_fittings$conc_ratio),
                        conc_2 = conc_mix / (1 + codilution_fittings$conc_ratio), fit_type = "by_codil")
      # avoid extrapolation
      capped_idx <- df_iso_codil$conc_1 > (max(drug1_axis$conc_1) * conc_margin) |
            df_iso_codil$conc_2 > (max(drug2_axis$conc_2) * conc_margin)
      df_iso_codil$conc_1[capped_idx] <- pmin(max(drug1_axis$conc_1) * conc_margin,
                  codilution_fittings$conc_ratio * (max(drug2_axis$conc_2) * conc_margin))[capped_idx]
      df_iso_codil$conc_2[capped_idx] <- pmin(max(drug2_axis$conc_2) * conc_margin,
                  (max(drug1_axis$conc_1) * conc_margin) / codilution_fittings$conc_ratio)[capped_idx]

      df_iso <- rbind(df_iso, df_iso_codil)
    }
    # remove low concentration values
    df_iso <- df_iso[!is.na(df_iso$conc_1) & !is.na(df_iso$conc_2), ]
    df_iso <- df_iso[(df_iso$conc_1 > drug1_axis$conc_1[2] / 2 | df_iso$fit_type == "by_row") &
                       (df_iso$conc_2 > drug2_axis$conc_2[2] / 2 | df_iso$fit_type == "by_col"), ]

    if (nrow(df_iso) < 5) {
      next
    }

    # perform the smoothing
    df_iso_curve <- data.frame(x1 = isobol_x1, x2_off = zoo::rollmean(
      rowMeans(do.call(cbind, lapply(names(which(table(df_iso$fit_type) > 1)), function(x)
        stats::approx(x = df_iso$x1[df_iso$fit_type == x], 
               y = df_iso$x2_off[df_iso$fit_type == x], xout = isobol_x1)$y)),
        na.rm = TRUE), 5, fill = NA))
    df_iso_curve <- df_iso_curve[!is.na(df_iso_curve$x2_off), ]
    df_iso_curve$x2 <- df_iso_curve$x2_off - abs(df_iso_curve$x1) * x2_extra_offset

    # rotate back the position
    df_iso_curve$pos_x <- (df_iso_curve$x1 + df_iso_curve$x2) / sqrt(2) + min(drug2_axis$pos_x)
    df_iso_curve$pos_y <- (-df_iso_curve$x1 + df_iso_curve$x2) / sqrt(2) + min(drug1_axis$pos_y)

    # calculate the reference (additive model in the rotated space)
    c2 <- ref_conc_2 / (1 + (ref_conc_2 / ref_conc_1) * (10 ^ (
                -(sqrt(2) * df_iso_curve$x1 + min(drug2_axis$pos_x) - min(drug1_axis$pos_y)))))
    df_iso_curve$x2_ref <- (log10(c2) - min(drug2_axis$pos_x) +
              (log10(ref_conc_1 * (1 - c2 / ref_conc_2)) - min(drug1_axis$pos_y))) / sqrt(2)

    # cap the concentrations for the reference
    over_edge <- pmax(0, (-df_iso_curve$x1 + df_iso_curve$x2_ref) / sqrt(2) +
                  min(drug1_axis$pos_y) -  (max(drug1_axis$pos_y) + conc_margin)) +
                pmax(0, (df_iso_curve$x1 + df_iso_curve$x2_ref) / sqrt(2) +
                  min(drug2_axis$pos_x) -  (max(drug2_axis$pos_x) + conc_margin))
    df_iso_curve$x2_ref_cap <- df_iso_curve$x2_ref - over_edge * sqrt(2)

    # rotate back the reference
    df_iso_curve$pos_x_ref <- (df_iso_curve$x1 + df_iso_curve$x2_ref) / sqrt(2) +
      min(drug2_axis$pos_x)
    df_iso_curve$pos_y_ref <- (-df_iso_curve$x1 + df_iso_curve$x2_ref) / sqrt(2) +
      min(drug1_axis$pos_y)
    df_iso_curve <- rbind(NA, df_iso_curve, NA)
    df_iso_curve[1, c("pos_x", "pos_x_ref")] <- min(drug2_axis$pos_x)
    df_iso_curve[1, c("pos_y", "pos_y_ref")] <- log10(ref_conc_1)
    df_iso_curve[1, "x1"] <-
                (df_iso_curve$pos_x[1] - min(drug2_axis$pos_x) -
                    (df_iso_curve$pos_y[1] - min(drug1_axis$pos_y))) / sqrt(2)
    df_iso_curve[1, c("x2", "x2_ref", "x2_ref_cap")] <-
                (df_iso_curve$pos_x[1] - min(drug2_axis$pos_x) +
                    (df_iso_curve$pos_y[1] - min(drug1_axis$pos_y))) / sqrt(2)

    df_iso_curve[nrow(df_iso_curve), c("pos_x", "pos_x_ref")] <- log10(ref_conc_2)
    df_iso_curve[nrow(df_iso_curve), c("pos_y", "pos_y_ref")] <- min(drug1_axis$pos_y)
    df_iso_curve[nrow(df_iso_curve), c("x2", "x2_ref", "x2_ref_cap")] <-
                (df_iso_curve$pos_x[nrow(df_iso_curve)] - min(drug2_axis$pos_x) +
                    (df_iso_curve$pos_y[nrow(df_iso_curve)] -
                    min(drug1_axis$pos_y))) / sqrt(2)
    df_iso_curve[nrow(df_iso_curve), "x1"] <-
                (df_iso_curve$pos_x[nrow(df_iso_curve)] - min(drug2_axis$pos_x) -
                    (df_iso_curve$pos_y[nrow(df_iso_curve)] -
                    min(drug1_axis$pos_y))) / sqrt(2)

    # calculate CI across range to concentration ratios (in the rotated space)
    df_iso_curve$log10_ratio_conc <- df_iso_curve$x1
    df_iso_curve$log2_CI <- zoo::rollmean(
        log2(10) * (df_iso_curve$x2 - df_iso_curve$x2_ref) / sqrt(2), 4,
        fill = c(0, 0, 0))

    # cap position for plotting the isobolograms
    df_iso$pos_y <- pmin(df_iso$pos_y, max(drug1_axis$pos_y) + log2_pos_offset)
    df_iso$pos_x <- pmin(df_iso$pos_x, max(drug2_axis$pos_x) + log2_pos_offset)
    df_iso_curve$pos_y <- pmin(df_iso_curve$pos_y, max(drug1_axis$pos_y) + log2_pos_offset)
    df_iso_curve$pos_x <- pmin(df_iso_curve$pos_x, max(drug2_axis$pos_x) + log2_pos_offset)
    df_iso_curve$pos_y_ref <- pmin(df_iso_curve$pos_y_ref,
        max(drug1_axis$pos_y) + log2_pos_offset)
    df_iso_curve$pos_x_ref <- pmin(df_iso_curve$pos_x_ref,
        max(drug2_axis$pos_x) + log2_pos_offset)

    range <- 2 # in log10 domain --> 100-fold range for calculating averaged CI 
    ratio_idx <- which(
        (df_iso_curve$log10_ratio_conc > (min(df_iso_curve$log10_ratio_conc) + range / 2)) &
        (df_iso_curve$log10_ratio_conc < (max(df_iso_curve$log10_ratio_conc) - range / 2)))

    if (length(ratio_idx) == 0) {
      ratio_idx <- round(nrow(df_iso_curve) / 2)
    }

    df_100x_AUC <- data.frame(log10_ratio_conc = df_iso_curve$log10_ratio_conc[ratio_idx],
      AUC_CI = vapply(ratio_idx, function(x) mean(df_iso_curve$log2_CI[
        (df_iso_curve$log10_ratio_conc > (df_iso_curve$log10_ratio_conc[x] - range / 2)) &
        (df_iso_curve$log10_ratio_conc <= (df_iso_curve$log10_ratio_conc[x] + range / 2))]),
        FUN.VALUE = double(1)
      )
    )

    all_iso[[as.character(isobol_value)]] <- list(
      df_iso = df_iso,
      df_iso_curve = df_iso_curve,
      ref_conc_1 = ref_conc_1,
      ref_conc_2 = ref_conc_2,
      max_log2CI = max(df_iso_curve$log2_CI),
      min_log2CI = min(df_iso_curve$log2_CI),
      CI_100x = min(df_100x_AUC$AUC_CI)
      )
    }
    all_iso <- all_iso[!vapply(all_iso, is.null, FUN.VALUE = logical(1))]

    df_CI_100x <- unlist(lapply(all_iso, "[[", "CI_100x"))
    df_CI_100x <- data.frame(level = as.numeric(names(df_CI_100x)), log2_CI = df_CI_100x)
}
