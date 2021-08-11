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


calculate_combo_mean <- function(row_fittings, col_fittings, codilution_fittings, nested_identifiers) {
  if (dim(row_fittings) != dim(col_fittings) || dim(col_fittings) != dim(codilution_fittings)) {
    stop("all input fittings must have the same dimensions")
  }

  merged <- merge(row_fittings, col_fittings, by = nested_identifiers)
  merged <- merge(merged, codilution_fittings, by = nested_identifiers)

  merged$average <- rowMeans(as.matrix(merged[, !colnames(merged) %in% nested_identifiers]))
  merged
}


#' @details HSA takes the minimum of the two single agents readouts.
calculate_HSA <- function(mean_matrix) {
  .calculate_matrix_metric(mean_matrix, pmin)
}


#' @details Bliss is the mulitplication of the single agent readouts.
calculate_Bliss <- function(mean_matrix) {
  .calculate_matrix_metric(mean_matrix, "*")
}


#' @keywords internal
# TODO: Fix me so I work off a DataFrame instead of a matrix.
.calculate_matrix_metric <- function(mean_matrix, FXN) {
  single_agent1 <- mean_matrix[-1, 1]
  single_agent2 <- mean_matrix[1, -1]

  m <- nrow(mean_matrix)
  n <- ncol(mean_matrix)

  mat <- matrix(NA, m - 1, n - 1)
  sa_mat1 <- t(matrix(t(single_agent2), n - 1, m - 1))
  sa_mat2 <- matrix(single_agent1, m - 1, n - 1)
  mat <- FXN(sa_mat1, sa_mat2)
  mat
}


calculate_Loewe <- function(mean_matrix) {
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
    df_fit <- all_fits[["by_row"]]
    df_fit <- df_fit[nrow(df_fit):1, ]
    df_iso <- cbind(df_fit[, "conc_1", drop = FALSE], data.frame(conc_2 =
      ifelse(df_fit$x_0 < isobol_value, 0, ifelse(df_fit$x_inf > isobol_value,
        Inf,
        df_fit$ec50 * ((((df_fit$x_0 - df_fit$x_inf) / (isobol_value - df_fit$x_inf)) - 1) ^
            (1 / pmax(df_fit$h, 0.01))))
          ),
          fit_type = "by_row"))

    # cutoff point by columns
    df_fit <- all_fits[["by_col"]]
    df_iso <- rbind(df_iso, cbind(data.frame(conc_1 =
      ifelse(df_fit$x_0 < isobol_value, 0, ifelse(df_fit$x_inf > isobol_value,
        Inf,
        df_fit$ec50 * ((((df_fit$x_0 - df_fit$x_inf) / (isobol_value - df_fit$x_inf)) - 1) ^
            (1 / pmax(df_fit$h, 0.01))))
          ),
      df_fit[, "conc_2", drop = FALSE],
        fit_type = "by_col")))


    df_iso$conc_1 <- pmin(df_iso$conc_1, max(drug1_axis$conc_1) * conc_margin)
    df_iso$conc_2 <- pmin(df_iso$conc_2, max(drug2_axis$conc_2) * conc_margin)

    ref_conc_1 <- pmin(df_iso$conc_1[df_iso$conc_2 == 0 & df_iso$fit_type == "by_col"],
                    max(drug1_axis$conc_1) * conc_margin)
    ref_conc_2 <- pmin(df_iso$conc_2[df_iso$conc_1 == 0 & df_iso$fit_type == "by_row"],
                    max(drug2_axis$conc_2) * conc_margin)

    # cutoff point by diagonal (co-dilution)
    # co-dil is given as concentration of drug 1
    df_fit <- all_fits[["by_codil"]]
    if (nrow(df_fit) > 1) {
      df_fit <- df_fit[nrow(df_fit):1, ]
      df_fit <- df_fit[df_fit$fit_type %in% "DRC3pHillFitModelFixS0", ]
      conc_mix <- ifelse(df_fit$x_0 < isobol_value, 0, ifelse(df_fit$x_inf > isobol_value,
        Inf,
        df_fit$ec50 * ((((df_fit$x_0 - df_fit$x_inf) / (isobol_value - df_fit$x_inf)) - 1) ^
            (1 / df_fit$h)))
          )
      df_iso_codil <- data.frame(conc_1 = conc_mix / (1 + 1 / df_fit$conc_ratio),
                        conc_2 = conc_mix / (1 + df_fit$conc_ratio), fit_type = "by_codil")
      # avoid extrapolation
      capped_idx <- df_iso_codil$conc_1 > (max(drug1_axis$conc_1) * conc_margin) |
            df_iso_codil$conc_2 > (max(drug2_axis$conc_2) * conc_margin)
      df_iso_codil$conc_1[capped_idx] <- pmin(max(drug1_axis$conc_1) * conc_margin,
                  df_fit$conc_ratio * (max(drug2_axis$conc_2) * conc_margin))[capped_idx]
      df_iso_codil$conc_2[capped_idx] <- pmin(max(drug2_axis$conc_2) * conc_margin,
                  (max(drug1_axis$conc_1) * conc_margin) / df_fit$conc_ratio)[capped_idx]

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
