#' @keywords internal
fit_combo_codilutions <- function(measured,
                                  series_identifiers,
                                  normalization_type,
                                  e_0 = 1,
                                  GR_0 = 1) {
  id <- series_identifiers[1]
  id2 <- series_identifiers[2]

  # Filter out all single-agents.
  single_agents <- measured[[id]] == 0 | measured[[id2]] == 0
  measured <- measured[!single_agents, , drop = FALSE]
  if (nrow(measured) < 4) {
    return(NULL)
  }

  # Filter only to what diagonals are valid for a fit (>4 points).
  measured$ratios <- round_concentration(measured[[id2]] / measured[[id]], ndigit = 1)
  ratios <- S4Vectors::split(measured, measured$ratios)
  keep <- unlist(lapply(ratios, function(x) {
    nrow(x) > 4
  }))
  valid <- ratios[keep]

  fits <- vector("list", length(valid))
  for (i in seq_along(fits)) {
    fits[[i]] <- fit_codilution_series(valid[[i]], id, id2, e_0 = e_0, GR_0 = GR_0, normalization_type)
  }

  out <- do.call("rbind", fits)
  out
}


#' @keywords internal
fit_codilution_series <- function(measured, series_1, series_2, e_0, GR_0, normalization_type) {
  ratio <- unique(round_concentration(measured[[series_2]] / measured[[series_1]], ndigit = 1))
  if (length(ratio) != 1L) {
    stop("more than one ratio between 'series_2' and 'series_1' detected")
  }

  measured$summed_conc <- measured[[series_1]] + measured[[series_2]]
  keep <- setdiff(colnames(measured), c(series_1, series_2))
  codilution_fit <- gDRutils::fit_curves(
    df_ = measured[, keep, drop = FALSE],
    series_identifiers = "summed_conc",
    e_0 = e_0,
    GR_0 = GR_0,
    force_fit = TRUE,
    cap = 0.2,
    normalization_type = normalization_type
  )

  codilution_fit$ratio <- ratio 
  codilution_fit
}
