#' @keywords internal
fit_codilutions <- function(measured, nested_identifiers, normalization_type) {
  id <- nested_identifiers[1]
  id2 <- nested_identifiers[2]

  # Filter out all single-agents.
  single_agents <- measured[[id]] == 0 | measured[[id2]] == 0
  measured <- measured[!single_agents, , drop = FALSE]

  # Filter only to what diagonals are valid for a fit (>4 points).
  measured$ratios <- measured[[id2]]/measured[[id]]
  ratios <- S4Vectors::split(measured, measured$ratios)
  keep <- unlist(lapply(ratios, function(x) {nrow(x) > 4}))
  valid <- ratios[keep]

  fits <- vector("list", length(valid))
  for (i in seq_along(fits)) {
    fits[[i]] <- fit_codilution_series(valid[[i]], id, id2, e_0 = 1, GR_0 = 1, normalization_type)
  }

  out <- do.call("rbind", fits)
  out
}


#' @keywords internal
fit_codilution_series <- function(measured, series_1, series_2, e_0, GR_0, normalization_type, summed_col = "summed_codil_conc") {
  if (length(unique(measured[[series_2]] / measured[[series_1]])) != 1L) {
    stop("more than one ratio between 'series_2' and 'series_1' detected")
  }

  measured[summed_col] <- measured[[series_1]] + measured[[series_2]]
  keep <- setdiff(colnames(measured), c(series_1, series_2))
  codilution_fit <- gDRutils::fit_curves(
    df_ = measured[, keep, drop = FALSE],
    series_identifiers = summed_col,
    e_0 = e_0,
    GR_0 = GR_0,
    force_fit = TRUE,
    cap = 0.2,
    normalization_type = normalization_type
  )

  out <- DataFrame(IRanges::DataFrameList(measured[, c(series_1, series_2)]), codilution_fit)
  out
}
