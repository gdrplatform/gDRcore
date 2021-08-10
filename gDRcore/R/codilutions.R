#' @keywords internal
fit_codilutions <- function(measured, nested_identifiers, normalization_type) {
  id <- nested_identifiers[1]
  id2 <- nested_identifiers[2]

  # Filter out all single-agents.
  single_agents <- measured[[id]] == 0 | measured[[id2]] == 0
  measured <- measured[!single_agents, , drop = FALSE]

  # Calculate concentration ratios.
  measured$ratios <- measured[[id2]]/measured[[id]] # TODO: Fix to cacluation the same ratios as below?
  ##n_conc_codil <- min(ncol(mx_response), nrow(mx_response)) - 1
  ##conc_ratio <- 2 ^ unique(round(log2(
  ##  as.numeric(colnames(mx_response))[ncol(mx_response):(ncol(mx_response) - n_conc_codil + 1)] / 
  ##    as.numeric(rownames(mx_response))[nrow(mx_response):(nrow(mx_response) - n_conc_codil + 1)]), 1))

  # Split by ratios.


  # Filter only to what diagonals are valid for a fit. (>4 points)
  valid <- lapply(, function(x) {length(x) > 4})

  fits <- vector("list", length(valid))
  for (i in seq_along(fits)) {
    # Fit the diagonal.
    # TODO: Call the wrapper function to fit_curves which also binds the IRanges object.
    fit_res <- gDRutils::fit_curves(
      df_,
      series_identifiers = series_identifiers,
      e_0 = 1,
      GR_0 = 1,
      normalization_type = normalization_type,
      force_fit = TRUE
    )
    fits[[i]] <- fit_res
  }

  out <- do.call("rbind", fits)
  out
}
