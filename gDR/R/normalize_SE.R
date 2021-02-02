#' normalize_SE
#'
#' Normalize drug response data from treated and untreated pairings.
#'
#' @param se \code{BumpyMatrix} object with assays \code{"treated"} and \code{Controls}.
#' @param nDigits_rounding integer specifying number of digits of rounding during calculations.
#' Defaults to \code{4}.
#'
#' @return \code{BumpyMatrix} object with a new assay named \code{"Normalized"} containing \code{DataFrame}s 
#' holding \code{RelativeViability} and \code{GRvalue} values.
#'
#' @export
#'
normalize_SE <- function(se, nDigits_rounding = 4) {
  # Assertions
  checkmate::assert_number(nDigits_rounding)

  duration_col <- gDRutils::get_identifier("duration") # TODO: See normalization call below.

  refs <- SummarizedExperiment::assays(se)[["Controls"]]
  trt <- SummarizedExperiment::assays(se)[["treated"]]

  # TODO: Create empty BM? 
  # bm <- create_empty_bm()
  for (i in rownames(se)) {
    for (j in colnames(se)) {
      ref_df <- refs[i, j][[1]]
      trt_df <- trt[i, j][[1]]

      if (nrow(trt_df) == 0L) {
	next # skip if no data
      }

      if (length(ref_df) == 0L) {
	futile.logger::flog.warn("Missing control data. Treatment Id: '%s' Cell_line Id: '%s'", i, j)
	next
      }

      # BumpyMatrix object has unmutable set of columns in DataFrame for each [i,j]
      # thus se assay was initialized with NAs for "GRvalue" and "RelativeViability"
      normalized <- normalize_trt_to_ref(trt_df, ref_df, ndigits_rounding = nDigits_rounding) # TODO: Potentially add the duration_col as an arg.
      # TODO: Put the normalized data.frame into the bumpy matrix. 
      #bm[i, j] <- normalized
    }
  }

  # TODO: Put the bumpy matrix assay back into the SE as a new "Normalized" assay.
  return(se)
}


#' Normalize a set of treated drug response values to its corresponding reference drug response values. 
#'
#' Normalize treated and reference pairings to obtain key metrics on a per-concentration basis.
#'
#' @param trt_df data.frame of the treated elements.
#' @param ref_df data.frame of the reference elements to normalize against.
#' @param ndigits_rounding integer of the number of digits to round to during normalization calculations.
#'
#' @return DataFrame of containing the \code{RelativeViability}, \code{GRvalue}, adn \code{CorrectedReadout}.
#'
#' @export
#'
normalize_trt_to_ref <- function(trt_df, ref_df, ndigits_rounding) {
  ## Normalize to the untreated readout, which should have the highest readout,
  ## to get relative viability values between [0, 1]. 
  normalized <- S4Vectors::DataFrame()
  normalized$RelativeViability <- round(trt_df$CorrectedReadout/ref_df$UntrtReadout, ndigits_rounding)
  normalized$GRvalue <- calculate_GR_value(trt_df, ndigits_rounding)
  normalized$CorrectedReadout <- trt_df$CorrectedReadout # TODO: Within average_SE(), 
							 # it looks like we need: c("GRvalue", "RelativeViability","CorrectedReadout"). 
							 # I don't actually see the "CorrectedReadout" being used though.
  normalized
}
