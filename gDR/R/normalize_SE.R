#' normalize_SE
#'
#' Normalize drug response data from treated and untreated pairings.
#'
#' @param se \code{BumpyMatrix} object with assays \code{"RawTreated"} and \code{"Controls"}.
#' @param nDigits_rounding integer specifying number of digits of rounding during calculations.
#' Defaults to \code{4}.
#'
#' @return \code{BumpyMatrix} object with a new assay named \code{"Normalized"} containing \code{DataFrame}s 
#' holding \code{RelativeViability}, \code{GRvalue}, \code{RefRelativeViability}, code{RefGRvalue}, and \code{DivisionTime} values.
#'
#' @export
#'
normalize_SE <- function(se, nDigits_rounding = 4) {
  # Assertions
  checkmate::assert_number(nDigits_rounding)

  duration_col <- gDRutils::get_identifier("duration") # TODO: See normalization call below.

  refs <- SummarizedExperiment::assays(se)[["Controls"]]
  trt <- SummarizedExperiment::assays(se)[["RawTreated"]]

  # TODO: Remove looping and just use the BumpyMatrix to do all of the arithmetic.  
  # This is simple for the relative viability, but we will need to think harder about the GR value calculations.
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

      controls <- normalized <- S4Vectors::DataFrame()

      # Normalized treated.
      normalized$RelativeViability <- round(trt_df$CorrectedReadout/ref_df$UntrtReadout, ndigits_rounding)
      normalized$GRvalue <- calculate_GR_value(trt_df, ndigits_rounding)

      # Normalized references.
      normalized$RefRelativeViability <- round(ref_df$RefReadout/normalized$UntrtReadout, nDigits_rounding)
      normalized$RefGRvalue <- calculate_GR_value(ref_df, ndigits_rounding)
      normalized$DivisionTime <- round(rdata[i, duration_col] / log2(ref_df$UntrtReadout/ref_df$Day0Readout), nDigits_rounding)

      # TODO: Put the normalized data.frame into the bumpy matrix. 
      norm[i, j] <- normalized
    }
  }

  # TODO: Put the bumpy matrix assay back into the SE as a new "Normalized" assay.
  assays(se)[["Normalized"]] <- norm
  return(se)
}
