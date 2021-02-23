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
  trt <- SummarizedExperiment::assays(se)[["RawTreated"]]

  # TODO: Create empty BM? 
  # bm <- create_empty_bm()
  # TODO: Remove looping and just use the BumpyMatrix to do all of the arithmetic.  
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
      normalized$RelativeViability <- round(trt_df$CorrectedReadout/ref_df$UntrtReadout, nDigits_rounding)
      normalized$GRvalue <- calculate_GR_value(trt_df, nDigits_rounding)
      #normalized$CorrectedReadout <- trt_df$CorrectedReadout # TODO: Within average_SE(), 

      # Normalized references.
      controls$RefRelativeViability <- round(controls$RefReadout/controls$UntrtReadout, nDigits_rounding)
      controls$RefGRvalue <- calculate_GR_value(controls)
      controls$DivisionTime <- round(rdata[i, duration_col] / log2(controls$UntrtReadout/controls$Day0Readout), nDigits_rounding)

      # TODO: Put the normalized data.frame into the bumpy matrix. 
      #bm[i, j] <- normalized
    }
  }

  # TODO: Put the bumpy matrix assay back into the SE as a new "Normalized" assay.
  return(se)
}
