#' Calculate a GR value.
#'
#' Calculate a GR value for a given set of dose response values.
#'
#' @param df_ data.frame of dose response values containing \code{"CorrectedReadout"}.
#' @param cellline_md named list of metadata for the cell line of interest.
#' @param day0_readout numeric vector 
#' @param ndigit_rounding integer
#'
#' @return 
#'
#' @export
#'
calculate_GR_value <- function(df_, cellline_md, day0_readout, ndigit_rounding) {
  if (any(is.na(day0_readout))) {
    ## Back-calculate the day0_readout using the reference doubling time and the duration of treatment.
    # TODO: These column values should later be read directly from the SE itself. 
    ref_div_time <- as.numeric(cellline_md[[gDRutils::get_identifier("cellline_ref_div_time")]])
    cl_name <- cellline_md[[gDRutils::get_identifier("cellline_name")]]
    duration <- cellline_md[[gDRutils::get_identifier("duration")]]

    if (is.null(ref_div_time) || is.na(ref_div_time)) {
      futile.logger::flog.warn(paste("no day 0, no reference doubling time, so GR values are NA for cell line", cl_name))
    } else if (ref_div_time > 1.5 * duration) {
      futile.logger::flog.warn(sprintf("reference doubling time for cell line '%s'='%s', too long for GR calculation with assay duration ('%s'), setting GR values to NA", 
        cl_name, ref_div_time, duration))
    } else {
      futile.logger::flog.warn(paste("no day 0, calculating GR value based on reference doubling time for ", cl_name))
    }

    GRvalue <- calculate_endpt_GR_value(df_, duration, ref_div_time, ndigit_rounding = nDigits_rounding)
  } else {
    GRvalue <- calculate_time_dep_GR_value(df_$CorrectedReadout, df_$Day0Readout, df_$UntrtReadout, nDigits_rounding)
  }
  GRvalue
}


#' @export
#'
calculate_time_dep_GR_value <- function(ref_readout, day0_readout, untrt_readout, ndigit_rounding) {
  round(2 ^ (log2(ref_readout/day0_readout) / log2(untrt_readout/day0_readout)), ndigit_rounding) - 1
}


#' @export
#'
calculate_endpt_GR_value <- function(rel_viability, duration, ref_div_time, cap = 1.25, ndigit_rounding) {
  round(2 ^ (1 + (log2(pmin(cap, rel_viability)) / (duration/ref_div_time))), ndigit_rounding) - 1
}
