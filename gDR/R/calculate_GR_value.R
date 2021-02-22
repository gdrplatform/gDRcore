#' Calculate a GR value.
#'
#' Calculate a GR value for a given set of dose response values.
#'
#' @param df_ data.frame of dose response values containing columns named
#' \code{"CorrectedReadout"}, \code{"Day0Readout"}, and \code{"UntrtReadout"}.
#' @param duration numeric value specifying the length of time the cells were treated.
#' @param ref_div_time numeric value specifying the reference division time for the cell line of interest.
#' @param cl_name character string specifying the name for the cell line of interest.
#' @param ndigit_rounding integer specifying the number of digits to use for calculation rounding.
#'
#' @return 
#'
#' @export
#'
calculate_GR_value <- function(df_, ndigit_rounding, duration, ref_div_time, cl_name) {
  if (any(is.na(df_$Day0_readout))) {
    ## Back-calculate the day0_readout using the reference doubling time and the duration of treatment.
    if (is.null(ref_div_time) || is.na(ref_div_time)) {
      futile.logger::flog.warn(sprintf("no day 0, no reference doubling time, so GR values are NA for cell line '%s'", cl_name))
    } else if (ref_div_time > 1.5 * duration) {
      futile.logger::flog.warn(
        sprintf("reference doubling time for cell line '%s'='%s', too long for GR calculation with assay duration ('%s'), setting GR values to NA", 
          cl_name, ref_div_time, duration))
    } else {
      futile.logger::flog.warn(sprintf("no day 0 data, calculating GR value based on reference doubling time for '%s'", cl_name))
    }

    GRvalue <- calculate_endpt_GR_value(df_, duration, ref_div_time, ndigit_rounding = ndigit_rounding)
  } else {
    GRvalue <- calculate_time_dep_GR_value(df_$CorrectedReadout, df_$Day0Readout, df_$UntrtReadout, ndigit_rounding)
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
