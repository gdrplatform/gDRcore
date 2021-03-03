#' Calculate a GR value.
#'
#' Calculate a GR value for a given set of dose response values.
#'
#' @param rel_viability numeric vector representing the Relative Viability.
#' @param corrected_readout numeric vector containing the corrected readout.
#' @param day0_readout numeric vector containing the day 0 readout.
#' @param untrt_readout numeric vector containing the untreated readout.
#' @param ndigit_rounding integer specifying the number of digits to use for calculation rounding.
#' @param duration numeric value specifying the length of time the cells were treated.
#' @param ref_div_time numeric value specifying the reference division time for the cell line of interest.
#' @param cl_name character string specifying the name for the cell line of interest.
#'
#' @return numeric GR value.
#'
#' @export
#'
calculate_GR_value <- function(rel_viability, 
                               corrected_readout, 
                               day0_readout, 
                               untrt_readout, 
                               ndigit_rounding, 
                               duration, 
                               ref_div_time, 
                               cl_name = "cell line") {

  if (any(is.na(day0_readout))) {
    ## Back-calculate the day0_readout using the reference doubling time and the duration of treatment.
    if (is.null(ref_div_time) || is.na(ref_div_time)) {
      futile.logger::flog.warn(
        sprintf("no day 0, no reference doubling time, so GR values are NA for '%s'", 
          cl_name))
    } else if (ref_div_time > 1.5 * duration) {
      futile.logger::flog.warn(
        sprintf("reference doubling time for '%s'='%s', too long for GR calculation with assay duration ('%s'), setting GR values to NA", 
          cl_name, ref_div_time, duration))
    } else {
      futile.logger::flog.warn(
        sprintf("no day 0 data, calculating GR value based on reference doubling time for '%s'", 
          cl_name))
    }
    
    GRvalue <- calculate_endpt_GR_value(rel_viability = rel_viability, 
      duration = duration, 
      ref_div_time = ref_div_time, 
      ndigit_rounding = ndigit_rounding)
  } else {
    GRvalue <- calculate_time_dep_GR_value(corrected_readout, 
      day0_readout,  
      untrt_readout, 
      ndigit_rounding)
  }
  GRvalue
}


# TODO: Add to same documenatation family as above. 
#' @export
#'
calculate_time_dep_GR_value <- function(trt_readout, 
                                        day0_readout, 
                                        untrt_readout, 
                                        ndigit_rounding) {
  round(2 ^ (log2(trt_readout/day0_readout) / log2(untrt_readout/day0_readout)), ndigit_rounding) - 1
}


# TODO: Add to same documenatation family as above. 
#' @export
#'
calculate_endpt_GR_value <- function(rel_viability, 
                                     duration, 
                                     ref_div_time, 
                                     cap = 1.25, 
                                     ndigit_rounding) {
  round(2 ^ (1 + (log2(pmin(cap, rel_viability)) / (duration/ref_div_time))), ndigit_rounding) - 1
}
