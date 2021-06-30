#' Calculate a GR value.
#'
#' Calculate a GR value for a given set of dose response values.
#'
#' @param rel_viability numeric vector representing the Relative Viability.
#' @param corrected_readout numeric vector containing the corrected readout.
#' @param day0_readout numeric vector containing the day 0 readout.
#' @param untrt_readout numeric vector containing the untreated readout.
#' @param ndigit_rounding integer specifying the number of digits to use for calculation rounding.
#' @param duration numeric value specifying the length of time the cells were treated (in hours).
#' @param ref_div_time numeric value specifying the reference division time for the cell line in the experiment.
#' @param cap numeric value representing the value to cap the highest allowed relative viability at.
#' @param cl_name character string specifying the name for the cell line of interest.
#'
#' @return numeric vector containing GR values, one value for each element of the input vectors.
#' 
#' @details Note that this function expects that all vectorized numeric vectors should be of the same length. 
#' \code{calculate_GR_value} will try to greedily calculate a GR value. If no day 0 readouts are available, 
#' the \code{duration} and \code{ref_div_time} will be used to try to back-calculate a day 0 value in order 
#' to produce a GR value.
#'
#' In the case of calculating the reference GR value from multiple reference readout values, the vectorized
#' calculation is performed and then the resulting vector should be averaged outside of this function. 
#' The \code{cl_name} is used purely for warning messages and will default to \code{"cell line"}. 
#'
#' Note that it is expected that the \code{ref_div_time} and \code{duration} are reported in the same units.
#'
#' @seealso normalize_SE2
#' @name calculate_GR_value
NULL


#' @export
#' @rdname calculate_GR_value
#'
calculate_GR_value <- function(rel_viability, 
                               corrected_readout, 
                               day0_readout, 
                               untrt_readout, 
                               ndigit_rounding, 
                               duration, 
                               ref_div_time, 
                               cap = 1.25,
                               cl_name = "cell line") {

  # Assertions.
  args_to_validate <- list(rel_viability, corrected_readout, day0_readout, untrt_readout)
  args_to_validate <- args_to_validate[!is.na(args_to_validate)]
  if (length(unique(vapply(args_to_validate, FUN = length, numeric(1)))) != 1L) {
    stop("unequal vector lengths: rel_viability, corrected_readout, day0_readout, untrt_readout")
  }

  # TODO: Is it correct to put the 'any' here?
  if (any(is.na(day0_readout))) {
    ## Back-calculate the day0_readout using the reference doubling time and the duration of treatment.
    if (is.null(ref_div_time) || is.na(ref_div_time)) {
      warning(
        sprintf("no day 0, no reference doubling time, so GR values are NA for '%s'", 
          cl_name))
      GRvalue <- rep(NA, length(rel_viability))
    } else if (ref_div_time > 1.5 * duration) {
      warning(
        sprintf(paste0("reference doubling time for '%s' is '%s', too long for GR calculation ", 
                       "with assay duration ('%s'), setting GR values to NA"), 
          cl_name, ref_div_time, duration))
      GRvalue <- rep(NA, length(rel_viability))
    } else {
      warning(
        sprintf("no day 0 data, calculating GR value based on reference doubling time for '%s'", 
          cl_name))
      GRvalue <- calculate_endpt_GR_value(rel_viability = rel_viability, 
  duration = duration, 
  ref_div_time = ref_div_time, 
  cap = cap,
  ndigit_rounding = ndigit_rounding)
    }
  } else {
    GRvalue <- calculate_time_dep_GR_value(corrected_readout, 
      day0_readout,  
      untrt_readout, 
      ndigit_rounding)
  }
  GRvalue
}


#' @rdname calculate_GR_value
#' @export
#'
calculate_time_dep_GR_value <- function(corrected_readout, 
                                        day0_readout, 
                                        untrt_readout, 
                                        ndigit_rounding) {
  round(2 ^ (log2(corrected_readout / day0_readout) / log2(untrt_readout / day0_readout)), ndigit_rounding) - 1
}


#' @rdname calculate_GR_value
#' @export
#'
calculate_endpt_GR_value <- function(rel_viability, 
                                     duration, 
                                     ref_div_time, 
                                     cap = 1.25,
                                     ndigit_rounding) {
  round(2 ^ (1 + (log2(pmin(cap, rel_viability)) / (duration / ref_div_time))), ndigit_rounding) - 1
}
