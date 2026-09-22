####
# Re-exports of the generic fitting layer
####
#
# apply_fit(), the fit profile registry and the reference fit function moved to
# gDRutils, which owns the fit configuration and the fit math they are built on.
# They are re-exported here for one release cycle so that existing code — including
# the Incucyte report template in gDR, which calls get_fit_profile() unqualified —
# keeps working. Call them as gDRutils::<fn>() in new code; these re-exports go away
# in the next cycle.

#' @importFrom gDRutils apply_fit
#' @export
gDRutils::apply_fit

#' @importFrom gDRutils apply_fit_to_se
#' @export
gDRutils::apply_fit_to_se

#' @importFrom gDRutils apply_fits
#' @export
gDRutils::apply_fits

#' @importFrom gDRutils fit_drug_response_metrics
#' @export
gDRutils::fit_drug_response_metrics

#' @importFrom gDRutils fit_drug_response_metrics_4p
#' @export
gDRutils::fit_drug_response_metrics_4p

#' @importFrom gDRutils get_fit_profile
#' @export
gDRutils::get_fit_profile

#' @importFrom gDRutils get_fit_profiles
#' @export
gDRutils::get_fit_profiles

#' @importFrom gDRutils register_fit_profile
#' @export
gDRutils::register_fit_profile
