#' Fit curves
#'
#' Fit GR and RV curves from a data.table.
#'
#' @param df_ data.table containing data to fit. See details.
#' @param series_identifiers character vector of the column names in \code{data.table}
#' whose combination represents a unique series for which to fit curves.
#' @param e_0 numeric value representing the \code{x_0} value for the RV curve.
#' Defaults to \code{1}.
#' @param GR_0 numeric value representing the \code{x_0} value for the GR curve.
#' Defaults to \code{1}.
#' @param n_point_cutoff integer of how many points should be considered the minimum required to try to fit a curve.
#' Defaults to \code{4}.
#' @param range_conc numeric vector of length 2 indicating the lower and upper concentration ranges.
#' Defaults to \code{c(5e-3, 5)}. See details.
#' @param force_fit boolean indicating whether or not to force a constant fit.
#' Defaults to \code{FALSE}.
#' @param pcutoff numeric of pvalue significance threshold above or equal to which to use a constant fit.
#' Defaults to \code{0.05}.
#' @param cap numeric value capping \code{norm_values} to stay below (\code{x_0} + cap).
#' Defaults to \code{0.1}.
#' @param normalization_type character vector of types of curves to fit.
#' Defaults to \code{c("GR", "RV")}.
#' @keywords fit_curves
#'
#' @return data.table of fit parameters as specified by the \code{normalization_type}.
#'
#' @details
#' The \code{df_} expects the following columns:
#'  \itemize{
#'   \item{RelativeViability} normalized relative viability values (if \code{normalization_type} includes \code{"RV"})
#'   \item{GRvalue} normalized GR values (if \code{normalization_type} includes \code{"GR"})
#'  }
#'
#' The \code{range_conc} is used to calculate the \code{x_AOC_range} statistic.
#' The purpose of this statistic is to enable comparison across different experiments with slightly
#' different concentration ranges.
#'
#' @examples 
#' df_ <- data.table::data.table(Concentration = c(0.001, 0.00316227766016838, 
#' 0.01, 0.0316227766016838),
#' x_std = c(0.1, 0.1, 0.1, 0.1), normalization_types = c("RV", "RV", "RV", "RV"),
#' x = c(0.9999964000144, 0.999964001439942, 0.999640143942423, 0.996414342629482))
#' 
#' fit_curves(df_, "Concentration", normalization_type = "RV")
#'
#' @export
#'
fit_curves <- function(df_,
                       series_identifiers,
                       e_0 = 1,
                       GR_0 = 1,
                       n_point_cutoff = 4,
                       range_conc = c(5e-3, 5),
                       force_fit = FALSE,
                       pcutoff = 0.05,
                       cap = 0.1,
                       normalization_type = c("GR", "RV")) {
  if (length(series_identifiers) != 1L) {
    stop("gDR does not yet support multiple series_identifiers, feature coming soon")
  }
  stopifnot(any(inherits(df_, "data.table"), inherits(df_, "DFrame")))
  if (any(bad_normalization_type <- ! normalization_type %in% c("GR", "RV"))) {
    stop(sprintf("unknown curve type: '%s'",
                 paste0(normalization_type[bad_normalization_type], collapse = ", ")))
  }
  
  req_fields <- series_identifiers
  opt_fields <- NULL
  
  req_fields <- c(req_fields, "x")
  opt_fields <- "x_std"
  
  if (!all(req_fields %in% colnames(df_))) {
    stop(sprintf("missing one of the required fields: '%s'", paste(req_fields, collapse = ", ")))
  }
  
  if (length(setdiff(opt_fields, colnames(df_))) > 0L) {
    df_[, setdiff(opt_fields, colnames(df_))] <- NA
  }
  df_metrics <- .applyLogisticFit(df_, normalization_type, series_identifiers, e_0, GR_0, range_conc, force_fit, 
                                  pcutoff, cap, n_point_cutoff)
  
  is_unique_normalization_type_and_fit_source <- 
    nrow(unique(df_metrics[, c("normalization_type", "fit_source")])) == nrow(df_metrics)
  if (!is_unique_normalization_type_and_fit_source) {
    stop("'normalization_type' and 'fit_source' columns do not create unique combinations") 
  }
  rownames(df_metrics) <- paste0(df_metrics$normalization_type, "_", df_metrics$fit_source)
  
  concsNA <- all(is.na(unique(df_[[series_identifiers]])))
  if (concsNA) df_metrics[] <- NA
  df_metrics
}

#' @keywords internal
.applyLogisticFit <- function(df_, normalization_type, series_identifiers, e_0, GR_0, range_conc, force_fit, 
                              pcutoff, cap, n_point_cutoff) {
  
  df_metrics <- NULL
  concs <- unique(df_[[series_identifiers]])
  med_concs <- stats::median(concs)
  min_concs <- min(concs)
  
  concsNA <- all(is.na(concs))
  if (concsNA) concs[] <- 0
  
  if ("RV" %in% normalization_type) {
    df_metrics <- logisticFit(
      concs,
      df_$x[df_$normalization_type == "RV"],
      df_$x_std[df_$normalization_type == "RV"],
      priors = c(2, 0.4, 1, med_concs),
      lower = c(0.1, 0, 0, min_concs / 10),
      x_0 = e_0,
      range_conc = range_conc,
      force_fit = force_fit,
      pcutoff = pcutoff,
      cap = cap,
      n_point_cutoff = n_point_cutoff
    )
    df_metrics$normalization_type <- "RV"
  }
  
  if ("GR" %in% normalization_type) {
    df_gr <- logisticFit(
      concs,
      df_$x[df_$normalization_type == "GR"],
      df_$x_std[df_$normalization_type == "GR"],
      priors = c(2, 0.1, 1, med_concs),
      lower = c(0.1, -1, -1, min_concs / 10),
      x_0 = GR_0,
      range_conc = range_conc,
      force_fit = force_fit,
      pcutoff = pcutoff,
      cap = cap,
      n_point_cutoff = n_point_cutoff
    )
    df_gr$normalization_type <- "GR"
    df_metrics <- data.table::rbindlist(list(df_metrics, df_gr), fill = TRUE)
  }
  
  df_metrics$fit_source <- "gDR"
  df_metrics
}


#' Logistic fit
#'
#' Fit a logistic curve to drug response data.
#'
#' @param concs concentrations that have not been transformed into log space.
#' @param norm_values normalized response values (Untreated = 1).
#' @param std_norm_values std of values.
#' @param x_0 upper limit.
#' Defaults to \code{1}. For co-treatments, this value should be set to \code{NA}.
#' @param priors numeric vector containing starting values for all.
#' mean parameters in the model. Overrules any self starter function.
#' @param lower numeric vector of lower limits for all parameters in a 4-param model.
#' @param range_conc range of concentration for calculating AOC_range.
#' @param force_fit boolean indicating whether or not to force a parameter-based fit.
#' @param pcutoff numeric of pvalue significance threshold above or equal to which to use a constant fit.
#' @param cap numeric value capping \code{norm_values} to stay below (\code{x_0} + cap).
#' @param n_point_cutoff integer indicating number of unique concentrations required to fit curve.
#' @param capping_fold Integer value of the fold number to use for capping IC50/GR50. Default is \code{5}.
#'
#' @examples
#' logisticFit(
#' c(0.001, 0.00316227766016838, 0.01, 0.0316227766016838),
#' c(0.9999964000144, 0.999964001439942, 0.999640143942423, 0.996414342629482),
#' rep(0.1, 4),
#' priors = c(2, 0.4, 1, 0.00658113883008419)
#' )
#'
#' @return data.table with metrics and fit parameters.
#'
#' @details
#' Implementation of the genedata approach for curve fit:
#' https://screener.genedata.com/documentation/display/DOC21/Business-Rules-for-Dose-Response-Curve-Fitting,-Model-Selection,-and-Fit-Validity.html #nolint
#'
#' The output parameter names correspond to the following definitions:
#' \describe{
#'  \item{x_mean}{The mean of a given dose-response metric}
#'  \item{x_AOC_range}{The range of the area over the curve}
#'  \item{x_AOC}{The area over the GR curve or, respectively, under the relative
#'  cell count curve, averaged over the range of concentration values}
#'  \item{xc50}{The concentration at which the effect reaches a value of 0.5 based
#'  on interpolation of the fitted curve}
#'  \item{x_max}{The maximum effect of the drug}
#'  \item{ec50}{The drug concentration at half-maximal effect}
#'  \item{x_inf}{The asymptotic value of the sigmoidal fit to the dose-response
#'  data as concentration goes to infinity}
#'  \item{x_0}{The asymptotic metric value corresponding to a concentration of 0
#'  for the primary drug}
#'  \item{h}{The hill coefficient of the fitted curve, which reflects how steep
#'  the dose-response curve is}
#'  \item{r2}{The goodness of the fit}
#'  \item{x_sd_avg}{The standard deviation of GR/IC}
#'  \item{fit_type}{This will be given by one of the following:
#'   \itemize{
#'    \item{"DRC4pHillFitModel"} Successfully fit with a 4-parameter model
#'    \item{"DRC3pHillFitModelFixS0"} Successfully fit with a 3-parameter model
#'    \item{"DRCConstantFitResult"} Successfully fit with a constant fit
#'    \item{"DRCTooFewPointsToFit"} Not enough points to run a fit
#'    \item{"DRCInvalidFitResult"} Fit was attempted but failed
#'   }
#'  }
#'  \item{maxlog10Concentration}{The highest log10 concentration}
#'  \item{N_conc}{Number of unique concentrations}
#' }
#' @keywords fit_curves
#'
#' @export
#'
logisticFit <-
  function(concs,
           norm_values,
           std_norm_values = NA,
           x_0 = 1,
           priors = NULL,
           lower = NULL,
           range_conc = c(5e-3, 5),
           force_fit = FALSE,
           pcutoff = 0.05,
           cap = 0.1,
           n_point_cutoff = 4,
           capping_fold = 5) {
    if (length(concs) != length(norm_values)) {
      stop("unequal vector lengths for 'conc' and 'norm_values'")
    }
    # Check that values have not been logged yet. 
    if (any(concs < 0)) {
      stop("logisticFit accepts only unlogged concentrations, negative concentrations are detected")
    }
    
    out <- .setup_metric_output()
    out$maxlog10Concentration <- log10(max(concs))
    out$N_conc <- length(unique(concs))
    out$x_sd_avg <- mean(std_norm_values, na.rm = TRUE)
    # Cap norm_values at (x_0 + cap) so as not to throw off the fit.
    limit <- if (is.na(x_0)) {
      1 + cap
    } else {
      x_0 + cap
    }
    
    norm_values <- pmin(norm_values, limit)
    df_ <- data.table::data.table(concs = concs, norm_values = norm_values)
    if (has_dups(df_$concs)) {
      warning("duplicates were found, averaging values")
      df_ <- average_dups(df_, "concs")
    }
    
    mean_norm_value <- mean(df_$norm_values, na.rm = TRUE)
    out$x_mean <- mean_norm_value
    out$x_AOC <- .calculate_complement(mean_norm_value)
    out$x_max <- .calculate_x_max(df_)
    ## Perform a 3-param or 4-param fit.
    ## Fit type is determined based on number of free variables available.
    fit_param <- c("h", "x_inf", "x_0", "ec50")
    controls <- drc::drmc(relTol = 1e-06, errorm = FALSE, noMessage = TRUE, rmNA = TRUE)
    
    out <-
      .setLogisticFit(out = out, df_ = df_, n_point_cutoff = n_point_cutoff, fit_param = fit_param,
                      priors = priors, lower = lower, force_fit = force_fit, x_0 = x_0, cap = cap, 
                      concs = concs, controls = controls, range_conc = range_conc, pcutoff = pcutoff, 
                      capping_fold = capping_fold, mean_norm_value = mean_norm_value)
    
    data.table::setDT(out)
    out
  }

#' @keywords internal
.setLogisticFit <- function(out, df_, n_point_cutoff, fit_param, priors, lower, force_fit, x_0, cap, concs, controls,
                            range_conc, pcutoff, capping_fold, mean_norm_value) {
  out <- tryCatch({
    .checkNonNaAvgNorm(df_ = df_, n_point_cutoff = n_point_cutoff)
    fit_model <- .prepareFitModel(out, df_, x_0, fit_param, priors, lower, cap, concs, controls)
    out$fit_type <- "DRC4pHillFitModel"
    if (!is.na(x_0)) {
      fit_param <- fit_param[-3]
      priors <- priors[-3]
      lower <- lower[-3]
      out$fit_type <- "DRC3pHillFitModelFixS0"
      out$x_0 <- x_0
    }
    out <- .set_model_fit_params(out, fit_model, fit_param)
    out$x_mean <- .predict_mean_from_model(fit_model, min(df_$concs), max(df_$concs))
    out$x_AOC <- 1 - out$x_mean
    out$x_AOC_range <- 1 - .predict_mean_from_model(fit_model, range_conc[1], range_conc[2])
    # F-test for the significance of the sigmoidal fit.
    out$rss <- RSS2 <- sum(stats::residuals(fit_model) ^ 2, na.rm = TRUE)
    RSS1 <- sum((df_$norm_values - mean(df_$norm_values, na.rm = TRUE)) ^ 2, na.rm = TRUE)
    out$r2 <- 1 - RSS2 / RSS1
    out$xc50 <- .calculate_xc50(ec50 = out$ec50, x0 = out$x_0, xInf = out$x_inf, h = out$h)
    # Test the significance of the fit and replace with flat function if required.
    nparam <- 3 + (is.na(x_0) * 1) # N of parameters in the growth curve; if (x0 = NA) {4}
    df1 <- nparam - 1 # (N of parameters in the growth curve) - (F-test for the models)
    df2 <- length(stats::na.omit(df_$norm_values)) - nparam + 1
    out$p_value <- f_pval <- .calculate_f_pval(df1, df2, RSS1, RSS2)
    if (all((!force_fit), 
           any(all(exists("f_pval"), !is.na(f_pval), f_pval >= pcutoff), is.na(out$ec50)))) {
      stop(fitting_handler(
        "constant_fit",
        message = sprintf("fit is not statistically significant (p=%.2f), setting constant fit", f_pval)
      ))
    }
    .extendWithXc50(x = out, concs = concs, capping_fold = capping_fold)
  }, too_few_fit = function(e) {
    warning(e$message)
    out <- .set_too_few_fit_params(out, df_$norm_values)
  }, constant_fit = function(e) {
    if (!is.na(x_0)) {
      warning(sprintf("overriding original x_0 argument '%s' with '%s' (%s)", x_0, mean_norm_value, e$message))
    }
    out <- set_constant_fit_params(out, mean_norm_value)
  }, invalid_fit = function(e) {
    warning(sprintf("fitting failed with: '%s'", e))
    out <- .set_invalid_fit_params(out, df_$norm_values)
  }, error = function(e) {
    stop(e)
  })
  out
}

#' @keywords internal
.checkNonNaAvgNorm <- function(df_, n_point_cutoff) {
  non_na_avg_norm <- !is.na(df_$norm_values)
  if (sum(non_na_avg_norm) < n_point_cutoff) {
    stop(fitting_handler(
      "too_few_fit",
      message = sprintf(
        "not enough data points (%i < %i) to perform fitting",
        sum(non_na_avg_norm),
        n_point_cutoff
      )
    ))
  }
  if (length(unique(df_$norm_values[non_na_avg_norm])) == 1L) {
    stop(fitting_handler("constant_fit", message = "only 1 normalized value detected, setting constant fit"))
  }
}

#' @keywords internal
.extendWithXc50 <- function(x, concs, capping_fold) {
  # Add xc50 = +/-Inf for any curves that do not reach RV/GR = 0.5.
  if (is.na(x$xc50)) {
    x$xc50 <- .estimate_xc50(x$x_inf)
  } else {
    # set the xc50 to Inf if the value is extrapolated beyond to 5-fold above/below the 
    # max/min tested concentrations (default)
    x$xc50 <- cap_xc50(
      x$xc50, 
      max_conc = 10 ^ x$maxlog10Concentration, 
      min_conc = min(concs[concs > 0]), 
      capping_fold = capping_fold
    )
  }
  x
}

#' @keywords internal
.prepareFitModel <- function(x, df_, x_0, fit_param, priors, lower, cap, concs, controls) {
  if (!is.na(x_0)) {
    # Override existing params for removal of x_0 parameter (i.e. cotreatments).
    fit_param <- fit_param[-3]
    priors <- priors[-3]
    lower <- lower[-3]
    
    fct <- drc::LL.3u(upper = x_0, names = fit_param)
    upperl <- c(5, min(x_0 + cap, 1), max(concs) * 10)
    
    x$fit_type <- "DRC3pHillFitModelFixS0"
    x$x_0 <- x_0
  } else {
    fct <- drc::LL.4(names = fit_param)
    upperl <- c(5, 1, 1 + cap, max(concs) * 10)
    
    x$fit_type <- "DRC4pHillFitModel"
  }
  
  drc::drm(
    norm_values ~ concs,
    data = df_,
    logDose = NULL,
    fct = fct,
    start = priors,
    lowerl = lower,
    upperl = upperl,
    control = controls,
    na.action = stats::na.omit
  )
}

#' Predict efficacy values given fit parameters and a concentration.
#'
#' Predict efficacy values given fit parameters and a concentration.
#'
#' @param c Numeric vector representing concentrations to predict efficacies for.
#' @param x_inf Numeric vector representing the asymptotic value of the sigmoidal fit to the dose-response
#'  data as concentration goes to infinity.
#' @param x_0 Numeric vector representing the asymptotic metric value corresponding to a concentration of 0
#'  for the primary drug.
#' @param ec50 Numeric vector representing the drug concentration at half-maximal effect.
#' @param h Numeric vector representing the hill coefficient of the fitted curve, which reflects how steep
#'  the dose-response curve is.
#' @keywords fit_curves
#'
#' @return Numeric vector representing predicted efficacies from given concentrations and fit parameters.
#'
#' @details The inverse of this function is \code{predict_conc_from_efficacy}.
#' @seealso predict_conc_from_efficacy
#' 
#' @examples 
#' predict_efficacy_from_conc(c = 1, x_inf = 0.1, x_0 = 1, ec50 = 0.5, h = 2)
#' 
#' @export
predict_efficacy_from_conc <- function(c, x_inf, x_0, ec50, h) {
  checkmate::assert_numeric(c)
  checkmate::assert_numeric(x_inf)
  checkmate::assert_numeric(x_0)
  checkmate::assert_numeric(ec50)
  checkmate::assert_numeric(h)
  assert_equal_input_len(outlier = c, x_inf, x_0, ec50, h)

  as.numeric(ifelse(c > 0,
                    x_inf + (x_0 - x_inf) / (1 + (c / ec50) ^ h),
                    x_0)) # avoid issues with c=0 for DRCConstantFitResult
}


#' Predict a concentration for a given efficacy with fit parameters.
#'
#' Predict a concentration for a given efficacy with fit parameters.
#'
#' @details The inverse of this function is \code{predict_efficacy_from_conc}.
#' 
#' @param efficacy Numeric vector representing efficacies to predict concentrations for.
#' @param x_inf Numeric vector representing the asymptotic value of the sigmoidal fit to the dose-response
#'  data as concentration goes to infinity.
#' @param x_0 Numeric vector representing the asymptotic metric value corresponding to a concentration of 0
#'  for the primary drug.
#' @param ec50 Numeric vector representing the drug concentration at half-maximal effect.
#' @param h Numeric vector representing the hill coefficient of the fitted curve, which reflects how steep
#' @keywords fit_curves
#'
#' @return Numeric vector representing predicted concentrations from given efficacies and fit parameters.
#'
#' @examples 
#' predict_conc_from_efficacy(efficacy = c(1, 1.5), x_inf = 0.1, x_0 = 1, ec50 = 0.5, h = 2)
#' 
#' @seealso predict_efficacy_from_conc .calculate_x50
#' @export
predict_conc_from_efficacy <- function(efficacy, x_inf, x_0, ec50, h) {
  assert_equal_input_len(outlier = efficacy, x_inf, x_0, ec50, h)
  ifelse(efficacy > x_0,
    0,
    ifelse(efficacy < x_inf,
      Inf,
      ec50 * ((x_0 - x_inf) / (efficacy - x_inf) - 1) ^ (1 / h)
    )
  )
}


logistic_metrics <- function(c, x_metrics) {
  metrics <- c("x_inf", "x_0", "h", "ec50")
  if (all(metrics %in% names(x_metrics))) {
    DRC_metrics <- as.vector(x_metrics[metrics])
  } else {
    gr_metric_cols <- get_header("GR_metrics")
    rv_metric_cols <- get_header("RV_metrics")
    if (all(gr_metric_cols %in% names(x_metrics))) {
      metric_cols <- gr_metric_cols
    } else if (all(rv_metric_cols %in% names(x_metrics))) {
      metric_cols <- rv_metric_cols
    } else {
      stop("wrong input parameters")
    }
    DRC_metrics <- as.vector(x_metrics[metric_cols[metrics]])
    names(DRC_metrics) <- metrics
  }

  DRC_metrics$x_inf + (DRC_metrics$x_0 - DRC_metrics$x_inf) /
                          (1 + (c / DRC_metrics$ec50) ^ DRC_metrics$h)
}


##################
# Helper functions
##################

#' @keywords fit_curves
#' @export
.setup_metric_output <- function() {
  resp_metric_all_cols <- get_header("response_metrics")
  # remove cols ending with "_sd"
  # they are not present in the primary assays 
  # but only with the assays followed by averaging of biological replicates
  resp_metric_cols <- resp_metric_all_cols[!endsWith(resp_metric_all_cols, "_sd")]
  
  out <- as.list(rep(NA, length(resp_metric_cols)))
  names(out) <- resp_metric_cols
  out
}


# Estimate values for undefined IC/GR50 values.
#' @keywords internal
.estimate_xc50 <- function(param) {
  if (all(is.na(param))) {
    NA
  } else if (all(param > 0.5, na.rm = TRUE)) {
    Inf
  } else if (all(param <= 0.5, na.rm = TRUE)) {
    -Inf
  } else {
    NA
  }
}


#' @keywords internal
has_dups <- function(vec) {
  freq <- table(vec)
  any(freq != 1L)
}


#' @keywords internal
average_dups <- function(dt, col) {
  checkmate::assert_subset(col, names(dt))
  checkmate::assert_data_table(dt, types = "numeric")
  dt[, lapply(.SD, mean, na.rm = TRUE), by = col]
}

############
# Setters
############

#' @keywords fit_curves
#' @export
.set_mean_params <- function(out, mean_norm_value) {
  checkmate::assert_number(mean_norm_value)

  out$xc50 <- .estimate_xc50(mean_norm_value)

  out$x_0 <- out$x_inf <- out$x_mean <- mean_norm_value
  out$x_AOC_range <- out$x_AOC <- .calculate_complement(mean_norm_value)
  out
}


#' @keywords internal
.set_model_fit_params <- function(out, model, fit_param) {
  for (p in fit_param) {
    # drm will output model with the ":(Intercept)" term concatenated at end.
    out[[p]] <- stats::coef(model)[[paste0(p, ":(Intercept)")]]
  }
  out
}


#' @keywords internal
.set_too_few_fit_params <- function(out, norm_values) {
  out$fit_type <- "DRCTooFewPointsToFit"
  out$xc50 <- .estimate_xc50(norm_values)
  out
}


#' Set fit parameters for a constant fit.
#'
#' Replace values for flat fits: ec50 = 0, h = 0.0001 and xc50 = +/- Inf
#'
#' @param out Named list of fit parameters.
#' @param mean_norm_value Numeric value that be used to set all parameters
#' that can be calculated from the mean.
#' @keywords fit_curves
#' @return Modified named list of fit parameters.
#' 
#' @examples 
#' na <- list(x_0 = NA)
#' set_constant_fit_params(na, mean_norm_value = 0.6)
#' 
#' @export
set_constant_fit_params <- function(out, mean_norm_value) {
  out$fit_type <- "DRCConstantFitResult"
  out$ec50 <- 0
  out$h <- 0.0001
  out$r2 <- 0

  out <- .set_mean_params(out, mean_norm_value)
  out
}


#' Set fit parameters for an invalid fit.
#' @param out Named list of fit parameters.
#' @param norm_values Numeric vector used to estimate an \code{xc50} value.
#' @keywords fit_curves
#' @return Modified named list of fit parameters.
#' 
#' @examples 
#' .set_invalid_fit_params(list(), norm_values = rep(0.3, 6))
#' 
#' @export
.set_invalid_fit_params <- function(out, norm_values) {
  out$fit_type <- "DRCInvalidFitResult"
  out$r2 <- NA
  out$xc50 <- .estimate_xc50(norm_values)
  out
}

#################
# Calculations
#################

#' @keywords internal
.predict_mean_from_model <- function(model, min, max, intervals = 100) {
  lg_min_con <- log10(min)
  lg_max_con <- log10(max)
  inputs <- data.table::data.table(concs = 10 ^ (seq(lg_min_con, lg_max_con, (lg_max_con - lg_min_con) / intervals)))
  mean(stats::predict(model, inputs), na.rm = TRUE)
}


#' @keywords internal
.calculate_xc50 <- function(ec50, x0, xInf, h) {
  predict_conc_from_efficacy(efficacy = 0.5, x_inf = xInf, x_0 = x0, ec50 = ec50, h = h)
}


# 'x_max' can be considered either the lowest readout (max efficacy)
# or the efficacy at the max concentration. We take the min
# of the two highest concentrations as a compromise.
#' @keywords internal
.calculate_x_max <- function(df) {
  stopifnot(c("concs", "norm_values") %in% colnames(df))

  # Ascending order for position norm_values subsetting.
  data.table::setorder(df, concs)
  norm_values <- df$norm_values

  if (all(is.na(norm_values))) {
    x_max <- NA
  } else {
    norm_values <- norm_values[!is.na(norm_values)]
    n <- length(norm_values)
    x_max <- min(norm_values[c((n - 1):n)])
  }
  x_max
}


#' @keywords internal
.calculate_f_pval <- function(df1, df2, RSS1, RSS2) {
  f_value <- ((RSS1 - RSS2) / df1) / (RSS2 / df2)
  f_pval <- stats::pf(f_value, df1, df2, lower.tail = FALSE)
  f_pval
}


#' @export
#' @keywords fit_curves
.calculate_complement <- function(x) {
  1 - x
}

#' Cap XC50 value.
#' 
#' Set IC50/GR50 value to \code{Inf} or \code{-Inf} based on upper and lower limits.
#'
#' @details 
#' Note: \code{xc50} and \code{max_conc} should share the same units.
#' Ideally, the \code{lower_cap} should be based on the lowest tested concentration.
#' However, since we don't record that, it is set 5 orders of magnitude below the highest dose.
#' 
#' @param xc50 Numeric value of the IC50/GR50 to cap. 
#' @param max_conc Numeric value of the highest concentration in a dose series used to calculate the \code{xc50}.
#' @param min_conc Numeric value of the lowest concentration in a dose series used to calculate the \code{xc50}. 
#' If \code{NA} (default), using \code{max_conc/1e5} instead.
#' @param capping_fold Integer value of the fold number to use for capping. Defaults to \code{5}.
#' @keywords fit_curves
#'
#' @return Capped IC50/GR50 value.
#'
#' @examples 
#' cap_xc50(xc50 = 1, max_conc = 2)
#' cap_xc50(xc50 = 2, max_conc = 5, min_conc = 1)
#' cap_xc50(xc50 = 26, max_conc = 5, capping_fold = 5)
#'
#' @export
cap_xc50 <- function(xc50, max_conc, min_conc = NA, capping_fold = 5) {
  checkmate::assert_numeric(capping_fold)
  checkmate::assert_number(xc50)
  checkmate::assert_number(max_conc)
  checkmate::assert_number(min_conc, na.ok = TRUE)
  
  upper_cap <- max_conc * capping_fold
  lower_cap <- if (!is.na(min_conc)) {
    min_conc / capping_fold
  } else {
    max_conc / (capping_fold * 1e5)
  }
  if (xc50 > upper_cap) {
    xc50 <- Inf
  } else if (xc50 < lower_cap) {
    xc50 <- -Inf
  }
  xc50
} 

#################
# Error handling
#################

#' @keywords internal
fitting_handler <- function(subclass, message, call = sys.call(-1), ...) {
  structure(
    class = c(subclass, "error", "condition"),
    list(message = message, call = call, ...)
  )
}
