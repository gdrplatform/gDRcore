#' Actual fitting function
#'
#' \code{ICGRfits} returns fit parameters
#'
#' returns fit parameters
#'
#' @import reshape2
#' @param log10concs concentrations
#' @param RelativeViability values
#' @param GRvalues values
#' @param e_0 = 1 by default
#' @param GR_0 = 1 by default
#' @param force use signifcance or not
#' @param cap enforce e_0 and GR_0
#' @return vector of parameters
#' @examples
#' sum(1:10)
#' @importFrom drc drm drmc LL.3u LL.4
#' @export
ICGRfits <- function(df_,
                     e_0 = 1,
                     GR_0 = 1,
                     force = FALSE) {
  if (sum(!is.na(df_$RelativeViability))<5) {
    df_IC = NA
  } else {
    df_IC = gDR::logisticFit(
      df_$Concentration,
      df_$RelativeViability,
      x_0 = e_0,
      curve_type = "IC",
      force = force
    )
  }
  if (sum(!is.na(df_$GRvalue))<5) {
    df_GR = NA
  } else {
    df_GR = gDR::logisticFit(
      df_$Concentration,
      df_$GRvalue,
      x_0 = GR_0,
      curve_type = "GR",
      force = force
    )
  }

  df_metrics = rbind(df_IC, df_GR)
  rownames(df_metrics) <- c("IC", "GR")

  return(df_metrics)
}




#' Actual fitting function
#'
#' \code{logisticFit} returns fit parameters
#'
#' returns fit parameters
#'
#' @import reshape2
#' @param log10concs log10 of concentrations
#' @param normValues values
#' @param x_0 upper limit; =1 by default )
#' @param curve_type response curve: either IC ([0,1]) or GR([-1,1])
#' @param force use signifcance or not
#' @param cap enforce x_0
#' @return vector of parameters
#' @examples
#' sum(1:10)
#' @importFrom drc drm drmc LL.3u
#' @export
logisticFit <-
  function(concs,
           normValues,
           x_0 = 1,
           curve_type = c("IC", "GR"),
           force = FALSE,
           cap = 0.1) {
    # Implementation of the genedata approach for curve fit: https://screener.genedata.com/documentation/display/DOC15/Business+Rules+for+Dose-Response+Curve+Fitting+Model+Selection+and+Fit+Validity
    #

    # define variables and prepare data
    log10concs <- log10(concs)
    df_ <- data.frame(log10conc = log10concs,
                     normValues = pmin(normValues, normValues + cap))

    fit_para <- c("h", "x_inf", "c50")

    out <- array(NA, length(get_header("response_metrics")))
    names(out) <- get_header("response_metrics")
    out["maxlog10Concentration"] <- max(log10concs)
    out["N_conc"] <- length(unique(log10concs))

    # fit parameters and boundaries
    if (curve_type == "IC") {
      priors <- c(2, 0.4, median(concs))
      lower <- c(.1, 0, min(concs) / 10)
    } else if (curve_type == "GR") {
      priors <- c(2, 0.1, median(concs))
      lower <- c(.1,-1, min(concs) / 10)
    }

    controls <- drc::drmc()
    controls$relTol <- 1e-06
    controls$errorm <- FALSE
    controls$noMessage <- TRUE
    controls$rmNA <- TRUE

    ######################################
    # IC curve fitting
    if (!is.na(x_0)) {
      out["x_0"] <- x_0
      upper <- c(5, min(x_0 + .1, 1), max(concs) * 10)

      output_model_new = try(drc::drm(
        normValues ~ log10conc,
        data = df_,
        logDose = 10,
        fct = drc::LL.3u(upper = x_0, names = fit_para),
        start = priors,
        lowerl = lower,
        upperl = upper,
        control = controls,
        na.action = na.omit
      ))
    } else {
      fit_para <- c("h", "x_inf", 'x_0', "c50")
      output_model_new = try(drc::drm(
        normValues ~ log10conc,
        data = df_,
        logDose = 10,
        fct = drc::LL.4(names = fit_para),
        start = c(priors[1:2], 1, priors[3]),
        lowerl = lower[c(1:2, 2, 3)],
        upperl = c(5, 1, 1, max(concs) * 10),
        control = controls,
        na.action = na.omit
      ))

    }
    # assuming proper fit result
    if (class(output_model_new) != "try-error") {
      for (p in fit_para) {
        out[p] = stats::coef(output_model_new)[paste0(p, ":(Intercept)")]
      }
      # F-test for the significance of the sigmoidal fit
      Npara <- 3 + (is.na(x_0)*1) # N of parameters in the growth curve
      Npara_flat <- 1 # F-test for the models
      RSS2 <- sum(stats::residuals(output_model_new) ^ 2, na.rm = TRUE)
      RSS1 <- sum((df_$normValues - mean(df_$normValues,
                                        na.rm = TRUE)) ^ 2, na.rm = TRUE)
      df1 <- (Npara - Npara_flat)
      df2 <- (length(na.omit(df_$normValues)) - Npara + 1)
      f_value <- ((RSS1 - RSS2) / df1) / (RSS2 / df2)
      f_pval <- stats::pf(f_value, df1, df2, lower.tail = FALSE)
      out["r2"] <- 1 - RSS2 / RSS1
    }


    # non-fitted metrics
    xAvg <- aggregate(
      df_$normValues,
      by = list(log10conc = df_$log10conc),
      FUN = function(x)
        mean(x, na.rm = T)
    )
    colnames(xAvg)[2] <- "normValues"
    l <- dim(xAvg)[1]

    out["x_max"] <- min(xAvg$normValues[c(l, l - 1)], na.rm = TRUE)

    out["x_mean"] <- mean(xAvg$normValues)
    out["x_AOC"] <- 1 - mean(xAvg$normValues)

    # analytical solution for ic50
    out["xc50"] <- out["c50"] * ((x_0 - out["x_inf"]) / (0.5 - out["x_inf"]) - 1) ^
      (1 / out["h"])

    # testing the significance of the fit and replacing with flat function if required
    pcutoff <- ifelse(force, 1, .05)
    if (!is.na(f_pval)) {
      out["flat_fit"] <- ifelse(f_pval >= pcutoff |
                                 is.na(out["c50"]), 1, 0)
    } else {
      out["flat_fit"] <- ifelse(is.na(out["c50"]), 1, 0)
    }

    # Replace values for flat fits: c50 = 0, h = 0.01 and xc50 = +/- Inf
    if (out["flat_fit"] == 1) {
      out["c50"] <- 0
      out["h"] <- 0.0001
      out["xc50"] <- ifelse(mean(xAvg$normValues) > .5, Inf,-Inf)
      out["x_inf"] <- mean(xAvg$normValues)
    }

    # Add xc50 = +/-Inf for any curves that don"t reach RelativeViability = 0.5
    if (is.na(out["xc50"])) {
      out["xc50"] <- ifelse(out["x_inf"] > .5, Inf,-Inf)
    }
    return(out)
  }


#' Actual fitting function
#'
#' \code{ICGRlogisticFit} returns fit parameters
#'
#' returns fit parameters
#'
#' @import reshape2
#' @param log10concs concentrations
#' @param RelativeViability values
#' @param GRvalues values
#' @param e_0 =1 by default
#' @param GR_0 =1 by default
#' @param force use signifcance or not
#' @param cap enforce e_0 and GR_0
#' @return vector of parameters
#' @examples
#' sum(1:10)
#' @importFrom drc drm drmc LL.3u
#' @export
ICGRlogisticFit <-
  function(log10concs,
           RelativeViability,
           GRvalues,
           e_0 = 1,
           GR_0 = 1,
           force = FALSE,
           cap = FALSE) {
    # Implementation of the genedata approach for curve fit: https://screener.genedata.com/documentation/display/DOC15/Business+Rules+for+Dose-Response+Curve+Fitting+Model+Selection+and+Fit+Validity
    #

    # define variables and prepare data
    IC_data_exp <-
      data.frame(log10conc = log10concs, RelativeViability = RelativeViability)
    GR_data_exp <-
      data.frame(log10conc = log10concs, GRvalue = GRvalues)
    concs <- 10 ** log10concs
    ICfit_parameters <- c("h_ic", "e_inf", "ec50")
    GRfit_parameters <- c("h_GR", "GRinf", "GEC50")

    out <- array(NA, length(get_header("metrics_results")))
    names(out) <- get_header("metrics_results")
    out["maxlog10Concentration"] <- max(log10concs)
    out["N_conc"] <- length(unique(log10concs))
    out["e_0"] <- e_0
    out["GR_0"] <- GR_0

    # fit parameters and boundaries
    ICpriors <- c(2, 0.4, median(concs))
    GRpriors <- c(2, 0.1, median(concs))
    IClower <- c(.1, 0, min(concs) / 10)
    GRlower <- c(.1, -1, min(concs) / 10)
    upper <- c(5, 1, max(concs) * 10)

    controls <- drc::drmc()
    controls$relTol <- 1e-06
    controls$errorm <- FALSE
    controls$noMessage <- TRUE
    controls$rmNA <- TRUE

    ######################################
    # IC curve fitting
    output_model_new <- try(drc::drm(
      RelativeViability ~ log10conc,
      data = IC_data_exp,
      logDose = 10,
      fct = drc::LL.3u(upper = e_0, names = ICfit_parameters),
      start = ICpriors,
      lowerl = IClower,
      upperl = upper,
      control = controls,
      na.action = na.omit
    ))

    # assuming proper fit result
    if (class(output_model_new) != "try-error") {
      for (p in ICfit_parameters) {
        out[p] <- stats::coef(output_model_new)[paste0(p, ":(Intercept)")]
      }
      # F-test for the significance of the sigmoidal fit
      Npara <- 3 # N of parameters in the growth curve
      Npara_flat <- 1 # F-test for the models
      RSS2 <-
        sum(stats::residuals(output_model_new) ^ 2, na.rm = TRUE)
      RSS1 <-
        sum((
          IC_data_exp$RelativeViability - mean(IC_data_exp$RelativeViability,
                                               na.rm = TRUE)
        ) ^ 2, na.rm = TRUE)
      df1 <- (Npara - Npara_flat)
      df2 <-
        (length(na.omit(IC_data_exp$RelativeViability)) - Npara + 1)
      f_value <- ((RSS1 - RSS2) / df1) / (RSS2 / df2)
      f_pval <- stats::pf(f_value, df1, df2, lower.tail = FALSE)
      out["ic_r2"] <- 1 - RSS2 / RSS1
    }


    # non-fitted metrics
    ICavg <-
      aggregate(
        IC_data_exp$RelativeViability,
        by = list(log10conc = IC_data_exp$log10conc),
        FUN = mean
      )
    colnames(ICavg)[2] <- "RelativeViability"
    l <- dim(ICavg)[1]

    out["e_max"] <-
      min(ICavg$RelativeViability[c(l, l - 1)], na.rm = TRUE)

    out["mean_viability"] <- mean(ICavg$RelativeViability)

    # analytical solution for ic50
    out["ic50"] <-
      out["ec50"] * ((e_0 - out["e_inf"]) / (0.5 - out["e_inf"]) - 1) ^ (1 /
                                                                           out["h_ic"])

    # testing the significance of the fit and replacing with flat function if required
    pcutoff <- ifelse(force, 1, .05)
    if (!is.na(f_pval)) {
      out["flat_fit_ic"] <- ifelse(f_pval >= pcutoff |
                                     is.na(out["ec50"]), 1, 0)
    } else {
      out["flat_fit_ic"] <- ifelse(is.na(out["ec50"]), 1, 0)
    }

    # Replace values for flat fits: ec50 = 0, h_ic = 0.01 and ic50 = +/- Inf
    if (out["flat_fit_ic"] == 1) {
      out["ec50"] <- 0
      out["h_ic"] <- 0.0001
      out["ic50"] <-
        ifelse(mean(ICavg$RelativeViability) > .5, Inf, -Inf)
      out["e_inf"] <- mean(ICavg$RelativeViability)
    }

    # Add ic50 = +/-Inf for any curves that don"t reach RelativeViability = 0.5
    if (is.na(out["ic50"])) {
      out["ic50"] <- ifelse(out["e_inf"] > .5, Inf, -Inf)
    }




    ######################################
    # GR curve fitting
    output_model_new <- try(drc::drm(
      GRvalue ~ log10conc,
      data = GR_data_exp,
      logDose = 10,
      fct = drc::LL.3u(upper = GR_0, names = GRfit_parameters),
      start = GRpriors,
      lowerl = GRlower,
      upperl = upper,
      control = controls,
      na.action = na.omit
    ))

    # assuming proper fit result
    if (class(output_model_new) != "try-error") {
      for (p in GRfit_parameters) {
        out[p] <- stats::coef(output_model_new)[paste0(p, ":(Intercept)")]
      }
      # F-test for the significance of the sigmoidal fit
      Npara <- 3 # N of parameters in the growth curve
      Npara_flat <- 1 # F-test for the models
      RSS2 <-
        sum(stats::residuals(output_model_new) ^ 2, na.rm = TRUE)
      RSS1 <-
        sum((
          GR_data_exp$GRvalue - mean(GR_data_exp$GRvalue, na.rm = TRUE)
        ) ^ 2,
        na.rm = TRUE)
      df1 <- (Npara - Npara_flat)
      df2 <- (length(na.omit(GR_data_exp$GRvalue)) - Npara + 1)
      f_value <- ((RSS1 - RSS2) / df1) / (RSS2 / df2)
      f_pval <- stats::pf(f_value, df1, df2, lower.tail = FALSE)
      out["GR_r2"] <- 1 - RSS2 / RSS1
    }

    # non-fitted metrics
    GRavg <-
      aggregate(
        GR_data_exp$GRvalue,
        by = list(log10conc = GR_data_exp$log10conc),
        FUN = mean
      )
    colnames(GRavg)[2] <- "GRvalue"
    l <- dim(GRavg)[1]

    out["GRmax"] <- min(GRavg$GRvalue[c(l, l - 1)], na.rm = TRUE)

    out["GR_AOC"] <-
      mean(1 - GRavg$GRvalue) # use mean for consistency with mean.viability
    # ## Alternative version: Trapezoid rule for integration of GR_AOC
    # diff_vector <- diff(GRavg$log10conc, lag = 1)
    # conc_range <- GRavg$log10conc[l] - GRavg$log10conc[1]
    # out["GR_AOC"] <- sum((1 - (GRavg$GRvalue[1:(l-1)]+GRavg$GRvalue[2:l])/2)*
    #                diff_vector, na.rm = TRUE)/conc_range

    # analytical solution for GR50
    out["GR50"] <-
      out["GEC50"] * ((GR_0 - out["GRinf"]) / (0.5 - out["GRinf"]) - 1) ^ (1 /
                                                                             out["h_GR"])

    # testing the significance of the fit and replacing with flat function if required
    pcutoff <- ifelse(force, 1, .05)
    if (!is.na(f_pval)) {
      out["flat_fit_GR"] <- ifelse(f_pval >= pcutoff |
                                     is.na(out["GEC50"]), 1, 0)
    } else {
      out["flat_fit_GR"] <- ifelse(is.na(out["GEC50"]), 1, 0)
    }

    # Replace values for flat fits: GEC50 = 0, h_GR = 0.01 and GR50 = +/- Inf
    if (out["flat_fit_GR"] == 1) {
      out["GEC50"] <- 0
      out["h_GR"] <- 0.0001
      out["GR50"] <- ifelse(mean(GRavg$GRvalue) > .5, Inf, -Inf)
      out["GRinf"] <- mean(GRavg$GRvalue)
    }

    # Add GR50 = +/-Inf for any curves that don't reach GR = 0.5
    if (is.na(out["GR50"])) {
      out["GR50"] <- ifelse(out["GRinf"] > .5, Inf, -Inf)
    }

    out
  }


range_mv <- function(c50, x_inf, hillCoefficient, x_0 = 1, linear_conc_range_uM = c(5e-3, 5)) {
   # c50 in uM
   # linear_conc_range_uM in uM (from 5nM to 5uM by default)

  conc_values = 10 ** seq(log10(linear_conc_range_uM[1]), log10(linear_conc_range_uM[2]), .01)

  viab = sapply(conc_values, function(x) logistic_4parameters(x, x_inf, x_0, c50, hillCoefficient))

  RV = (viab + 100)/100
  MV = mean(RV)

  return(MV)
}


# logistic function (not used in the file but useful for plotting externally)
#' @export
logistic_4parameters <- function(c, Vinf, V0, EC50, h) {
  Vinf + (V0 - Vinf) / (1 + (c / EC50) ^ h)
}
