library(reshape2)

#' Actual fitting function
#'
#' \code{GRlogisticFit} returns fit parameters
#'
#' returns fit parameters
#'
#' @param log10concs concentrations
#' @param GRvalues values
#' @param upper_GR =1 by default
#' @param force use signifcance or not
#' @param cap enforce upper_GR
#' @return vector of values
#' @examples
#' sum(1:10)
#' @export
GRlogisticFit <- function(log10concs, GRvalues, upper_GR = 1, force = FALSE, cap = FALSE) {
    # TODO: test properly and match algortihm to GENEDATA

    # Implementation of the genedata approach for curve fit: https://screener.genedata.com/documentation/display/DOC15/Business+Rules+for+Dose-Response+Curve+Fitting+Model+Selection+and+Fit+Validity
    #

    # define variables and prepare data
    data_exp = data.frame(log10conc=log10concs, GRvalue=GRvalues)
    concs = 10**log10concs
    fit_parameters = c('h_GR','GRinf','GEC50')
    metrics = c('GR50', 'GRmax', 'GR_AOC', 'R_square_GR', 'pval_GR',
        'flat_fit_GR', 'maxlog10Concentration', 'N_conc', 'log10_conc_step', 'upper_GR')
    out = array(NA, length(fit_parameters) + length(metrics))
    names(out) = c(fit_parameters, metrics)
    out['maxlog10Concentration'] = max(log10concs)
    out['upper_GR'] = upper_GR

    # fit parameters and boundaries
    priors = c(2, 0.1, median(concs))
    lower = c(.1, -1, min(concs)/10)
    upper = c(5, 1, max(concs)*10)

    controls = drc::drmc()
    controls$relTol = 1e-06
    controls$errorm = FALSE
    controls$noMessage = TRUE
    controls$rmNA = TRUE

    # GR curve fitting
    output_model_new = try(drc::drm(
            GRvalue ~ log10conc,
            data=data_exp,
            logDose = 10,
            fct=drc::LL.3u(upper = upper_GR, names = fit_parameters),
            start = priors, lowerl = lower, upperl = upper, control = controls,
            na.action = na.omit))

    # assuming proper fit result
    if(class(output_model_new)!="try-error") {
        for (p in fit_parameters) {
            out[p] = stats::coef(output_model_new)[paste0(p, ':(Intercept)')]
        }
        # F-test for the significance of the sigmoidal fit
        Npara = 3 # N of parameters in the growth curve
        Npara_flat = 1 # F-test for the models
        RSS2 = sum(stats::residuals(output_model_new)^2, na.rm = TRUE)
        RSS1 = sum((data_exp$GRvalue - mean(data_exp$GRvalue, na.rm = TRUE))^2,
                   na.rm = TRUE)
        df1 = (Npara - Npara_flat)
        df2 = (length(na.omit(data_exp$GRvalue)) - Npara + 1)
        f_value = ((RSS1-RSS2)/df1)/(RSS2/df2)
        f_pval = stats::pf(f_value, df1, df2, lower.tail = FALSE)
        out['pval_GR'] = f_pval
        out['R_square_GR'] = 1 - RSS2/RSS1
    }

    # non-fitted metrics
    GRavg = aggregate(data_exp$GRvalue, by=list(log10conc = data_exp$log10conc),
            FUN = mean)
    colnames(GRavg)[2] = 'GRvalue'
    l = dim(GRavg)[1]
    out['N_conc'] = l

    out['GRmax'] = min(GRavg$GRvalue[c(l,l-1)], na.rm = TRUE)
    diff_vector = diff(GRavg$log10conc, lag = 1)
    out['log10_conc_step'] = mean(diff_vector)
    conc_range = GRavg$log10conc[l] - GRavg$log10conc[1]

    out['GR_AOC'] = mean(1 - GRavg$GRvalue) # use mean for consistency with mean.viability
    # ## Alternative version: Trapezoid rule for integration of GR_AOC
    # out['GR_AOC'] = sum((1 - (GRavg$GRvalue[1:(l-1)]+GRavg$GRvalue[2:l])/2)*
    #                diff_vector, na.rm = TRUE)/conc_range

    # analytical solution for GR50
    out['GR50'] = out['GEC50']*((upper_GR-out['GRinf'])/(0.5-out['GRinf']) - 1)^(1/out['h_GR'])

    # testing the significance of the fit and replacing with flat function if required
    pcutoff = ifelse(force, 1, .05)
    if(!is.na(out['pval_GR'])) {
        out['flat_fit_GR'] = ifelse(out['pval_GR'] >= pcutoff |
                                 is.na(out['GEC50']), 1, 0)
    } else {
        out['flat_fit_GR'] = ifelse(is.na(out['GEC50']), 1, 0)
    }

    # Replace values for flat fits: GEC50 = 0, h_GR = 0.01 and GR50 = +/- Inf
    if(out['flat_fit_GR'] == 1) {
          out['GEC50'] = 0
          out['h_GR'] = 0.0001
          out['GR50'] = ifelse(mean(GRavg$GRvalue) > .5, Inf, -Inf)
          out['GRinf'] = mean(GRavg$GRvalue)
        }

    # Add GR50 = +/-Inf for any curves that don't reach GR = 0.5
    if(is.na(out['GR50'])) {
        out['GR50'] = ifelse(out['GRinf'] > .5, Inf, -Inf)
    }

    out
}

# logistic function (not used in the file but useful for plotting externally)
#' @export
GRlogistic_4parameters <- function(c, GRinf, GR0, GEC50, h_GR) {
  GRinf + (GR0 - GRinf)/(1 + (c/GEC50)^h_GR)
}
