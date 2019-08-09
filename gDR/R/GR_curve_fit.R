
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
ICGRlogisticFit <- function(log10concs, RelativeViability, GRvalues, e_0 = 1, GR_0 = 1,
                    force = FALSE, cap = FALSE) {

    # Implementation of the genedata approach for curve fit: https://screener.genedata.com/documentation/display/DOC15/Business+Rules+for+Dose-Response+Curve+Fitting+Model+Selection+and+Fit+Validity
    #

    # define variables and prepare data
    IC_data_exp = data.frame(log10conc=log10concs, RelativeViability = RelativeViability)
    GR_data_exp = data.frame(log10conc=log10concs, GRvalue=GRvalues)
    concs = 10**log10concs
    ICfit_parameters = c("h_ic", "e_inf", "ec50")
    GRfit_parameters = c('h_GR', 'GRinf', 'GEC50')

    out = array(NA, length(get_header('metrics_results')))
    names(out) = get_header('metrics_results')
    out['maxlog10Concentration'] = max(log10concs)
    out['N_conc'] = length(unique(log10concs))
    out['e_0'] = e_0
    out['GR_0'] = GR_0

    # fit parameters and boundaries
    ICpriors = c(2, 0.4, median(concs))
    GRpriors = c(2, 0.1, median(concs))
    IClower = c(.1, 0, min(concs)/10)
    GRlower = c(.1, -1, min(concs)/10)
    upper = c(5, 1, max(concs)*10)

    controls = drc::drmc()
    controls$relTol = 1e-06
    controls$errorm = FALSE
    controls$noMessage = TRUE
    controls$rmNA = TRUE

    ######################################
    # IC curve fitting
    output_model_new = try(drc::drm(
            RelativeViability ~ log10conc,
            data=IC_data_exp,
            logDose = 10,
            fct=drc::LL.3u(upper = e_0, names = ICfit_parameters),
            start = ICpriors, lowerl = IClower, upperl = upper, control = controls,
            na.action = na.omit))

    # assuming proper fit result
    if(class(output_model_new)!="try-error") {
        for (p in ICfit_parameters) {
            out[p] = stats::coef(output_model_new)[paste0(p, ':(Intercept)')]
        }
        # F-test for the significance of the sigmoidal fit
        Npara = 3 # N of parameters in the growth curve
        Npara_flat = 1 # F-test for the models
        RSS2 = sum(stats::residuals(output_model_new)^2, na.rm = TRUE)
        RSS1 = sum((IC_data_exp$RelativeViability - mean(IC_data_exp$RelativeViability,
                        na.rm = TRUE))^2, na.rm = TRUE)
        df1 = (Npara - Npara_flat)
        df2 = (length(na.omit(IC_data_exp$RelativeViability)) - Npara + 1)
        f_value = ((RSS1-RSS2)/df1)/(RSS2/df2)
        f_pval = stats::pf(f_value, df1, df2, lower.tail = FALSE)
        out['ic_r2'] = 1 - RSS2/RSS1
    }


    # non-fitted metrics
    ICavg = aggregate(IC_data_exp$RelativeViability, by=list(log10conc = IC_data_exp$log10conc),
            FUN = mean)
    colnames(ICavg)[2] = 'RelativeViability'
    l = dim(ICavg)[1]

    out['e_max'] = min(ICavg$RelativeViability[c(l,l-1)], na.rm = TRUE)

    out['mean_viability'] = mean(ICavg$RelativeViability)

    # analytical solution for ic50
    out['ic50'] = out['ec50']*((e_0-out['e_inf'])/(0.5-out['e_inf']) - 1)^(1/out['h_ic'])

    # testing the significance of the fit and replacing with flat function if required
    pcutoff = ifelse(force, 1, .05)
    if(!is.na(f_pval)) {
        out['flat_fit_ic'] = ifelse(f_pval >= pcutoff |
                                 is.na(out['ec50']), 1, 0)
    } else {
        out['flat_fit_ic'] = ifelse(is.na(out['ec50']), 1, 0)
    }

    # Replace values for flat fits: ec50 = 0, h_ic = 0.01 and ic50 = +/- Inf
    if(out['flat_fit_ic'] == 1) {
          out['ec50'] = 0
          out['h_ic'] = 0.0001
          out['ic50'] = ifelse(mean(ICavg$RelativeViability) > .5, Inf, -Inf)
          out['e_inf'] = mean(ICavg$RelativeViability)
        }

    # Add ic50 = +/-Inf for any curves that don't reach RelativeViability = 0.5
    if(is.na(out['ic50'])) {
        out['ic50'] = ifelse(out['e_inf'] > .5, Inf, -Inf)
    }




    ######################################
    # GR curve fitting
    output_model_new = try(drc::drm(
            GRvalue ~ log10conc,
            data=GR_data_exp,
            logDose = 10,
            fct=drc::LL.3u(upper = GR_0, names = GRfit_parameters),
            start = GRpriors, lowerl = GRlower, upperl = upper, control = controls,
            na.action = na.omit))

    # assuming proper fit result
    if(class(output_model_new)!="try-error") {
        for (p in GRfit_parameters) {
            out[p] = stats::coef(output_model_new)[paste0(p, ':(Intercept)')]
        }
        # F-test for the significance of the sigmoidal fit
        Npara = 3 # N of parameters in the growth curve
        Npara_flat = 1 # F-test for the models
        RSS2 = sum(stats::residuals(output_model_new)^2, na.rm = TRUE)
        RSS1 = sum((GR_data_exp$GRvalue - mean(GR_data_exp$GRvalue, na.rm = TRUE))^2,
                   na.rm = TRUE)
        df1 = (Npara - Npara_flat)
        df2 = (length(na.omit(GR_data_exp$GRvalue)) - Npara + 1)
        f_value = ((RSS1-RSS2)/df1)/(RSS2/df2)
        f_pval = stats::pf(f_value, df1, df2, lower.tail = FALSE)
        out['GR_r2'] = 1 - RSS2/RSS1
    }

    # non-fitted metrics
    GRavg = aggregate(GR_data_exp$GRvalue, by=list(log10conc = GR_data_exp$log10conc),
            FUN = mean)
    colnames(GRavg)[2] = 'GRvalue'
    l = dim(GRavg)[1]

    out['GRmax'] = min(GRavg$GRvalue[c(l,l-1)], na.rm = TRUE)

    out['GR_AOC'] = mean(1 - GRavg$GRvalue) # use mean for consistency with mean.viability
    # ## Alternative version: Trapezoid rule for integration of GR_AOC
    # diff_vector = diff(GRavg$log10conc, lag = 1)
    # conc_range = GRavg$log10conc[l] - GRavg$log10conc[1]
    # out['GR_AOC'] = sum((1 - (GRavg$GRvalue[1:(l-1)]+GRavg$GRvalue[2:l])/2)*
    #                diff_vector, na.rm = TRUE)/conc_range

    # analytical solution for GR50
    out['GR50'] = out['GEC50']*((GR_0-out['GRinf'])/(0.5-out['GRinf']) - 1)^(1/out['h_GR'])

    # testing the significance of the fit and replacing with flat function if required
    pcutoff = ifelse(force, 1, .05)
    if(!is.na(f_pval)) {
        out['flat_fit_GR'] = ifelse(f_pval >= pcutoff |
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
logistic_4parameters <- function(c, Vinf, V0, EC50, h) {
  Vinf + (V0 - Vinf)/(1 + (c/EC50)^h)
}
