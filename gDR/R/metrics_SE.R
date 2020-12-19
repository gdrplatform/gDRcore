metrics_SE <- function(avgSE, studyConcThresh = 4) {

    # Assertions:
    checkmate::assert_class(avgSE, "SummarizedExperiment")
    checkmate::assert_number(studyConcThresh)

    stopifnot(is.numeric(studyConcThresh))
    # this is not used as we enforce the same conditions as the input SE; not collapsing allowed
    # if (is.null(DoseRespKeys)) {
    #     if ("Keys" %in% names(metadata(avgSE))) DoseResp = metadata(avgSE)$Keys$DoseResp
    #     else DoseRespKeys = identify_keys(avgSE)$DoseResp
    # } else {
    #     metadata(avgSE)$Keys$DoseResp = DoseRespKeys
    # }

    metricsSE <- avgSE
    SummarizedExperiment::assay(metricsSE, "Metrics") <- SummarizedExperiment::assay(metricsSE, "Averaged")

    # temporary optimization (use 'normSE_n' and 'normSE_c' to avoid using 'assay<-` in a foor loops)
    # TODO: refactor this part of code once we switch to DataFrameMatrix class
    mSE_m <- SummarizedExperiment::assay(metricsSE, "Metrics")
    a_SE = SummarizedExperiment::assay(metricsSE, "Averaged")
    aCtrl_SE = SummarizedExperiment::assay(metricsSE, "Avg_Controls")
    for (i in rownames(metricsSE)) {
        for (j in colnames(metricsSE)) {
            df_ <- a_SE[[i, j]]
            if (!is.null(df_) && all(dim(df_) > 0)) { # studyConcThresh is embeded in RVGRfits
                mSE_m[[i, j]] <- DataFrame(gDRutils::RVGRfits(df_,
                    e_0 = aCtrl_SE[[i, j]]$RefRelativeViability,
                    GR_0 = aCtrl_SE[[i, j]]$RefGRvalue,
                    n_point_cutoff = studyConcThresh))
            } else {
                out <- DataFrame(matrix(NA, 0, length(gDRutils::get_header("response_metrics"))+2))
                colnames(out) <- c(gDRutils::get_header("response_metrics"), "maxlog10Concentration", "N_conc")
                mSE_m[[i, j]] <- out
            }
        }
    }
    SummarizedExperiment::assay(metricsSE, "Metrics") <- mSE_m
    return(metricsSE)
}
