#' average_SE
#'
#' Avereage normalized SummarizedExperiment of DR data
#'
#' @param normSE a SummarizedExperiment with normalized DR data
#' @param TrtKeys a vector of keys used for averaging (NULL by default)
#'
#' @return a SummarizedExperiment with additional assay with averaged DR data
#' @export
#'

average_SE <- function(normSE, TrtKeys = NULL, include_masked = F) {

  # Assertions:
  checkmate::assert_class(normSE, "SummarizedExperiment")
  checkmate::assert_vector(TrtKeys, null.ok = TRUE)


    avgSE <- normSE
    if (is.null(TrtKeys)) {
        if ("Keys" %in% names(metadata(normSE))) {
          TrtKeys <- metadata(normSE)$Keys$Trt
          TrtKeys <- setdiff(TrtKeys, metadata(normSE)$Keys$discard_keys)
        } else {
          TrtKeys <- identify_keys(normSE)$Trt
        }
    }
    metadata(normSE)$Keys$Trt <- TrtKeys

    SummarizedExperiment::assay(avgSE, "Averaged") <- SummarizedExperiment::assay(avgSE, "Normalized")
    avgSE <- aapply(avgSE, function(x) {
        # bypass 'masked' filter
        x$masked = x$masked & !include_masked

        subKeys <- intersect(TrtKeys, colnames(x))
        if (sum(!x$masked) >= 1) {
            df_av <- aggregate(x[ !x$masked ,
                                  c("GRvalue", "RelativeViability","CorrectedReadout")],
                            by = as.list(x[ !x$masked , subKeys, drop = FALSE]),
                            FUN = function(y) mean(y, na.rm = TRUE))
            df_std <- aggregate(x[!x$masked, c("GRvalue", "RelativeViability")],
                                by = as.list(x[ !x$masked, subKeys, drop = FALSE]),
                                FUN = function(x) sd(x, na.rm = TRUE))
            colnames(df_std)[colnames(df_std) %in% c("GRvalue", "RelativeViability")] =
                paste0("std_",
                    colnames(df_std)[colnames(df_std) %in% c("GRvalue", "RelativeViability")])
            return( merge(df_av, df_std, by = subKeys) )
        } else { # case: (nrow(x) == 0 || all(x$masked))
            df_ = as.data.frame(matrix(0,0,length(subKeys)+5))
            colnames(df_) = c(subKeys,
                  c("GRvalue", "RelativeViability","CorrectedReadout"),
                  paste0("std_", c("GRvalue", "RelativeViability")))
            return(df_)
        } 
    }, "Averaged")

    SummarizedExperiment::assay(avgSE, "Avg_Controls") <- SummarizedExperiment::assay(avgSE, "Controls")
    avgSE <- aapply(avgSE, function(x) {
        if (nrow(x) > 1) {
            subKeys <- intersect(TrtKeys, colnames(x))
            df_av <- DataFrame(lapply(x[, c("Day0Readout", "UntrtReadout",
                    "RefGRvalue", "RefRelativeViability",
                    "RefReadout", "DivisionTime")], FUN = function(y) mean(y, na.rm = TRUE)))
            return( df_av )
        } else return(x)
    }, "Avg_Controls")

    return(avgSE)
}
