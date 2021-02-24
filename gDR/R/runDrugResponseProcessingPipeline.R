# TODO: Add in all other parameters. 
#' @export
#'
runDrugResponseProcessingPipeline <- function(df_) {
  se <- create_SE2(df_)
  se <- normalize_SE2(se)
  se <- average_SE2(se)
  se <- fit_SE2(se)
  se
}
