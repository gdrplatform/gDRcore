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

#' @export
#'
runDrugResponseProcessingPipeline1 <- function(df_) {
  se <- normalize_SE(df_)
  se <- average_SE(se)
  se <- fit_SE(se)
  se
}
