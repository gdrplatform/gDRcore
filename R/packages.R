#' @import S4Vectors
#' @import checkmate
#' @importFrom data.table ":=" ".SD"
#' @useDynLib gDRcore

# data.table awareness
.datatable.aware <- TRUE

#' onload function
#'
#' @param libname library name
#' @param pkgname package name
#' @noRd
.onLoad <- function(libname, pkgname) {
  # scientific notation was disabled due to the problem with unit tests
  options(scipen = 999) 
  
  cores <- Sys.getenv("NUM_CORES")
  if (cores != "") {
    BiocParallel::register(
      BiocParallel::MulticoreParam(workers = as.numeric(cores)), 
      default = TRUE
    )
  }
}

# Prevent R CMD check from complaining about the use of pipe expressions
# standard data.table variables
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      "..codrug_cols",
      "..col_intersect",
      "..colnames",
      "..cols",
      "..cols_pairs",
      "..conc_ids",
      "..conc_idx",
      "..cotrt_id",
      "..drug_id",
      "..drugs_cotrt_ids",
      "..drugs_ids",
      "..duration_col",
      "..k",
      "..keep",
      "..measured_col",
      "..metric_col",
      "..out_col_name",
      "..overridden",
      "..present_ref_cols",
      "..row_order_col",
      "..selected_columns",
      "..series_id",
      "..series_id1",
      "..series_id2",
      "..series_identifiers",
      "..standardized_conc_col",
      "..untrt_cols",
      "..valid",
      "..x",
      "..y",
      "normalization_type",
      "cotrt_value",
      "x"
    ), 
    utils::packageName())
}

