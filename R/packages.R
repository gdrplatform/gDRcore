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
      "normalization_type",
      "cotrt_value",
      "x",
      "Duration",
      "isDay0",
      "record_id",
      "ratio"
    ), 
    utils::packageName())
}

