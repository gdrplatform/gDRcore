#' 'gDRcore - a package with core functionality'
#' 
#' gDRcore uses `BiocParallel` package to make computation parallel using all available cores.
#' It is possible to switch into sequential mode by modifying gDRcore settings.
#' If you want to use a sequential approach you can easily run `options(MulticoreParam = SerialParam())`
#' to turn off parallel computing mode. For other available settings please visit the documentation
#' of `BiocParallel` package:
#' https://bioconductor.org/packages/release/bioc/html/BiocParallel.html
#'
#' @docType package
#' @name gDRcore
#' @import S4Vectors
#' @import checkmate
#' @importFrom gDRutils get_synthetic_data
#' @useDynLib gDRcore
NULL
