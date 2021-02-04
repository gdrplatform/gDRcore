#' getMetaData
#'
#' Divide the columns of an input data.frame into treatment metadata, condition metadata, 
#' experiment metadata, and core data for further analysis. This will most commonly be used
#' to identify the different components of a \linkS4class{SummarizedExperiment} object.
#'
#' @param data data.frame with drug-response data
#' @param discard_keys character vector of keys to discard from the row or column data, 
#' and to leave in the matrix data. See details.
#'
#' @return named list containing different elements of a \linkS4class{SummarizedExperiment};
#' see details.
#' 
#' @details
#' Named list containing the following elements:
#' \itemize{
#'  \item{treatment_md}{treatment metadata}
#'  \item{condition_md}{condition metadata}
#'  \item{data_fields}{all data.frame column names corresponding to fields represented within a BumpyMatrix cell}
#'  \item{experiment_md}{metadata that is constant for all entries of the data.frame}
#' }
#'
#' The \code{discard_keys} provides the user the opportunity to specify that they would not 
#' like to use that metadata field as a differentiator of the treatments, and instead, incorporate it
#' into the \code{DataFrame} in the BumpyMatrix.
#'
#' @export
#'
getMetaData <- function(data, discard_keys = NULL) {
  # Assertions.
  stopifnot(any(inherits(data, "data.frame"), inherits(data, "DataFrame")))
  checkmate::assert_character(discard_keys, null.ok = TRUE)

  data <- methods::as(data, "DataFrame")
  all_cols <- colnames(data)

  ## Identify all known fields.
  data_fields <- c(gDRutils::get_header("raw_data"),
    gDRutils::get_header("normalized_results"),
    gDRutils::get_header("averaged_results"),
    gDRutils::get_header("metrics_results"),
    gDRutils::get_identifier("well_position"),
    "Barcode",
    "Template",
    "Concentration"
  )
  data_cols <- data_fields[data_fields %in% all_cols]

  cell_id <- gDRutils::get_identifier("cellline")
  cell_fields <- c(cell_id, gDRutils::get_header("add_clid"))
  cell_cols <- cell_fields[cell_fields %in% all_cols]

  drug_fields <- c(gDRutils::get_identifier("drug"),
    gDRutils::get_identifier("drugname"),
    gDRutils::get_identifier("duration"),
    gDRutils::get_identifier("drug_moa")
  )
  drug_cols <- drug_fields[drug_fields %in% all_cols]

  ## Separate metadata from data fields.
  # Get data columns.
  data_cols <- c(data_cols, discard_keys) # Assign discarded keys to data.

  meta_cols <- setdiff(all_cols, data_cols) 
  md <- unique(data[, meta_cols]) 

  trt_cols <- setdiff(meta_cols, cell_cols)
  singletons <- vapply(trt_cols,
    function(x) {nrow(unique(md[, x, drop = FALSE])) == 1L},
    logical(1))

  # Get experiment columns.
  constant_cols <- setdiff(trt_cols[singletons], drug_cols) # Protect drug fields.
  exp_md <- unique(md[, constant_cols, drop = FALSE])

  remaining <- setdiff(meta_cols, constant_cols)

  ## Ensure drugs are protected and not accidentally pulled into the cell metadata.
  pattern <- sprintf("^%s*|^%s*|^%s$|^%s$", 
    drug_fields[1], drug_fields[2], drug_fields[3], drug_fields[4])
  remaining <- remaining[!grepl(pattern, remaining)]

  ## Identify cellline properties by checking what columns have only a 1:1 mapping for each cell line.
  cl_entries <- cell_id
  for (j in setdiff(remaining, cell_id)) {
    cl_prop <- split(md[[j]], as.factor(md[, cell_id]))
    cl_prop <- lapply(cl_prop, function(grp) {length(unique(grp)) == 1L})
    if (all(unlist(cl_prop))) {
      cl_entries <- c(cl_entries, j)
    }
  }

  if (!all(present <- cell_cols %in% cl_entries)) {
    warning(sprintf("'%s' not metadata for unique cell line identifier column: '%s'", 
      paste(cell_cols[!present], collapse=", "), cell_id))
  }

  ## Condition metadata.
  condition_md <- unique(md[, cl_entries, drop = FALSE])
  condition_md$col_id <- seq_len(nrow(condition_md))
  rownames(condition_md) <- apply(condition_md, 1, function(x) {paste(x, collapse = "_")})
  condition_md <- condition_md[! names(condition_md) %in% c('col_id')]

  ## Treatment metadata.
  trt_cols <- setdiff(meta_cols, c(cl_entries, constant_cols))
  treatment_md <- unique(md[, trt_cols, drop = FALSE])
  treatment_md$row_id <- seq_len(nrow(treatment_md))
  rownames(treatment_md) <- apply(treatment_md, 1, function(x) {paste(x, collapse = "_")})
  treatment_md <- treatment_md[! names(treatment_md) %in% c('row_id')]

  return(list(
    condition_md = condition_md,
    treatment_md = treatment_md,
    data_fields = data_cols,
    experiment_md = exp_md
  ))
}
