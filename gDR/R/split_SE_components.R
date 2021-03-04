#' split_SE_components
#'
#' Divide the columns of an input data.frame into treatment metadata, condition metadata, 
#' experiment metadata, and core data for further analysis. This will most commonly be used
#' to identify the different components of a \linkS4class{SummarizedExperiment} object.
#'
#' @param df_ data.frame with drug-response data
#' @param nested_keys character vector of keys to discard from the row or column data, 
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
#' The \code{nested_keys} provides the user the opportunity to specify that they would not 
#' like to use that metadata field as a differentiator of the treatments, and instead, incorporate it
#' into the \code{DataFrame} in the BumpyMatrix.
#'
#' In the event that if any of the \code{nested_keys} are constant throughout the whole data.frame, 
#' they will still be included in the DataFrame of the BumpyMatrix and not in the experiment_metadata.
#'
#' @export
#'
split_SE_components <- function(df_, nested_keys = NULL) {
  # Assertions.
  stopifnot(any(inherits(df_, "data.frame"), inherits(df_, "DataFrame")))
  checkmate::assert_character(nested_keys, null.ok = TRUE)

  df_ <- S4Vectors::DataFrame(df_)
  all_cols <- colnames(df_)

  ## Identify all known fields.
  data_fields <- c(gDRutils::get_header("raw_data"),
    gDRutils::get_header("normalized_results"),
    gDRutils::get_header("averaged_results"),
    gDRutils::get_header("metrics_results"),
    gDRutils::get_identifier("well_position"),
    "Template",
    "Concentration",
    nested_keys
  )
  data_fields <- unique(data_fields)
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
  meta_cols <- setdiff(all_cols, data_cols) 
  md <- unique(df_[, meta_cols]) 

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
  condition_md <- condition_md[sort(colnames(condition_md))] # TODO: can delete me later. Sorting for simple comparison. 
  rownames(condition_md) <- apply(condition_md, 1, function(x) {paste(x, collapse = "_")})
  condition_md <- condition_md[! names(condition_md) %in% c('col_id')]

  ## Treatment metadata.
  trt_cols <- setdiff(meta_cols, c(cl_entries, constant_cols))
  # Adjust the order of the columns.
  drug_trt_cols <- intersect(drug_cols, trt_cols)
  trt_cols <- c(drug_trt_cols, setdiff(trt_cols, drug_trt_cols))
  treatment_md <- unique(md[, trt_cols, drop = FALSE])
  treatment_md$row_id <- seq_len(nrow(treatment_md))
  treatment_md <- treatment_md[sort(colnames(treatment_md))] # TODO: can delete me later. Sorting for simple comparison. 
  rownames(treatment_md) <- apply(treatment_md, 1, function(x) {paste(x, collapse = "_")})
  treatment_md <- treatment_md[! names(treatment_md) %in% c('row_id')]

  return(list(
    condition_md = condition_md,
    treatment_md = treatment_md,
    data_fields = data_cols,
    experiment_md = exp_md
  ))
}
