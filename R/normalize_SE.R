#' @rdname runDrugResponseProcessingPipelineFxns
#' @examples 
#' td <- gDRimport::get_test_data()
#' l_tbl <- gDRimport::load_data(
#'   manifest_file = gDRimport::manifest_path(td), 
#'   df_template_files = gDRimport::template_path(td), 
#'   results_file = gDRimport::result_path(td)
#' )
#' imported_data <- merge_data(
#'   l_tbl$manifest, 
#'   l_tbl$treatments, 
#'   l_tbl$data
#' )
#'
#' inl <- prepare_input(imported_data)
#' se <- create_SE(
#'  inl$df_list[["single-agent"]],
#'  data_type = "single-agent",
#'  nested_confounders = inl$nested_confounders)
#'  
#' normalize_SE(se, data_type = "single-agent")
#' @export
#'
normalize_SE <- function(se, 
                         data_type,
                         nested_identifiers = NULL,
                         nested_confounders = 
                           gDRutils::get_SE_identifiers(
                             se, 
                             "barcode", 
                             simplify = TRUE
                           ),
                         control_mean_fxn = function(x) {
                           mean(x, trim = 0.25)
                         },
                         control_assay = "Controls", 
                         raw_treated_assay = "RawTreated", 
                         normalized_assay = "Normalized",
                         ndigit_rounding = 4) {
  
  # Assertions
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_character(data_type)
  checkmate::assert_character(nested_identifiers, null.ok = TRUE)
  checkmate::assert_character(nested_confounders, null.ok = TRUE)
  checkmate::assert_function(control_mean_fxn)
  checkmate::assert_string(control_assay)
  checkmate::assert_string(raw_treated_assay)
  checkmate::assert_string(normalized_assay)
  checkmate::assert_number(ndigit_rounding)
  
  gDRutils::validate_se_assay_name(se, control_assay)
  gDRutils::validate_se_assay_name(se, raw_treated_assay)
  
  
  if (is.null(nested_identifiers)) {
    nested_identifiers <- get_default_nested_identifiers(
      se, 
      data_model(data_type)
    )
  }
  
  # Keys
  nested_keys <- c(nested_identifiers, nested_confounders)
  masked_key <- gDRutils::get_SE_identifiers(se, "masked_tag", simplify = TRUE)
  trt_keys <- gDRutils::get_SE_keys(se, key_type = "Trt")
  
  cdata <- SummarizedExperiment::colData(se)
  rdata <- SummarizedExperiment::rowData(se)
  
  rn_cdata <- rownames(cdata)
  rn_rdata <- rownames(rdata)
  
  cdata <- data.table::as.data.table(cdata)
  rdata <- data.table::as.data.table(rdata)
  
  cdata[, rn := rn_cdata]
  rdata[, rn := rn_rdata]
  
  data.table::setkey(cdata, rn)
  data.table::setkey(rdata, rn)
  
  cl_names <- gDRutils::get_SE_identifiers(
    se, 
    "cellline_name", 
    simplify = TRUE
  )
  
  cell_ref_div_col <- gDRutils::get_SE_identifiers(se, 
                                                   "cellline_ref_div_time", 
                                                   simplify = TRUE
                                                   )
  if (!any(cell_ref_div_col == names(cdata))) {
    cdata[[cell_ref_div_col]] <- NA
  }

  duration_name <- gDRutils::get_SE_identifiers(
    se, 
    "duration", 
    simplify = TRUE
  )
  
  refs <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assays(se)[[control_assay]]
  ))
  data.table::setkey(refs, row, column)
  
  trt <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assays(se)[[raw_treated_assay]]
  ))
  data.table::setkey(trt, row, column)
  
  if (any("swap_sa" == names(trt))) {
    conc <- gDRutils::get_env_identifiers("concentration")
    conc2 <- gDRutils::get_env_identifiers("concentration2")
    swap_idx <- !is.na(trt$swap_sa)
    if (any(swap_idx)) {
    trt[swap_idx, c(conc, conc2) := .(get(conc2), get(conc))]
    }
  }
  
  refs$record_id <- trt$record_id <- trt$swap_sa <- NULL
  
  # Extract common nested_confounders shared by trt_df and ref_df
  nested_confounders <- Reduce(intersect, list(nested_confounders,
                                               names(trt),
                                               names(refs)))
  
  norm_cols <- c("RV", "GR")
  out <- vector("list", nrow(se) * ncol(se))
  ref_rel_viability <- ref_GR_value <- div_time <- 
    matrix(NA, nrow = nrow(se), ncol = ncol(se), dimnames = dimnames(se))
  msgs <- NULL
  iterator <- unique(rbind(refs[, c("column", "row")],
                           trt[, c("column", "row")]))
  
  
  # Column major order, so go down first.

  out <- gDRutils::loop(seq_len(nrow(iterator)), function(row) {
    
    x <- iterator[row, ]
    i <- x[["row"]]
    j <- x[["column"]]
    cdata_subset <- cdata[j]
    rdata_subset <- rdata[i]
    
    cl_name <- cdata_subset[[cl_names]]
    ref_div_time <- cdata_subset[[cell_ref_div_col]]

    duration <- rdata_subset[[duration_name]]

    ref_df <- refs[.(i, j)]
    trt_df <- trt[.(i, j)]

    all_readouts_df <- merge_trt_with_ref(ref_df,
                                          trt_df,
                                          nested_confounders,
                                          control_mean_fxn)
    normalized <- data.table::data.table(
      matrix(NA, nrow = nrow(trt_df), ncol = length(norm_cols))
    )
    colnames(normalized) <- c(norm_cols)

    # Normalized treated.
    normalized$RV <- round(
      all_readouts_df$CorrectedReadout / all_readouts_df$UntrtReadout, 
      ndigit_rounding
    )
    
    
    normalized$GR <- calculate_GR_value(
      rel_viability = normalized$RV,
      corrected_readout = all_readouts_df$CorrectedReadout,
      day0_readout = all_readouts_df$Day0Readout,
      untrt_readout = all_readouts_df$UntrtReadout,
      ndigit_rounding = ndigit_rounding,
      duration = duration,
      ref_div_time = ref_div_time
    )
    
    if (any(is.na(all_readouts_df$Day0Readout))) {
      msgs <- c(msgs,
                sprintf("could not calculate GR values for '%s'", cl_name))
    }
    
    # Carry over present treated keys.
    keep <- intersect(
      c(nested_identifiers, trt_keys, masked_key), colnames(all_readouts_df)
    )
    normalized <- cbind(all_readouts_df[, keep, with = FALSE], normalized) 
    normalized$row_id <- i
    normalized$col_id <- j
    normalized$id <- as.character(seq_len(nrow(normalized)))
    normalized <- data.table::melt(normalized,
                                   measure.vars = norm_cols,
                                   variable.name = "normalization_type",
                                   value.name = "x")
    rownames <- paste(
      normalized$id, 
      normalized$normalization_type, 
      sep = "_"
    )
    normalized$id <- NULL
    S4Vectors::DataFrame(normalized, row.names = rownames)
  })
  
  if (!is.null(msgs)) {
    futile.logger::flog.warn(paste0(msgs, collapse = "\n"))
  }
  
  out <- S4Vectors::DataFrame(do.call("rbind", out))
  
  norm <- BumpyMatrix::splitAsBumpyMatrix(
    out[!colnames(out) %in% c("row_id", "col_id")], 
    row = out$row_id, 
    col = out$col_id
  )
  
  SummarizedExperiment::assays(se)[[normalized_assay]] <- norm
  se
}
  
#' @keywords internal
aggregate_ref <- function(ref_df, control_mean_fxn) {
  
  checkmate::assert_data_table(ref_df)
  data_columns <- setdiff(colnames(ref_df), c("row", "column", "masked", "isDay0"))
  corr_readout <- "CorrectedReadout"
  ref_cols <- ref_df[, data_columns, with = FALSE]
  group_cols <- setdiff(names(ref_cols), corr_readout)
  additional_cov <- setdiff(group_cols, "control_type")
  aggregate_formula <- stats::reformulate("control_type",
                                          ifelse(length(additional_cov) == 0, ".", additional_cov))
  
  ref_df_aggregate <- unique(ref_cols[, (control_mean_fxn(get(corr_readout))), by = eval(group_cols)])
  ref_df_dcast <- data.table::dcast(ref_df_aggregate,
                                    aggregate_formula,
                                    value.var = "V1")
  cleaned_df <- ref_df_dcast[!rowSums(is.na(ref_df_dcast)) >=
                               length(setdiff(names(ref_df_dcast), group_cols)), ]
  fill_NA_by_mean(cleaned_df, cleaned_df, unique(ref_df$control_type))
}

#' @keywords internal
fill_NA_by_mean <- function(dt, ref_df, cols) {
  dt2 <- data.table::copy(dt)
  for (col in cols) {
    dt2[is.na(dt2[[col]]), (col) := mean(ref_df[[col]], na.rm = TRUE)]
  }
  dt2
}

#' @keywords internal
merge_trt_with_ref <- function(ref_df,
                               trt_df,
                               nested_confounders,
                               control_mean_fxn) {
  
  # pad the ref_df for missing values based on nested_keys 
  # (uses mean across all available values)
  control_types <- unique(ref_df$control_type)
  ref_df <- aggregate_ref(ref_df, control_mean_fxn = control_mean_fxn)

  # Merge to ensure that the proper discard_key values are mapped.
  all_readouts_df <- if (length(nested_confounders)) {
    ref_df[trt_df, on = nested_confounders]
  } else {
    cbind(trt_df, ref_df)
  }
  
  # Backfill missing values when the `nested_keys` are not matching with an average. 
  # This is necessary if a control is only present on another plate
  fill_NA_by_mean(all_readouts_df, ref_df, control_types)
}
