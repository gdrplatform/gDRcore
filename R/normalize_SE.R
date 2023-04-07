#' @rdname runDrugResponseProcessingPipelineFxns
#' @examples 
#' td <- gDRimport::get_test_data()
#' l_tbl <- gDRimport::load_data(td$m_file, td$t_files, td$r_files)
#' imported_data <- merge_data(
#'   l_tbl$manifest, 
#'   l_tbl$treatments, 
#'   l_tbl$data
#' )
#'
#' se <- create_SE(imported_data, data_type = "single-agent")
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
                         control_assay = "Controls", 
                         raw_treated_assay = "RawTreated", 
                         normalized_assay = "Normalized",
                         ndigit_rounding = 4) {
  
  # Assertions
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_character(data_type)
  checkmate::assert_character(nested_identifiers, null.ok = TRUE)
  checkmate::assert_character(nested_confounders, null.ok = TRUE)
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

  
  cl_names <- cdata[, 
    gDRutils::get_SE_identifiers(
      se, 
      "cellline_name", 
      simplify = TRUE
    ), 
    drop = FALSE
  ]
  
  cell_ref_div_col <- gDRutils::get_SE_identifiers(se, 
                                                   "cellline_ref_div_time", 
                                                   simplify = TRUE
                                                   )
  if (!cell_ref_div_col %in% names(cdata)) {
    cdata[[cell_ref_div_col]] <- NA
  }
  
  cl_ref_div_times <- cdata[,
                            cell_ref_div_col,
                            drop = FALSE
                            ]
  durations <- rdata[, 
    gDRutils::get_SE_identifiers(
      se, 
      "duration", 
      simplify = TRUE
    ), 
    drop = FALSE
  ]

  refs <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assays(se)[[control_assay]]
  )
  trt <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assays(se)[[raw_treated_assay]]
  )
  
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
    cl_name <- cl_names[j, ]
    ref_div_time <- cl_ref_div_times[j, ]

    duration <- durations[i, ]

    ref_df <- refs[refs$row == i & refs$column == j, ]
    trt_df <- trt[trt$row == i & trt$column == j, ]

    # pad the ref_df for missing values based on nested_keys 
    # (uses mean across all available values)
    if (!is.null(nested_keys) && length(nested_keys) > 0) {
      ref_df <- fill_NA_ref(ref_df, nested_keys)
    }
    
    # Merge to ensure that the proper discard_key values are mapped.
    all_readouts_df <- merge(trt_df, 
      ref_df, 
      by = nested_confounders,
      all.x = TRUE)

    normalized <- S4Vectors::DataFrame(
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
    normalized <- cbind(all_readouts_df[keep], normalized) 
    normalized$row_id <- i
    normalized$col_id <- j
    normalized$id <- as.character(seq_len(nrow(normalized)))
    normalized <- reshape2::melt(as.data.frame(normalized),
                                 measure.vars = norm_cols,
                                 variable.name = "normalization_type",
                                 value.name = "x")
    rownames(normalized) <- paste(
      normalized$id, 
      normalized$normalization_type, 
      sep = "_"
    )
    normalized$id <- NULL
    S4Vectors::DataFrame(normalized)
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
fill_NA_ref <- function(ref_df, nested_keys) {
  data_columns <- setdiff(colnames(ref_df), c(nested_keys, "row", "column"))
  ref_cols <- data.frame(ref_df[, data_columns, drop = FALSE])
  
  if (any(!is.na(ref_cols))) {
    ref_df_mean <- colMeans(
      data.frame(ref_df[, data_columns, drop = FALSE]), na.rm = TRUE
    )
    for (col in data_columns) {
      ref_df[is.na(ref_df[, col]), col] <- ref_df_mean[[col]]
    }
  }
  S4Vectors::DataFrame(ref_df)
}