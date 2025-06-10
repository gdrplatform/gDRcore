#' @examples
#' 
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
#' se <- purrr::quietly(create_SE)(imported_data, data_type = "single-agent")
#' 
#' 
#' @export
#' @rdname runDrugResponseProcessingPipelineFxns
#'
create_SE <- function(df_,
                      data_type,
                      readout = "ReadoutValue",
                      nested_identifiers = NULL,
                      nested_confounders = intersect(
                        names(df_),
                        gDRutils::get_env_identifiers("barcode")
                      ),
                      override_untrt_controls = NULL) {
  # Assertions:
  checkmate::assert_multi_class(df_, c("data.table", "DataFrame"))
  checkmate::assert_string(data_type)
  checkmate::assert_string(readout)
  checkmate::assert_character(nested_identifiers, null.ok = TRUE)
  checkmate::assert_character(nested_confounders, null.ok = TRUE)
  checkmate::assert_vector(override_untrt_controls, null.ok = TRUE)
  
  
  if (length(nested_confounders) == 0) {
    nested_confounders <- NULL
  }
  
  if (is.null(nested_identifiers)) {
    nested_identifiers <-
      get_default_nested_identifiers(df_)[[data_model(data_type)]]
  }
  
  untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")
  
  nested_keys <- intersect(c(nested_identifiers, nested_confounders, "record_id"),
                           names(df_))
  identifiers <- gDRutils::get_env_identifiers()
  Keys <- identify_keys(df_, nested_keys, override_untrt_controls, identifiers)

  if (identifiers$masked_tag %chin% colnames(df_)) {
    df_ <- df_[!get(identifiers$masked_tag)]
    df_[, (identifiers$masked_tag) := NULL]
  }

  # Remove background value from readout (at least 1e-10 to avoid artefactual 
  # normalized values).
  if (any("BackgroundValue" == colnames(df_))) {
    df_[, CorrectedReadout := pmax(.SD[[readout]] - df_$BackgroundValue, 1e-10)]
  } else {
    df_[, CorrectedReadout := .SD[[readout]]]
  }
  
  # overwrite "drug", "drug_name", "drug_moa"
  # with "untreated" if "concentration2" == 0
  if (any(gDRutils::get_env_identifiers("concentration2") == colnames(df_))) {
    single_agent_idx <-
      df_[[gDRutils::get_env_identifiers("concentration2")]] %in% 0
    drug_cols <- c("drug", "drug_name", "drug_moa")
    drug2_var <- intersect(
      unlist(gDRutils::get_env_identifiers(
        paste0(drug_cols, "2"), simplify = FALSE)
      ),
      colnames(df_)
    )
    df_[single_agent_idx, (drug2_var) := untreated_tag[1]]
  }
  
  
  df_[, (names(df_)) := lapply(.SD, function(x) {
    if (is.character(x)) {
      gsub(paste(untreated_tag, collapse = "|"), untreated_tag[1], x)
    } else {
      x
    }
  }), .SDcols = names(df_)]

  # Identify treatments, conditions, and experiment metadata.
  md <- gDRutils::split_SE_components(df_, nested_keys = Keys$nested_keys)
  
  coldata <- md$condition_md
  rowdata <- md$treatment_md
  exp_md <- md$experiment_md

  mapping_entries <- .create_mapping_factors(rowdata, coldata)

  refs <- .map_references(mapping_entries, rowData_colnames = colnames(rowdata))
  emptyRefs <- all(is.null(unlist(refs)))
  trt_conditions <- names(refs)
  sa_conditions <- unique(unname(unlist(refs)))
  
  treated <- mapping_entries[as.numeric(trt_conditions), ]
  # not all entries are mapped on trt_conditions or sa_conditions 
  #   --> untreated should be mapped explicitely to avoid issues with treatments
  #       being considered as untreated in specific combination cases
  untreated <- mapping_entries[map_untreated(mapping_entries), ]

  ## Map controls.
  
  controls <- list(untrt_Endpoint = "untrt_Endpoint", Day0 = "Day0")
  
  ctl_maps <- gDRutils::loop(controls, function(ctl_type) {
    map_df(
      treated, 
      untreated, 
      override_untrt_controls = override_untrt_controls,
      ref_cols = Keys[[ctl_type]], 
      ref_type = ctl_type
    )
  })
  
  ## Combine all controls with respective treatments.
  # Merge raw data back with groupings.
  dfs <- mapping_entries[df_, on = c(colnames(rowdata), colnames(coldata)), allow.cartesian = TRUE]
  data.table::setkey(dfs, rn)

  ## The mapping_entries contain all exhaustive combinations of treatments 
  ## and cells. Not all conditions will actually exist in the data, so filter 
  ## out those that do not exist. 
  treated_rows <- which(treated$rn %in% dfs$rn)
  treated <- treated[treated_rows, ]
  untreated <- dfs[dfs$rn %in% unique(unlist(ctl_maps)), ]
  
  data_fields <- c(md$data_fields, "row_id", "col_id", "swap_sa")

  
  out <- vector("list", length = nrow(treated))
  out <- gDRutils::loop(seq_len(nrow(treated)), function(i) {
    trt <- treated$rn[i]
    
    trt_df <- dfs[trt]
    row_id <- unique(trt_df[["row_id"]])
    refs_df <- dfs[refs[[trt]]]
    
    trt_df <- 
      validate_mapping(trt_df, refs_df, nested_confounders)
    selected_columns <- names(trt_df) %in% data_fields
    trt_df <- trt_df[, selected_columns, with = FALSE]

    # Override the row_id of the references.
    trt_df$row_id <- row_id
    
    untrt_ref <- ctl_maps[["untrt_Endpoint"]][[trt]]
    untrt_cols <- intersect(c("CorrectedReadout", "record_id", nested_confounders), names(dfs))
    untrt_df <- untreated[untreated$rn == untrt_ref[1], untrt_cols, with = FALSE]
    untrt_df <- if (nrow(untrt_df) == 0) {
      data.table::data.table(CorrectedReadout = NA, isDay0 = FALSE)
    } else {
      untrt_df
    }
    untrt_df$isDay0 <- FALSE
    
    day0_ref <- ctl_maps[["Day0"]][[trt]]
    day0_df <- untreated[untreated$rn %chin% day0_ref]
    isDay0 <- day0_df[[gDRutils::get_env_identifiers("duration")]] == 0
    
    day0_df <- day0_df[isDay0, untrt_cols, with = FALSE]
    day0_df <- if (nrow(day0_df) == 0) {
      data.table::data.table(CorrectedReadout = NA, isDay0 = FALSE)
    } else {
      df <- day0_df
      df$isDay0 <- TRUE
      df
    } 
    
    ctl_df <- data.table::rbindlist(list(UntrtReadout = untrt_df,
                                         Day0Readout = day0_df),
                                    fill = TRUE,
                                    idcol = "control_type")
    ctl_df[is.na(isDay0), isDay0 := FALSE]
    
    row_id <- unique(trt_df$row_id)
    col_id <- unique(trt_df$col_id)
    if (length(row_id) != 1L || length(col_id) != 1L) {
      stop(sprintf("non-unique row_ids: '%s' and col_ids: '%s'", 
                   paste0(row_id, collapse = ", "), paste0(col_id, collapse = ", ")))
    }
    ctl_df$row_id <- row_id
    ctl_df$col_id <- col_id

    list(ctl_df = ctl_df,
         trt_df = trt_df)
  })
  
  trt_out <- rbindParallelList(out, "trt_df")
  ctl_out <- rbindParallelList(out, "ctl_df")
  
  trt_keep <- !colnames(trt_out) %in% c("row_id", "col_id")
  ctl_keep <- !colnames(ctl_out) %in% c("row_id", "col_id")

  trt_mat <- BumpyMatrix::splitAsBumpyMatrix(
    trt_out[, trt_keep, drop = FALSE],
    row = trt_out$row_id, col = trt_out$col_id
  )
  ctl_mat <- BumpyMatrix::splitAsBumpyMatrix(
    ctl_out[, ctl_keep, drop = FALSE],
    row = ctl_out$row_id, col = ctl_out$col_id
  )
  matsL <- list(RawTreated = trt_mat, Controls = ctl_mat)

  # Filter out to 'treated' conditions only.
  trt_rowdata <- rowdata[rownames(trt_mat), , drop = FALSE] 
  stopifnot(nrow(trt_rowdata) > 0)
  stopifnot(nrow(trt_rowdata) == length(unique(trt_out$row_id)))
 
  se <- SummarizedExperiment::SummarizedExperiment(assays = matsL,
    colData = coldata[match(colnames(trt_mat), rownames(coldata)), ],
    rowData = trt_rowdata,
    metadata = list())
  
  # Capture important values in experiment metadata.
  se <- gDRutils::set_SE_identifiers(se, identifiers)
  se <- gDRutils::set_SE_experiment_metadata(se, exp_md)
  se <- gDRutils::set_SE_keys(se, lapply(Keys, sort))

  se
}

#' @keywords internal
validate_mapping <- function(trt_df, refs_df, nested_confounders) {
  
  if (!is.null(nested_confounders)) {
    matching_confounders <- refs_df[[nested_confounders]] %in% unique(trt_df[[nested_confounders]])
    if (any(matching_confounders)) {
      refs_df <- refs_df[matching_confounders, ]
    }
  }
  
  drug_id <- gDRutils::get_env_identifiers("drug")
  drug2_id <- gDRutils::get_env_identifiers("drug2")
  conc2 <- gDRutils::get_env_identifiers("concentration2")
  untrt_tag <- gDRutils::get_env_identifiers("untreated_tag")
  trt_df <- rbind(
    trt_df, 
    refs_df[
      refs_df[[drug_id]] %in% 
        c(unique(c(trt_df[[drug_id]], trt_df[[drug2_id]])), untrt_tag), 
    ]
  )
  
  # Swap concentrations for single-agent with drug2
  if (any(conc2 == colnames(trt_df))) {
    swap_idx <- trt_df[[drug_id]] %in%
      setdiff(unique(trt_df[[drug2_id]]), untrt_tag)
    trt_df[swap_idx, "swap_sa"] <- TRUE
  }
  trt_df
}
