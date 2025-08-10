#' Get default nested identifiers
#'
#' @param x data.table with raw data or \code{SummarizedExperiment} object 
#' with gDR assays
#' @param data_model single-agent vs combination
#' 
#' @examples 
#' get_default_nested_identifiers(data.table::data.table())
#' 
#' @return vector of nested identifiers
#' 
#' @keywords utils
#' @export
get_default_nested_identifiers <- function(x, data_model = NULL) {
  UseMethod("get_default_nested_identifiers")
}


#' @export
#' @keywords utils
#' @rdname get_default_nested_identifiers
get_default_nested_identifiers.data.table <- function(x, data_model = NULL) {
  
  checkmate::assert_data_table(x)
  checkmate::assert_choice(
    data_model,
    c("single-agent", "combination", "time-course"), 
    null.ok = TRUE
  )
  
  .get_default_nested_identifiers(se = NULL, data_model = NULL)
  
}

#' @export
#' @keywords utils
#' @rdname get_default_nested_identifiers
get_default_nested_identifiers.SummarizedExperiment <- function(
    x,
    data_model = NULL) {
  
  .get_default_nested_identifiers(x, data_model)
}

#' @keywords utils
.get_default_nested_identifiers <- function(se = NULL, data_model = NULL) {
 
  checkmate::assert_choice(
    data_model, c("single-agent", "time-course", "combination"),
    null.ok = TRUE
  )
 
  ml <- list(`single-agent` = .get_default_single_agent_nested_identifiers(se),
             combination = .get_default_combination_nested_identifiers(se),
             `time-course` = .get_default_combination_nested_identifiers(se))
  if (is.null(data_model)) {
    ml
  } else {
    ml[[data_model]]
  }
  
}

#' @keywords utils
.get_default_single_agent_nested_identifiers <- function(se = NULL) {
  if (is.null(se)) {
    gDRutils::get_env_identifiers("concentration")
  } else {
    gDRutils::get_SE_identifiers(se, "concentration")
  }
}

#' @keywords utils
.get_default_combination_nested_identifiers <- function(se = NULL) {
  identifiers <- if (is.null(se)) {
    gDRutils::get_env_identifiers(
      k = c("concentration", "concentration2"),
      simplify = FALSE
    )
  } else {
    gDRutils::get_SE_identifiers(
      se = se, 
      id_type = c("concentration", "concentration2"),
      simplify = FALSE
    )
  }
  
  as.character(unlist(identifiers))
}







#' Detect model of data in data.table
#'
#' @param x data.table of raw drug response data 
#'          containing both treated and untreated values. 
#'
#' @return string with the information of the raw data follows single-agent or 
#' combination data model
#' @keywords utils
#' @export
data_model.data.table <- function(x) {
  
  checkmate::assert_data_table(x)
  drug_ids <- unlist(
    gDRutils::get_env_identifiers(
      c("drug_name", "drug_name2"), 
      simplify = FALSE
    )
  )
  cl_id <- gDRutils::get_env_identifiers("cellline")
  conc2 <- gDRutils::get_env_identifiers("concentration2")
  if (all(.get_default_combination_nested_identifiers() %in% colnames(x))) {
    if (all(x[[conc2]]
            %in% gDRutils::get_env_identifiers("untreated_tag"))) {
      "single-agent"
    } else {
      "combination"
    }
  } else if (.get_default_single_agent_nested_identifiers() %in% colnames(x)) {
    "single-agent"
  } else {
    "time-course"
    # stop("Unknown data model")
  }
}

#' Detect model of data from experiment name
#'
#' @param x character with experiment name
#'
#' @return string with the information of the raw data follows single-agent or 
#' combination data model
#' @keywords utils
#' @export
data_model.character <- function(x) {
  checkmate::assert_subset(x,  c('single-agent','combination','co-dilution', 'time-course'))
 
  exp_v <- gDRutils::get_experiment_groups()
  exp_v$`time-course` = "time-course"
  names(exp_v[grep(x, exp_v)])
}




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

  data_fields <- c(md$data_fields, "Duration", "row_id", "col_id", "swap_sa")


  out <- vector("list", length = nrow(treated))
  out <- lapply(seq_len(nrow(treated)), function(i) {
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
    untrt_cols <- intersect(c("CorrectedReadout", "record_id", 'Duration', nested_confounders), names(dfs))
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
  # trt_keep[1] <- TRUE
  # trt_keep[2] <- TRUE
  # print(trt_out)
  trt_mat <- BumpyMatrix::splitAsBumpyMatrix(
    trt_out[, trt_keep, drop = FALSE],
    row = trt_out$row_id, col = trt_out$col_id
  )
  ctl_mat <- BumpyMatrix::splitAsBumpyMatrix(
    ctl_out[, ctl_keep, drop = FALSE],
    row = ctl_out$row_id, col = ctl_out$col_id
  )
  matsL <- list(RawTreated = trt_mat, Controls = ctl_mat)
  print(rownames(trt_mat))
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
    
  ta <- assay(se, 'RawTreated')
  
  stopifnot(all(ta[1,1,1] == 0))
  log_cell_counts_0 <- log(unlist(ta[1,1,6]))
  
  normalized_cell_counts = c()
  for (i in 1:length(ta)) {
    norm_log_cell_counts <- log(unlist(ta[i,1,6])) - log_cell_counts_0
    normalized_cell_counts <- c(normalized_cell_counts, norm_log_cell_counts)
  }
  ta <- unsplitAsDataFrame(ta)
  ta$NormReadoutValue <- normalized_cell_counts
  
  return(ta)
}


