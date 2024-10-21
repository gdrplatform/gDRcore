#' @noRd
#' @keywords internal
#'
.create_mapping_factors <- function(rowdata, coldata) {
  
  col_id <- rownames(coldata)
  row_id <- rownames(rowdata)
  
  mapping_grid <- data.table::CJ(
    col_id = col_id,
    row_id = row_id
  )
  
  rowdata_dt <- data.table::as.data.table(rowdata)
  coldata_dt <- data.table::as.data.table(coldata)
  
  rowdata_dt[, row_id := row_id]
  coldata_dt[, col_id := col_id]
  
  dt <- coldata_dt[rowdata_dt[mapping_grid, on = "row_id"], on = "col_id"]
  dt[, rn := as.character(.I)]
}


#' cleanup_metadata
#'
#' Cleanup a data.table with metadata
#'
#' @param df_metadata a data.table with metadata
#'
#' @examples
#'
#' df <- data.table::data.table(
#'   clid = "CELL_LINE",
#'   Gnumber = "DRUG_1",
#'   Concentration = c(0, 1),
#'   Duration = 72
#' )
#' cleanup_df <- cleanup_metadata(df)
#'
#' @return a data.table with cleaned metadata
#' @details Adds annotations and check whether user provided correct input data.
#' @keywords utils
#' @export
#'
cleanup_metadata <- function(df_metadata) {
  
  # Assertions:
  checkmate::assert_data_table(df_metadata)

  # Round duration to 6 values. 
  duration_id <- gDRutils::get_env_identifiers("duration")
  df_metadata[[duration_id]] <- round(as.numeric(df_metadata[[duration_id]], 6))
  
  if (!gDRutils::get_env_identifiers("cellline_name") %in% names(df_metadata)) {
    df_metadata <- add_CellLine_annotation(df_metadata)
  }
  
  drug_conc_cols <- unlist(
    gDRutils::get_env_identifiers(
      c(
        "drug", 
        "drug2", 
        "drug3", 
        "concentration", 
        "concentration2", 
        "concentration3"
      ),
      simplify = FALSE
    )
  )
  drug_conc_cols <- drug_conc_cols[drug_conc_cols %in% names(df_metadata)]
  split_idx <- stringr::str_extract(names(drug_conc_cols), "[0-9]")
  split_idx[is.na(split_idx)] <- 1
  drug_conc_cols_list <- split(drug_conc_cols, split_idx)
  drug_conc_cols_list <- lapply(drug_conc_cols_list, function(x) {
    names(x) <- gsub("[0-9]", "", names(x))
    x
  })
  
  idfs_len <- vapply(drug_conc_cols_list, length, FUN.VALUE = numeric(1))
  
  if (any(idfs_len != 2)) {
    df_metadata[, (intersect(unlist(gDRutils::get_env_identifiers(
      paste0(c("drug", "drug_name", "drug_moa", "concentration"),
             names(drug_conc_cols_list[which(idfs_len != 2)])),
      simplify = FALSE)), names(df_metadata))) := NULL]
    drug_conc_cols_list <- drug_conc_cols_list[-which(idfs_len != 2)]
  }
  
  # clean up concentration fields
  for (i in drug_conc_cols_list) {
    conc_id <- i[["concentration"]]
    drug_id <- i[["drug"]]
    untrt_id <- gDRutils::get_env_identifiers("untreated_tag")
    
    df_metadata[[conc_id]] <- ifelse(
      is.na(df_metadata[[conc_id]]),
      0,
      df_metadata[[conc_id]]
    )
    
    # set all untreated to 0
    df_metadata[df_metadata[[drug_id]] %in% untrt_id, conc_id] <- 0
    
    df_metadata[[conc_id]] <- 
      10 ^ round(log10(as.numeric(df_metadata[[conc_id]])), 6)
    df_metadata[[drug_id]] <- ifelse(
      is.na(df_metadata[[drug_id]]) & df_metadata[[conc_id]] == 0,
      untrt_id[[1]],
      df_metadata[[drug_id]]
    )
  }
  
  if (!gDRutils::get_env_identifiers("drug_name") %in% names(df_metadata)) {
    df_metadata <- add_Drug_annotation(df_metadata)
  }
  df_metadata
}


#' Order_result_df
#'
#' Order a data.table with results
#'
#' @param df_ a data.table with results
#'
#' @keywords utils
#' @return a ordered data.table with results
#'
order_result_df <- function(df_) {

  # Assertions:
  
  checkmate::assert_data_table(df_)

  ordered_1 <- gDRutils::get_header("ordered_1")
  ordered_2 <- gDRutils::get_header("ordered_2")
  cols <- c(
    ordered_1,
    setdiff(colnames(df_), c(ordered_1, ordered_2)),
    ordered_2
  )
  
  cols <- intersect(cols, colnames(df_))

  row_order_col <-
    unlist(intersect(
      c(
        gDRutils::get_env_identifiers(
          c("cellline_name", "duration", "drug_name"), 
          simplify = FALSE
        ),
        "Concentration",
        paste0(
          c(
            paste0(gDRutils::get_env_identifiers("drug_name"), "_"), 
            "Concentration_"
          ),
          sort(rep(2:10, 2))
        ),
        setdiff(colnames(df_), c(ordered_1, ordered_2))
      ),
      cols
    ))
  
  cols <- unlist(cols)
  data.table::setorderv(df_, row_order_col)
  df_ <- df_[, cols, with = FALSE]

  return(df_)
}

#' Detect model of data
#'
#' @param x data.table with raw data or SummarizedExperiment object 
#'          with gDR assays
#' 
#' @examples 
#' data_model("single-agent")
#' 
#' @return string with the information of the raw data follows single-agent or 
#' combination data model
#' @keywords utils
#' @export
data_model <- function(x) {
  UseMethod("data_model")
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
    stop("Unknown data model")
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
  checkmate::assert_subset(x, gDRutils::get_supported_experiments())
  
  exp_v <- gDRutils::get_experiment_groups()
  names(exp_v[grep(x, exp_v)])
}


#' Validate availability of data models
#' 
#' @param d_types character vector with experiment names in \code{MultiAssayExperiment} object
#' @param s_d_models character vector with names of supported experiment
#'
#' @keywords internal
validate_data_models_availability <- function(d_types, s_d_models) {
  checkmate::assert_character(d_types)
  checkmate::assert_character(s_d_models, null.ok = TRUE)
  
  dm_v <- gDRutils::get_experiment_groups()

  req_d_models <-
    unique(names(dm_v)[vapply(dm_v, function(x) any(d_types %in% x), logical(1))])
  f_models <- req_d_models[!req_d_models %in% s_d_models]
  if (length(f_models)) {
    msg1 <-
      sprintf(
        "'nested_identifiers_l' lacks information for the 
        following data model(s): '%s'",
        toString(f_models)
      )
    stop(msg1)
  }
}

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
    c("single-agent", "combination"), 
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
    data_model, c("single-agent", "combination"),
    null.ok = TRUE
  )
 
  ml <- list(`single-agent` = .get_default_single_agent_nested_identifiers(se),
             combination = .get_default_combination_nested_identifiers(se))
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

#' @keywords utils
.catch_warnings <- function(x) {
  warn  <- unlist(unique(x$warning))
  if (!is.null(warn)) {
    futile.logger::flog.warn(paste0(warn, collapse = "\n"))
  }
}

#' @keywords internal
#' @noRd
rbindParallelList <- function(x, name) {
  S4Vectors::DataFrame(
    do.call(
      rbind, 
      c(lapply(x, function(x) {
        dt <- data.table::as.data.table("[[" (x, name))
        data.table::setorder(dt)
        dt
      }), fill = TRUE)
    )
  )
}

#' Value Matching
#' 
#' Returns a lookup table or list of the positions of ALL matches of its first
#' argument in its second and vice versa. Similar to \code{\link{match}}, though
#' that function only returns the first match.
#' 
#' This behavior can be imitated by using joins to create lookup tables, but
#' \code{matches} is simpler and faster: usually faster than the best joins in
#' other packages and thousands of times faster than the built in
#' \code{\link{merge}}.
#' 
#' \code{all.x/all.y} correspond to the four types of database joins in the
#' following way:
#' 
#' \describe{ \item{left}{\code{all.x=TRUE}, \code{all.y=FALSE}} 
#' \item{right}{\code{all.x=FALSE}, \code{all.y=TRUE}} 
#' \item{inner}{\code{all.x=FALSE}, \code{all.y=FALSE}} 
#' \item{full}{\code{all.x=TRUE}, \code{all.y=TRUE}} }
#' 
#' Note that \code{NA} values will match other \code{NA} values.
#' 
#' @param x vector.  The values to be matched.  Long vectors are not currently
#'   supported.
#' @param y vector.  The values to be matched.  Long vectors are not currently
#'   supported.
#' @param all.x logical; if \code{TRUE}, then each value in \code{x} will be
#'   included even if it has no matching values in \code{y}
#' @param all.y logical; if \code{TRUE}, then each value in \code{y} will be
#'   included even if it has no matching values in \code{x}
#' @param list logical.  If \code{TRUE}, the result will be returned as a list
#'   of vectors, each vector being the matching values in y. If \code{FALSE},
#'   result is returned as a data.table with repeated values for each match.
#' @param indexes logical.  Whether to return the indices of the matches or the
#'   actual values.
#' @param nomatch the value to be returned in the case when no match is found.
#'   If not provided and \code{indexes=TRUE}, items with no match will be
#'   represented as \code{NA}.  If set to \code{NULL}, items with no match will
#'   be set to an index value of \code{length+1}.  If \code{indexes=FALSE}, they will
#'   default to \code{NA}.
#' @details Source of the function: 
#' https://github.com/cran/grr/blob/master/R/grr.R
#' 
#' @examples 
#' mat_elem <- data.table::data.table(
#'   DrugName = rep(c("untreated", "drugA", "drugB", "untreated"), 2),
#'   DrugName_2 = rep(c("untreated", "vehicle", "drugA", "drugB"), 2),
#'   clid = rep(c("C1", "C2"), each = 4)
#' )
#' untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")
#' ref_idx <- which(
#'   mat_elem$DrugName %in% untreated_tag |
#'    mat_elem$DrugName_2 %in% untreated_tag
#' )
#' ref <- mat_elem[ref_idx, ]
#' treated <- mat_elem[-ref_idx, ]
#' valid <- c("DrugName", "DrugName_2")
#' trt <- lapply(valid, function(x) {
#'   colnames <- c("clid", x) 
#'   treated[, colnames, with = FALSE]
#' })
#' trt <- do.call(paste, 
#'   do.call(rbind, lapply(trt, function(x) setNames(x, names(trt[[1]]))))
#' )
#' ref <- lapply(valid, function(x) {
#'   colnames <- c("clid", x) 
#'   ref[, colnames, with = FALSE]
#' })
#' ref <- do.call(paste, 
#'   do.call(rbind, lapply(ref, function(x) setNames(x, names(ref[[1]]))))
#' )
#' grr_matches(trt, ref, list = FALSE, all.y = FALSE)
#' 
#' @return data.table
#' 
#' @keywords utils
#' @export
grr_matches <- function(x, 
                    y, 
                    all.x = TRUE, 
                    all.y = TRUE, 
                    list = FALSE, 
                    indexes = TRUE, 
                    nomatch = NA) {
  result <- .Call("matches", x, y)
  result <- data.table::data.table(x = result[[1]], y = result[[2]])
  if (!all.y) {
    selected_rows <- which(result$x != length(x) + 1)
    result <- result[selected_rows, ]
  }
  if (!all.x) {
    selected_rows <- which(result$y != length(y) + 1)
    result <- result[selected_rows, ]
  }
  if (!indexes) {
    result$x <- x[result$x]
    result$y <- y[result$y]
  } else if (!is.null(nomatch)) {
    result$x[result$x == length(x) + 1] <- nomatch
    result$y[result$y == length(y) + 1] <- nomatch
  }
  if (list) {
    result <- tapply(result$y, result$x, function(z) z[!is.na(z)])
  }
  result
}


#' get info about created/present assays in SE at the given pipeline step
#' @param step string with pipeline step
#' @param data_model single-agent vs combination
#' @param status string return vector of assays created or present at the 
#' given step?
#' 
#' @keywords utils
#' @return assay
#' 
get_assays_per_pipeline_step <-
  function(step,
           data_model,
           status = c("created", "present")) {
    
    checkmate::assert_choice(step, get_pipeline_steps())
    checkmate::assert_choice(data_model, c("single-agent", "combination"))
    checkmate::assert_choice(status, c("created", "present"))
    
    fit_se <- if (data_model == "single-agent") {
      "Metrics"
    } else {
      c(
        "BlissExcess",
        "BlissScore",
        "CIScore_50",
        "CIScore_80",
        "HSAExcess",
        "HSAScore",
        "isobolograms",
        "Metrics",
        "SmoothMatrix"
      )
    }
    as_map <- list(
      "create_SE" = c("RawTreated", "Controls"),
      "normalize_SE" = "Normalized",
      "average_SE" = "Averaged",
      "fit_SE" = fit_se
    )
    
    ass <- if (status == "created") {
      as_map[[step]]
    } else {
      as.character(unlist(as_map[seq_len(which(names(as_map) == step))]))
    }
    ass
  }

#' add intermediate data (qs files) for given ma
#' @param mae mae with dose-response data
#' @param data_dir output directory
#' @param steps character vector with pipeline steps for which 
#'              intermediate data should be saved
#' 
#' @return \code{NULL}
#' 
#' @keywords internal
#' 
add_intermediate_data <- function(mae, data_dir, steps = get_pipeline_steps()) {
  
  checkmate::assert_class(mae, "MultiAssayExperiment")
  checkmate::assert_directory(data_dir, access = "rw")
  checkmate::assert_subset(steps, get_pipeline_steps())
  
  for (data_type in names(mae)) {
    for (step in steps) {
      as_names <-
        get_assays_per_pipeline_step(
          step, 
          data_model(data_type), 
          status = "present"
        )
      se <- mae[[data_type]]
      se_subset <- SummarizedExperiment::SummarizedExperiment(
        assays = SummarizedExperiment::assays(se)[as_names],
        rowData = SummarizedExperiment::rowData(se),
        colData = SummarizedExperiment::colData(se),
        metadata = S4Vectors::metadata(se)
      )
      save_intermediate_data(data_dir, step, data_type, se_subset)
    }
  }
}

#' get mae dataset from intermediate data
#' 
#' @param data_dir directory with intermediate data
#' 
#' @return MAE object
#' 
#' @keywords internal
#' 
get_mae_from_intermediate_data <- function(data_dir) {
  
  checkmate::assert_directory(data_dir)
  
  last_step <- tail(get_pipeline_steps(), n = 1)
  s_pattern <- paste0("__", last_step, ".qs")
  
  fpaths <- list.files(data_dir, pattern = s_pattern, full.names = TRUE)
  checkmate::assert_true(length(fpaths) > 0)
  
  sel <- list()
  
  for (fpath in fpaths) {
    exp_name <- sub(s_pattern, "", basename(fpath))
    sel[[exp_name]] <- qs::qread(fpath)
  }
  MultiAssayExperiment::MultiAssayExperiment(experiments = sel)
}
