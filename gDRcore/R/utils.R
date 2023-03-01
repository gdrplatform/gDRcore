#' @noRd
#' @keywords internal
#'
.create_mapping_factors <- function(rowdata, coldata) {
  mapping_cols <- c(colnames(rowdata), colnames(coldata))
  mapping_entries <- expand.grid(row_id = rownames(rowdata), col_id = rownames(coldata), stringsAsFactors = FALSE)

  mapping_entries <- base::merge(mapping_entries, rowdata, by.x = "row_id", by.y = 0, all.x = TRUE)
  mapping_entries <- base::merge(mapping_entries, coldata, by.x = "col_id", by.y = 0, all.x = TRUE)

  mapping_entries
}


#' cleanup_metadata
#'
#' Cleanup a dataframe with metadata
#'
#' @param df_metadata a dataframe with metadata
#'
#' @return a dataframe with cleaned metadata
#' @details Adds annotations and check whether user provided correct input data.
#' @export
#'
cleanup_metadata <- function(df_metadata) {

  # Assertions:
  stopifnot(inherits(df_metadata, "data.frame"))

  # Round duration to 6 values. 
  df_metadata[[gDRutils::get_env_identifiers("duration")]] <-
    round(as.numeric(df_metadata[[gDRutils::get_env_identifiers("duration")]], 6))

  df_metadata <- add_CellLine_annotation(df_metadata)
  
  drug_conc_cols <- unlist(gDRutils::get_env_identifiers(c("drug", "drug2", "drug3",
                                                 "concentration", "concentration2", "concentration3"),
                                                    simplify = FALSE))
  drug_conc_cols <- drug_conc_cols[drug_conc_cols %in% names(df_metadata)]
  split_idx <- stringr::str_extract(names(drug_conc_cols), "[0-9]")
  split_idx[is.na(split_idx)] <- 1
  drug_conc_cols_list <- split(drug_conc_cols, split_idx)
  drug_conc_cols_list <- lapply(drug_conc_cols_list, function(x) {
    names(x) <- gsub("[0-9]", "", names(x))
    x
    })
  
  # clean up concentration fields
  for (i in drug_conc_cols_list) {
    df_metadata[[i[["concentration"]]]] <- ifelse(is.na(df_metadata[[i[["concentration"]]]]),
                                                  0,
                                                  df_metadata[[i[["concentration"]]]])
    df_metadata[df_metadata[[i[["drug"]]]] %in%
                  gDRutils::get_env_identifiers("untreated_tag"), i[["concentration"]]] <- 0 # set all untreated to 0
    df_metadata[[i[["concentration"]]]] <- 10 ^ round(log10(as.numeric(df_metadata[[i[["concentration"]]]])), 6)
    df_metadata[[i[["drug"]]]] <- ifelse(is.na(df_metadata[[i[["drug"]]]]) & df_metadata[[i[["concentration"]]]] == 0,
           gDRutils::get_env_identifiers("untreated_tag")[[1]],
           df_metadata[[i[["drug"]]]])
  }
  
  df_metadata <- add_Drug_annotation(df_metadata)
  df_metadata
  return(df_metadata)
}


#' Order_result_df
#'
#' Order a dataframe with results
#'
#' @param df_ a dataframe with results
#'
#' @return a ordered dataframe with results
#' @export
#'
Order_result_df <- function(df_) {

  # Assertions:
  stopifnot(inherits(df_, "data.frame"))

  cols <- c(gDRutils::get_header("ordered_1"),
            setdiff(colnames(df_),
                    c(
                      gDRutils::get_header("ordered_1"), gDRutils::get_header("ordered_2")
                    )),
            gDRutils::get_header("ordered_2"))
  cols <- intersect(cols, colnames(df_))

  row_order_col <-
    intersect(
      c(
        gDRutils::get_env_identifiers(c("cellline_name", "duration", "drug_name"), simplify = FALSE),
        "Concentration",
        paste0(c(
          paste0(gDRutils::get_env_identifiers("drug_name"), "_"), "Concentration_"
        ),
        sort(rep(2:10, 2))),
        setdiff(colnames(df_), c(
          gDRutils::get_header("ordered_1"), gDRutils::get_header("ordered_2")
        ))
      ),
      cols
    )

  df_ <- df_[do.call(order, df_[, unlist(row_order_col)]), unlist(cols)]

  return(df_)
}

#' Detect model of data
#'
#' @param x data.frame with raw data or SummarizedExperiment object with gDR assays
#' 
#' @return string with the information of the raw data follows single-agent or combination data model
#' @export
data_model <- function(x) {
  UseMethod("data_model")
}


#' Detect model of data in data.frame
#'
#' @param x data.frame of raw drug response data containing both treated and untreated values. 
#'
#' @return string with the information of the raw data follows single-agent or combination data model
#' @export
data_model.data.frame <- function(x) {
  
  checkmate::assert_data_frame(x)
  drug_ids <- unlist(gDRutils::get_env_identifiers(c("drug_name", "drug_name2"), simplify = FALSE))
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


#'
#' 
get_data_type_to_data_model_mapping <- function() {
  c(
    `single-agent` = "single-agent",
    "cotreatment" = "single-agent",
    "co-dilution" = "single-agent",
    "matrix" = "combination"
  )
}

#' Detect model of data from experiment name
#'
#' @param x character with experiment name
#'
#' @return string with the information of the raw data follows single-agent or combination data model
#' @export
data_model.character <- function(x) {
  # TODO: switch to gDRutils::get_experiment_groups()
  # once we clean-up single-agent/combination assignemnts
  exp_v <- get_data_type_to_data_model_mapping()
  
  checkmate::assert_subset(x, names(exp_v))
  
  exp_v[[x]]
  
}

validate_data_models_availability <- function(d_types, s_d_models) {
  
  dm_v <- get_data_type_to_data_model_mapping()
  
  req_d_models <-
    unique(as.character(dm_v[names(dm_v) %in% d_types]))
  f_models <- req_d_models[!req_d_models %in% s_d_models]
  if (length(f_models)) {
    msg1 <-
      sprintf(
        "'nested_identifiers_l' lacks information for the following data model(s): '%s'",
        toString(f_models)
      )
    stop(msg1)
  }
}

#' Get default nested identifiers
#'
#' @param x data.frame with raw data or SummarizedExperiment object with gDR assays
#' @param data_model single-agent vs combination
#' @return vector of nested identifiers
#' @export
get_default_nested_identifiers <- function(x, data_model = NULL) {
  UseMethod("get_default_nested_identifiers")
}


#' @export
#' @rdname get_default_nested_identifiers
get_default_nested_identifiers.data.frame <- function(x, data_model = NULL) {
  
  checkmate::assert_data_frame(x)
  checkmate::assert_choice(data_model, c("single-agent", "combination"), null.ok = TRUE)
  
  .get_default_nested_identifiers(se = NULL, data_model = NULL)
  
}

#' @export
#' @rdname get_default_nested_identifiers
get_default_nested_identifiers.SummarizedExperiment <- function(x,
                                                                data_model = NULL) {
  .get_default_nested_identifiers(se, data_model)
}

.get_default_nested_identifiers <- function(se = NULL, data_model = NULL) {
 
  checkmate::assert_choice(data_model, c("single-agent", "combination"), null.ok = TRUE)
 
  ml <- list(`single-agent` = .get_default_single_agent_nested_identifiers(se),
             combination = .get_default_combination_nested_identifiers(se))
  if (is.null(data_model)) {
    ml
  } else {
    ml[[data_model]]
  }
  
}

.get_default_single_agent_nested_identifiers <- function(se = NULL) {
  if (is.null(se)) {
    gDRutils::get_env_identifiers("concentration")
  } else {
    gDRutils::get_SE_identifiers(se, "concentration")
  }
}

.get_default_combination_nested_identifiers <- function(se = NULL) {
  if (is.null(se)) {
    as.character(unlist(gDRutils::get_env_identifiers(c("concentration", "concentration2"),
                                         simplify = FALSE)))
  } else {
    as.character(unlist(gDRutils::get_SE_identifiers(se, c("concentration", "concentration2"),
                                        simplify = FALSE)))
  }
}

.catch_warnings <- function(x) {
  warn  <- unlist(unique(x$warning))
  if (!is.null(warn)) {
    futile.logger::flog.warn(paste0(warn, collapse = "\n"))
  }
}

#' Round concentration to ndigit significant digits
#'
#' @param x value to be rounded.
#' @param ndigit number of significant digits (default = 4).
#'
#' @return rounded x
#' @export
round_concentration <- function(x, ndigit = 3) {
  round(10 ^ (round(log10(x), ndigit)), ndigit - 1 - floor(log10(x)))
}


#' @keywords internal
#' @noRd
rbindList <- function(x) {
  S4Vectors::DataFrame(do.call(plyr::rbind.fill, x))
}

#' @keywords internal
#' @noRd
rbindParallelList <- function(x, name) {
  S4Vectors::DataFrame(do.call(plyr::rbind.fill, lapply(x, function(x) as.data.frame("[[" (x, name)))))
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
#'   result is returned as a data frame with repeated values for each match.
#' @param indexes logical.  Whether to return the indices of the matches or the
#'   actual values.
#' @param nomatch the value to be returned in the case when no match is found.
#'   If not provided and \code{indexes=TRUE}, items with no match will be
#'   represented as \code{NA}.  If set to \code{NULL}, items with no match will
#'   be set to an index value of \code{length+1}.  If {indexes=FALSE}, they will
#'   default to \code{NA}.
#' @details Source of the function: https://github.com/cran/grr/blob/master/R/grr.R
#' @export
matches <- function(x, y, all.x = TRUE, all.y = TRUE, list = FALSE, indexes = TRUE, nomatch = NA) {
  result <- .Call("matches", x, y)
  result <- data.frame(x = result[[1]], y = result[[2]])
  if (!all.y) {
    result <- result[result$x != length(x) + 1, ]
  }
  if (!all.x) {
    result <- result[result$y != length(y) + 1, ]
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
#' @param status string return vector of assays created or present at the given step?
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
#' @param data_dir  output directory
#' @param steps charvec with pipeline steps for which intermediate data should be saved
#' 
#' @export
add_intermediate_data <- function(mae, data_dir, steps = get_pipeline_steps()) {
  
  checkmate::assert_class(mae, "MultiAssayExperiment")
  checkmate::assert_directory(data_dir, access = "rw")
  checkmate::assert_subset(steps, get_pipeline_steps())
  
  for (data_type in names(mae)) {
    for (step in steps) {
      as_names <-
        get_assays_per_pipeline_step(step, data_model(data_type), status = "present")
      se <- mae[[data_type]]
      se_subset <- SummarizedExperiment(
        assays = assays(se)[as_names],
        rowData = rowData(se),
        colData = colData(se),
        metadata = metadata(se)
      )
      save_intermediate_data(data_dir, step, data_type, se_subset)
    }
  }
}

#' get mae dataset from intermediate data
#' 
#' @param data_dir directory with intermediate data
#' 
#' @export
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