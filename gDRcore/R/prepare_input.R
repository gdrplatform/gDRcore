#' Prepare input data common for all experiments
#'
#' Current steps
#' - refining nested confounders
#' - refining nested identifiers
#' - splitting df_ into (per experiment) df_list
#' 
#' @param x data.frame with raw data or MAE object with dose-reponse data
#' @param ... additional parameters
#' 
#' @return list of input data
#' 
#' @export
prepare_input <-
  function(x, ...) {
    UseMethod("prepare_input")
  }

#' Prepare input data common for all experiments
#'
#'  Current steps
#' - refining nested confounders
#' - refining nested identifiers
#' - splitting df_ into (per experiment) df_list
#' @param x data.frame with raw data
#' @param ... additional parameters
#' @param nested_identifiers_l list with the nested_identifiers(character vectors) 
#' for `single-agent` and (optionally) for `combination` data
#' @param nested_confounders Character vector of the nested_confounders for a given assay.
#' nested_keys is character vector of column names to include in the data.frames
#' in the assays of the resulting \code{SummarizedExperiment} object.
#' Defaults to the \code{nested_identifiers} and \code{nested_confounders} if passed through
#' 
#' @return list of input data
#' 
#' @export
prepare_input.data.frame <-
  function(x,
           nested_confounders = gDRutils::get_env_identifiers("barcode"),
           nested_identifiers_l = .get_default_nested_identifiers(),
           ...) {
    
    checkmate::assert_data_frame(x, min.rows = 1, min.cols = 1)
    checkmate::assert_character(nested_confounders, null.ok = TRUE)
    checkmate::assert_list(nested_identifiers_l, null.ok = TRUE)
    
    inl <- list(
      df_ = NULL,
      df_list = NULL,
      nested_confounders = NULL
    )
    
    inl$df_ <- identify_data_type(x)
    inl$df_list <- split_raw_data(inl$df_)
    inl$nested_identifiers_l <- .set_nested_identifiers(nested_identifiers_l)
    inl$nested_confounders <- .set_nested_confounders(
      nested_confounders = nested_confounders,
      df = x
    )
    inl$exps <- .set_exps(inl$df_list)
    
    validate_data_models_availability(names(inl$df_list), names(inl$nested_identifiers_l))
    
    inl
  }

#' Prepare input data common for all experiments
#'
#' Current steps
#' - refining nested confounders
#' - refining nested identifiers
#' - splitting df_ into (per experiment) df_list
#' @param x MAE object with dose-reponse data
#' @param ... additional parameters
#' @param nested_identifiers_l list with the nested_identifiers(character vectors) 
#' for `single-agent` and (optionally) for `combination` data
#' @param nested_confounders Character vector of the nested_confounders for a given assay.
#' nested_keys is character vector of column names to include in the data.frames
#' in the assays of the resulting \code{SummarizedExperiment} object.
#' Defaults to the \code{nested_identifiers} and \code{nested_confounders} if passed through
#' @param raw_data_field metadata field with raw data
#' @param split_data Boolean indicating need of splitting the data into experiment types
#' 
#' @return list of input data
#' 
#' @export
prepare_input.MultiAssayExperiment <-
  function(x,
           nested_confounders = gDRutils::get_SE_identifiers(x[[1]], "barcode"),
           nested_identifiers_l = .get_default_nested_identifiers(x[[1]]),
           raw_data_field = "experiment_raw_data",
           split_data = TRUE,
           ...) {
    
    checkmate::assert_true(inherits(x, "MultiAssayExperiment"))
    checkmate::assert_character(nested_confounders, null.ok = TRUE)
    checkmate::assert_list(nested_identifiers_l, null.ok = TRUE)
    
    inl <- list(
      df_list = NULL,
      nested_confounders = NULL,
      nested_identifiers_l = NULL
    )
    
    inl$df_list <-
      lapply(names(x), function(y) {
        md <- S4Vectors::metadata(x[[y]])
        if (is.null(md[[raw_data_field]])) {
          NULL
        }
        md[[raw_data_field]]
      })
    
    if (split_data) {
      inl$df_ <- lapply(inl$df_list, identify_data_type)
      if ("matrix" %in% names(x)) {
        inl$df_ <- inl$df_[grep("single-agent",
                                names(x), invert = TRUE)]
      }
      inl$df_list <- split_raw_data(unique(plyr::rbind.fill(inl$df_)))
    } else {
      names(inl$df_list) <- names(x)
    }
    
    inl$nested_confounders <- .set_nested_confounders(
      nested_confounders = nested_confounders,
      df = inl$df_list[[1]]
    )
    inl$nested_identifiers_l <- .set_nested_identifiers(
      nested_identifiers_l,
      args = list(se = x[[1]])
    )
    inl$exps <- .set_exps(x)
    
    validate_data_models_availability(names(x), names(inl$nested_identifiers_l))
    
    inl
  }

.set_nested_confounders <- function(nested_confounders, df) {
  x_names <- names(df)
  # Some experiment can have nested_confounders = NULL that is appropriate situation for
  # internal data
  if (!is.null(nested_confounders) &&
      any(!nested_confounders %in% x_names)) {
    
    confounders_intersect <- intersect(
      c(nested_confounders, gDRutils::get_env_identifiers("barcode")), 
      x_names
    )
    
    warning(
      sprintf(
        "'%s' nested confounder(s) is/are not present in the data.
    Switching into '%s' nested confounder(s).",
        setdiff(nested_confounders, x_names),
        confounders_intersect
      )
    )
    
    if (length(confounders_intersect) > 0) {
      confounders_intersect
    }
  } else {
    nested_confounders
  }
}

.set_nested_identifiers <- function(nested_identifiers_l, args = list()) {
  if (is.null(nested_identifiers_l)) {
    do.call(.get_default_nested_identifiers, args)
  } else {
    nested_identifiers_l
  }
}

.set_exps <- function(df_list) {
  df_names <- names(df_list)
  exps <- lapply(df_names, function(x) NULL)
  names(exps) <- df_names
  
  exps
}
