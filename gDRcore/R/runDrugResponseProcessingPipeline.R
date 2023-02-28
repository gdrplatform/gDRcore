#' Run drug response processing pipeline
#'
#' Run different components of the gDR drug response processing pipeline.
#' Either: create a SummarizedExperiment and normalize raw treated and control data (create_and_normalize_SE),
#' average data (average_SE), or fit the processed data (fit_SE). See details for more in-depth explanations.
#'
#' @param df_ data.frame of raw drug response data containing both treated and untreated values.
#' @param se \code{SummarizedExperiment} object.
#' @param readout string of the name containing the cell viability readout values.
#' @param control_mean_fxn function indicating how to average controls.
#' Defaults to \code{mean(x, trim = 0.25)}.
#' @param nested_identifiers_l list with the nested_identifiers(character vectors) 
#' for `single-agent` and (optionally) for `combination` data
#' @param nested_identifiers character vector with the nested_identifiers
#' for the given SE with a given data_type
#' @param nested_confounders Character vector of the nested_confounders for a given assay.
#' nested_keys is character vector of column names to include in the data.frames
#' in the assays of the resulting \code{SummarizedExperiment} object.
#' Defaults to the \code{nested_identifiers} and \code{nested_confounders} if passed through
#' \code{create_and_normalize_SE} or \code{runDrugResponseProcessingPipeline}.
#' @param series_identifiers character vector of identifiers in \code{measured} or \code{metric}
#' which define a unique data point.
#' @param override_untrt_controls named list containing defining factors in the treatments.
#' Defaults to \code{NULL}.
#' @param override_masked boolean indicating whether or not to override the masked wells
#' in the averaging and include all wells. 
#' Defaults to \code{FALSE}.
#' @param ndigit_rounding integer indicating number of digits to round to in calculations.
#' Defaults to \code{4}.
#' @param n_point_cutoff integer of how many points should be considered the minimum required to try to fit a curve.
#' Defaults to \code{4}.
#' @param control_assay string containing the name of the assay representing the controls in the \code{se}.
#' Defaults to \code{"Controls"}.
#' @param raw_treated_assay string containing the name of the assay representing the raw treated data in the \code{se}.
#' Defaults to \code{"RawTreated"}.
#' @param normalized_assay string of the assay name containing the normalized data.
#' Defaults to \code{"Normalized"}.
#' @param averaged_assay string of the name of the averaged assay in the \linkS4class{SummarizedExperiment}.
#' Defaults to \code{"Averaged"}.
#' @param metrics_assay string of the name of the metrics assay to output
#' in the returned \linkS4class{SummarizedExperiment}
#' Defaults to \code{"Metrics"}.
#' @param add_raw_data  boolean indicating whether or not to include raw data into experiment metadata.
#' @param range_conc vector of concetrations range values.
#' @param force_fit boolean indicating whether or not to force the fit.
#' @param pcutoff numeric cutoff value.
#' @param cap numeric value representing the value to cap the highest allowed relative viability at.
#' @param curve_type vector of curve type values.
#' @param data_dir string with the path to the directory with intermediate data of experiments (qs files).
#' If set to NULL (default) intermediate data is not saved/read in.
#' @param partial_run logical flag indicating if the pipeline should be run partially 
#' (from the step defined with `start_from`)
#' @param start_from string indicating the pipeline step from which partial run should be launched
#' @param selected_experiments character vector with experiments for which pipeline should be run.
#' This option works only for the pipeline being run partially (i.e. with `partial_run` flag set to `TRUE`)
#'
#' @details
#' \code{runDrugResponseProcessingPipeline} is made up of 3 separate steps:
#' \itemize{
#'  \item{"create_and_normalize_SE"}
#'  \item{"average_SE"}
#'  \item{"fit_SE"}
#'}
#'
#' For create_and_normalize_SE, this creates a SummarizedExperiment object from a data.frame, 
#' where the data.frame contains treatments on rows, and conditions on columns. 
#' A \linkS4class{SummarizedExperiment} object containing two asssays is created:
#' treated readouts will live in an assay called \code{"RawTreated"},
#' and reference readouts live in an assay called \code{"Controls"}. Subsequently, the treated
#' and control elements will be normalized to output two metrics: 
#'
#' For average_SE, take the normalized assay and average the nested \code{DataFrame}s across unique
#' \code{nested_identifiers}. 
#'
#' For fit_SE, take the averaged assay and fit curves to obtain metrics, one set of metrics for each
#' normalization type set.
#' 
#' Pipeline can be run partially with `partial_run` flag set to TRUE. The `start_from` string defines the step 
#' from which the pipeline will be launched. However, partial run of the pipeline is possible only if the whole
#' pipeline was launched at least once with defined `data_dir` and intermediate data was saved as qs files 
#' into `data_dir`. 
#' 
#' Pipeline can be run for the selected experiments by changing the default value of `selected_experiments` param`.
#' This scenario only works when `partial_run` is enabled.
#'
#' @name runDrugResponseProcessingPipelineFxns
#'
NULL


#' @rdname runDrugResponseProcessingPipelineFxns
#' @export
create_and_normalize_SE <- function(df_,
                                    data_type,
                                    readout = "ReadoutValue",
                                    control_mean_fxn = function(x) {
                                      mean(x, trim = 0.25)
                                    },
                                    nested_identifiers = NULL,
                                    nested_confounders = gDRutils::get_env_identifiers("barcode"), 
                                    override_untrt_controls = NULL,
                                    ndigit_rounding = 4,
                                    control_assay = "Controls",
                                    raw_treated_assay = "RawTreated",
                                    normalized_assay = "Normalized") {
  
  checkmate::assert_data_frame(df_)
  checkmate::assert_string(data_type)
  checkmate::assert_string(readout)
  checkmate::assert_function(control_mean_fxn)
  checkmate::assert_character(nested_identifiers, null.ok = TRUE)
  checkmate::assert_character(nested_confounders, null.ok = TRUE)
  checkmate::assert_vector(override_untrt_controls, null.ok = TRUE)
  checkmate::assert_numeric(ndigit_rounding)
  checkmate::assert_string(control_assay)
  checkmate::assert_string(raw_treated_assay)
  checkmate::assert_string(normalized_assay)
 
  se <- create_SE(df_ = df_, 
    data_type = data_type, 
    readout = readout, 
    control_mean_fxn = control_mean_fxn, 
    nested_identifiers = nested_identifiers,
    nested_confounders = nested_confounders,
    override_untrt_controls = override_untrt_controls)
  se <- normalize_SE(se = se, 
    data_type = data_type, 
    nested_identifiers = nested_identifiers,
    nested_confounders = nested_confounders,
    control_assay = control_assay, 
    raw_treated_assay = raw_treated_assay, 
    normalized_assay = normalized_assay,
    ndigit_rounding = ndigit_rounding)
  se
}


#' @rdname runDrugResponseProcessingPipelineFxns 
#' @export
runDrugResponseProcessingPipeline <- function(x,
                                              readout = "ReadoutValue",
                                              control_mean_fxn = function(x) {
                                                mean(x, trim = 0.25)
                                              },
                                              nested_identifiers_l = .get_default_nested_identifiers(),
                                              nested_confounders = gDRutils::get_env_identifiers("barcode"),
                                              override_untrt_controls = NULL,
                                              override_masked = FALSE,
                                              ndigit_rounding = 4,
                                              n_point_cutoff = 4,
                                              control_assay = "Controls",
                                              raw_treated_assay = "RawTreated",
                                              normalized_assay = "Normalized",
                                              averaged_assay = "Averaged",
                                              metrics_assay = "Metrics",
                                              add_raw_data = TRUE,
                                              data_dir = NULL,
                                              partial_run = FALSE,
                                              start_from = get_pipeline_steps()[1],
                                              selected_experiments = NULL) {
  
  checkmate::assert_multi_class(x, c("data.frame", "MultiAssayExperiment"))
  if (inherits(x, "data.frame")) {
    checkmate::assert_true(any(gDRutils::get_env_identifiers("untreated_tag") %in%
                                 x[[gDRutils::get_env_identifiers("drug")]]))
  }
  checkmate::assert_string(readout)
  checkmate::assert_function(control_mean_fxn)
  checkmate::assert_multi_class(nested_identifiers_l, c("list"), null.ok = TRUE)
  checkmate::assert_character(nested_confounders, null.ok = TRUE)
  checkmate::assert_vector(override_untrt_controls, null.ok = TRUE)
  checkmate::assert_flag(override_masked)
  checkmate::assert_numeric(ndigit_rounding)
  checkmate::assert_number(n_point_cutoff)
  checkmate::assert_string(control_assay)
  checkmate::assert_string(raw_treated_assay)
  checkmate::assert_string(normalized_assay)
  checkmate::assert_string(averaged_assay)
  checkmate::assert_string(metrics_assay)
  checkmate::assert_flag(add_raw_data)
  if (!is.null(data_dir)) {
    checkmate::assert_directory(data_dir)
  } 
  checkmate::assert_flag(partial_run)
  checkmate::assert_string(start_from)
  checkmate::assert_choice(start_from, get_pipeline_steps())
  checkmate::assert_character(selected_experiments, null.ok = TRUE)
  checkmate::assert_subset(selected_experiments,
                           names(gDRutils::get_experiment_groups()),
  )
  
  if (!is.null(selected_experiments) && !partial_run) {
    stop("Selected experiments are only supported with partial_run enabled")
  } 
  
  if (is.null(data_dir) && partial_run) {
    stop("Path for/to the intermediate data is required with partial_run enabled")
  } 
 
  inl <- prepare_input(x, nested_confounders, nested_identifiers_l)
  
  # sel - list with all experiments data
  # se - list with single experiment data 
  se <- list()
  sel <- list()
  
  for (experiment in names(inl$df_list)) {
    
    # skip pipeline for not selected experiments (optionally)
    # in detail: read intermediate data from the last step
    if (!is.null(selected_experiments) &&
        !(experiment %in% selected_experiments) &&
        partial_run) {
      
      message("Reading intermediate data for ", experiment)
    sel[[experiment]] <- read_intermediate_data(data_dir, tail(get_pipeline_steps(), n = 1), experiment)
      
    } else {
    
    message("Processing ", experiment)
    data_type <- experiment
    nested_identifiers <- inl$nested_identifiers_l[[data_model(data_type)]]
    
    if (!partial_run || !do_skip_step("create_SE", start_from)) {
      if (is.null(inl$df_list[[experiment]])) {
        
        msg1 <- sprintf("It's impossible to run pipeline from the first step ('%s') for experiment: '%s'. ",
                get_pipeline_steps()[1], experiment)
        msg2 <- "The pipeline has been run with 'add_raw_data' flag disabled? "
        msg3 <- sprintf("Consider running the pipeline from the second ('%s') step", get_pipeline_steps()[2])
        stop(c(msg1, msg2, msg3))
      }
    se <- purrr::quietly(create_SE)(df_ = inl$df_list[[experiment]],
                                                     data_type = data_type,
                                                     readout = readout,
                                                     control_mean_fxn = control_mean_fxn,
                                                     nested_identifiers = nested_identifiers,
                                                     nested_confounders = inl$nested_confounders,
                                                     override_untrt_controls = override_untrt_controls)
    
    if (add_raw_data) {
      se$result <-
        gDRutils::set_SE_experiment_raw_data(se$result, inl$df_list[[experiment]])
    }
    
      if (!is.null(data_dir)) {
        save_intermediate_data(data_dir, "create_SE", experiment, se$result)
      }
    } else {
       if (is_preceding_step("create_SE", start_from)) {
        se$result <- read_intermediate_data(data_dir, "create_SE", experiment)
       }
    }
    
    if (!partial_run || !do_skip_step("normalize_SE", start_from)) {
      se <- purrr::quietly(normalize_SE)(
        se = se$result,
        data_type = data_type,
        nested_identifiers = nested_identifiers,
        nested_confounders = inl$nested_confounders,
        control_assay = control_assay,
        raw_treated_assay = raw_treated_assay,
        normalized_assay = normalized_assay,
        ndigit_rounding = ndigit_rounding
      )
      
      if (!is.null(data_dir)) {
        save_intermediate_data(data_dir, "normalize_SE", experiment, se$result)
      }
    } else {
       if (is_preceding_step("normalize_SE", start_from)) {
        se$result <- read_intermediate_data(data_dir, "normalize_SE", experiment)
       }
    }
    
    if (!partial_run || !do_skip_step("average_SE", start_from)) {
    
      paste_warnings(se$warnings)
      se <- purrr::quietly(average_SE)(se = se$result, 
                                        data_type = data_type,
                                        series_identifiers = nested_identifiers,
                                        override_masked = override_masked, 
                                        normalized_assay = normalized_assay, 
                                        averaged_assay = averaged_assay)
      if (!is.null(data_dir)) {
        save_intermediate_data(data_dir, "average_SE", experiment, se$result)
      }
      paste_warnings(se$warnings)
    } else {
       if (is_preceding_step("average_SE", start_from)) {
        se$result <- read_intermediate_data(data_dir, "average_SE", experiment)
       }
    }
    
    if (!partial_run || !do_skip_step("fit_SE", start_from)) {
      se <- if (data_type == "matrix") {
        purrr::quietly(fit_SE.combinations)(se = se$result,
                            data_type = data_type,
                            series_identifiers = nested_identifiers,
                            averaged_assay = averaged_assay)
        
      } else {
        purrr::quietly(fit_SE)(se = se$result, 
               data_type = data_type,
               nested_identifiers = nested_identifiers,
               averaged_assay = averaged_assay, 
               metrics_assay = metrics_assay, 
               n_point_cutoff = n_point_cutoff)
      }
      if (!is.null(data_dir)) {
        save_intermediate_data(data_dir, "fit_SE", experiment, se$result)
      }
    }
    
    paste_warnings(se$warnings)
    sel[[experiment]] <- se$result
    se <- list()
    }
  }
  
  MultiAssayExperiment::MultiAssayExperiment(experiments = sel)
}


#' get pipeline steps
#' 
#' @keywords internal
get_pipeline_steps <-
  function() {
    c("create_SE",
      "normalize_SE",
      "average_SE",
      "fit_SE")
  }

#' check if the given step can be skipped if partial run is chosen
#' 
#' @param current_step, string with the step to be evaluated
#' @param start_from string indicating the pipeline step from which partial run should be launched
#' @param steps charvect with all available steps
#' 
#' @keywords internal
#' @return logical
#' 
do_skip_step <-
  function(current_step, start_from, steps = get_pipeline_steps()) {
    
    checkmate::assert_choice(current_step, steps)
    checkmate::assert_choice(start_from, steps, null.ok = TRUE)
    checkmate::assert_character(steps)
    
    if (is.null(start_from)) {
      FALSE
    } else {
      c_idx <- which(steps %in% current_step)
      s_idx <- which(steps %in% start_from)
      c_idx < s_idx
    }
  }

#' check if the given step is preceding the step chosen in the partial run
#' 
#' @param current_step, string with the step to be evaluated
#' @param start_from string indicating the pipeline step from which partial run should be launched
#' @param steps charvect with all available steps
#' 
#' @keywords internal
#' @return logical
#' 
is_preceding_step <-
  function(current_step, start_from, steps = get_pipeline_steps()) {
   
    checkmate::assert_choice(current_step, steps)
    checkmate::assert_choice(start_from, steps)
    checkmate::assert_character(steps)
    
    c_idx <- which(steps %in% current_step)
    s_idx <- which(steps %in% start_from)
    s_idx - c_idx == 1
  }

#' save intermediate data for the given experiment and step to qs file
#' 
#' @param path string with the save directory for the qs file 
#' @param step, string with the step name
#' @param experiment string with the experiment name
#' @param se output se 
#' 
#' @keywords internal
#' 
save_intermediate_data <- function(path, step, experiment, se) {
  
  checkmate::assert_directory(path, "rw")
  checkmate::assert_string(step)
  checkmate::assert_string(experiment)
  
  fpath <- file.path(path, paste0(experiment, "__", step, ".qs"))
  qs::qsave(se, fpath)
}

#' read intermediate data for the given experiment and step to qs file
#' 
#' @param path string with the input directory of the qs file 
#' @param step, string with the step name
#' @param experiment string with the experiment name
#' 
#' @keywords internal
#' @return se
#' 
read_intermediate_data <- function(path, step, experiment) {
  
  checkmate::assert_directory(path, "r")
  checkmate::assert_string(step)
  checkmate::assert_string(experiment)
  
  fpath <- file.path(path, paste0(experiment, "__", step, ".qs"))
  qs::qread(fpath)
}

#' @keywords internal
paste_warnings <- function(list, sep = "\n") {
  warning(paste0(list, sep = sep), call. = FALSE)
}


#' Prepare input data common for all experiments
#'
#' Current steps
#' - refining nested confounders
#' - refining nested identifiers
#' - splitting df_ into (per experiment) df_list
#' 
#' @param x data.frame with raw data or MAE object with dose-reponse data
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
#' 
#' @export
prepare_input.data.frame <-
  function(x,
           nested_confounders = gDRutils::get_env_identifiers("barcode"),
           nested_identifiers_l = .get_default_nested_identifiers()) {
    
    checkmate::assert_data_frame(x, min.rows = 1, min.cols = 1)
    checkmate::assert_character(nested_confounders, null.ok = TRUE)
    checkmate::assert_list(nested_identifiers_l, null.ok = TRUE)
    
    
    inl <- list(
      df_ = NULL,
      df_list = NULL,
      nested_confounders = NULL,
      nested_identifiers_l = nested_identifiers_l
    )
    
    nested_confounders <- if (!is.null(nested_confounders) &&
                              any(!nested_confounders %in% names(x))) {
      warning(
        sprintf(
          "'%s' nested confounder(s) is/are not present in the data.
    Switching into '%s' nested confounder(s).",
          setdiff(nested_confounders, names(x)),
          intersect(nested_confounders, names(x))
        )
      )
      confounders_intersect <-
        intersect(nested_confounders, names(x))
      if (length(confounders_intersect) == 0) {
        NULL
      } else {
        confounders_intersect
      }
    } else {
      nested_confounders
    }
    if (!is.null(nested_confounders)) {
      inl$nested_confounders <- nested_confounders
    }
    
    inl$df_ <- identify_data_type(x)
    inl$df_list <- split_raw_data(inl$df_)
    
    validate_data_models_availability(names(inl$df_list), names(nested_identifiers_l))
    
    inl$exps <- lapply(names(inl$df_list), function(x) {
      NULL
    })
    names(inl$exps) <- names(inl$df_list)
    
    inl
  }

#' Prepare input data common for all experiments
#'
#' Current steps
#' - refining nested confounders
#' - refining nested identifiers
#' - splitting df_ into (per experiment) df_list
#' @param x MAE object with dose-reponse data
#' 
#' @export
prepare_input.MultiAssayExperiment <-
  function(x,
           nested_confounders = gDRutils::get_SE_identifiers(x[[1]], "barcode"),
           nested_identifiers_l = .get_default_nested_identifiers(x[[1]]),
           raw_data_field = "experiment_raw_data") {
    
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
        md <- metadata(x[[y]])
        if (is.null(md[[raw_data_field]])) {
          NULL
        }
         md[[raw_data_field]]
      })
    names(inl$df_list) <- names(x)
    
    nested_confounders <- if (!is.null(nested_confounders) &&
                              any(!nested_confounders %in% names(inl$df_list[[1]]))) {
      warning(
        sprintf(
          "'%s' nested confounder(s) is/are not present in the data.
    Switching into '%s' nested confounder(s).",
          setdiff(nested_confounders, names(inl$df_list[[1]])),
          intersect(nested_confounders, names(inl$df_list[[1]]))
        )
      )
      confounders_intersect <-
        intersect(nested_confounders, names(inl$df_list[[1]]))
      if (length(confounders_intersect) == 0) {
        NULL
      } else {
        confounders_intersect
      }
    } else {
      nested_confounders
    }
    if (!is.null(nested_confounders)) {
      inl$nested_confounders <- nested_confounders
    }
    
    inl$nested_confounders <- if (is.null(nested_confounders)) {
      gDRutils::get_SE_identifiers(x[[1]], "barcode")
    } else {
      nested_confounders
    }
    
    inl$nested_identifiers_l <- if (is.null(nested_identifiers_l)) {
      .get_default_nested_identifiers(x[[1]])
    } else {
      nested_identifiers_l
    }
    
    
    validate_data_models_availability(names(x), names(nested_identifiers_l))
    
    inl$exps <- lapply(names(x), function(y) {
      NULL
    })
    names(inl$exps) <- names(x)
    
    inl
  }