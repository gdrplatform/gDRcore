#' Run drug response processing pipeline
#'
#' Run different components of the gDR drug response processing pipeline.
#' Either: create a SummarizedExperiment and normalize raw treated and control 
#' data (create_and_normalize_SE), average data (average_SE), or fit the 
#' processed data (fit_SE). See details for more in-depth explanations.
#'
#' @param x data.table of MAE with drug response data
#' @param df_ data.table of raw drug response data containing both treated and 
#' untreated values.
#' @param data_type single-agent vs combination
#' @param se \code{SummarizedExperiment} object.
#' @param readout string of the name containing the cell viability readout 
#' values.
#' @param control_mean_fxn function indicating how to average controls.
#' Defaults to \code{mean(x, trim = 0.25)}.
#' @param nested_identifiers_l list with the nested_identifiers(character v
#' ectors) for `single-agent` and (optionally) for `combination` data
#' @param nested_identifiers character vector with the nested_identifiers
#' for the given SE with a given data_type
#' @param nested_confounders Character vector of the nested_confounders for a 
#' given assay. nested_keys is character vector of column names to include in 
#' the data.tables in the assays of the resulting \code{SummarizedExperiment} 
#' object. Defaults to the \code{nested_identifiers} and 
#' \code{nested_confounders} if passed through \code{create_and_normalize_SE} 
#' or \code{runDrugResponseProcessingPipeline}.
#' @param series_identifiers character vector of identifiers in 
#' \code{measured} or \code{metric} which define a unique data point.
#' @param override_untrt_controls named list containing defining factors in 
#' the treatments. Defaults to \code{NULL}.
#' @param override_masked boolean indicating whether or not to override 
#' the masked wells in the averaging and include all wells. 
#' Defaults to \code{FALSE}.
#' @param ndigit_rounding integer indicating number of digits to round to 
#' in calculations. Defaults to \code{4}.
#' @param n_point_cutoff integer of how many points should be considered the 
#' minimum required to try to fit a curve. Defaults to \code{4}.
#' @param control_assay string containing the name of the assay representing 
#' the controls in the \code{se}. Defaults to \code{"Controls"}.
#' @param raw_treated_assay string containing the name of the assay 
#' representing the raw treated data in the \code{se}.
#' Defaults to \code{"RawTreated"}.
#' @param normalized_assay string of the assay name containing the 
#' normalized data. Defaults to \code{"Normalized"}.
#' @param averaged_assay string of the name of the averaged assay in the
#' \linkS4class{SummarizedExperiment}. Defaults to \code{"Averaged"}.
#' @param metrics_assay string of the name of the metrics assay to output
#' in the returned \linkS4class{SummarizedExperiment}
#' Defaults to \code{"Metrics"}.
#' @param range_conc vector of concetrations range values.
#' @param force_fit boolean indicating whether or not to force the fit.
#' @param pcutoff numeric cutoff value.
#' @param cap numeric value representing the value to cap the highest allowed 
#' relative viability at.
#' @param curve_type vector of curve type values.
#' @param data_dir string with the path to the directory with intermediate data 
#' of experiments (qs files). If set to NULL (default) intermediate data is not 
#' saved/read in.
#' @param partial_run logical flag indicating if the pipeline should be run 
#' partially (from the step defined with `start_from`)
#' @param start_from string indicating the pipeline step from which partial 
#' run should be launched
#' @param selected_experiments character vector with experiments for which 
#' pipeline should be run. This option works only for the pipeline being run 
#' partially (i.e. with `partial_run` flag set to `TRUE`)
#'
#' @details
#' \code{runDrugResponseProcessingPipeline} is made up of 3 separate steps:
#' \itemize{
#'  \item{"create_and_normalize_SE"}
#'  \item{"average_SE"}
#'  \item{"fit_SE"}
#'}
#'
#' For create_and_normalize_SE, this creates a SummarizedExperiment object 
#' from a data.table, where the data.table contains treatments on rows, and 
#' conditions on columns. 
#' A \linkS4class{SummarizedExperiment} object containing two asssays is 
#' created: treated readouts will live in an assay called \code{"RawTreated"},
#' and reference readouts live in an assay called \code{"Controls"}. 
#' Subsequently, the treated and control elements will be normalized to output 
#' two metrics: 
#'
#' For average_SE, take the normalized assay and average the nested 
#' \code{DataFrame}s across unique\code{nested_identifiers}. 
#'
#' For fit_SE, take the averaged assay and fit curves to obtain metrics, one 
#' set of metrics for each normalization type set.
#' 
#' Pipeline can be run partially with `partial_run` flag set to TRUE. The 
#' `start_from` string defines the step from which the pipeline will be 
#' launched. However, partial run of the pipeline is possible only if the whole
#' pipeline was launched at least once with defined `data_dir` and intermediate 
#' data was saved as qs files into `data_dir`. 
#' 
#' Pipeline can be run for the selected experiments by changing the default 
#' value of `selected_experiments` param. This scenario only works when 
#' `partial_run` is enabled.
#'
#' @name runDrugResponseProcessingPipelineFxns
#' 
#' @examples 
#' p_dir <- file.path(tempdir(), "pcheck")
#' dir.create(p_dir) 
#' td <- gDRimport::get_test_data()
#' l_tbl <- gDRimport::load_data(gDRimport::manifest_path(td), gDRimport::template_path(td), gDRimport::result_path(td))
#' imported_data <- merge_data(
#'   l_tbl$manifest, 
#'   l_tbl$treatments, 
#'   l_tbl$data
#' )
#' runDrugResponseProcessingPipeline(
#'   imported_data, 
#'   data_dir = p_dir
#' )
#' 
#' @return MAE object
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
                                    nested_confounders = intersect(
                                      names(df_), 
                                      gDRutils::get_env_identifiers("barcode")
                                    ), 
                                    override_untrt_controls = NULL,
                                    ndigit_rounding = 4,
                                    control_assay = "Controls",
                                    raw_treated_assay = "RawTreated",
                                    normalized_assay = "Normalized") {
  
  checkmate::assert_data_table(df_)
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
    nested_identifiers = nested_identifiers,
    nested_confounders = nested_confounders,
    override_untrt_controls = override_untrt_controls)
  se <- normalize_SE(se = se, 
    data_type = data_type, 
    nested_identifiers = nested_identifiers,
    nested_confounders = nested_confounders,
    control_mean_fxn = control_mean_fxn,
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
                                              nested_identifiers_l = NULL,
                                              nested_confounders = 
                                                gDRutils::get_env_identifiers(
                                                  "barcode"
                                                ),
                                              override_untrt_controls = NULL,
                                              override_masked = FALSE,
                                              ndigit_rounding = 4,
                                              n_point_cutoff = 4,
                                              control_assay = "Controls",
                                              raw_treated_assay = "RawTreated",
                                              normalized_assay = "Normalized",
                                              averaged_assay = "Averaged",
                                              metrics_assay = "Metrics",
                                              data_dir = NULL,
                                              partial_run = FALSE,
                                              start_from = 
                                                get_pipeline_steps()[1],
                                              selected_experiments = NULL) {
  
  checkmate::assert_multi_class(x, c("data.table", "MultiAssayExperiment"))
  if (inherits(x, "data.table")) {
    checkmate::assert_true(any(
      gDRutils::get_env_identifiers("untreated_tag") %in%
        x[[gDRutils::get_env_identifiers("drug")]]
    ))
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
    stop(
      "Path for/to the intermediate data is required with partial_run enabled"
    )
  } 

  # when running pipeline with x = MAE the identifiers from MAE's metadata 
  # should be restored
  if (inherits(x, "MultiAssayExperiment")) {
    m_idfs <- S4Vectors::metadata(x[[1]])[["identifiers"]]
    for (idx in seq_along(m_idfs)) {
      gDRutils::set_env_identifier(names(m_idfs)[idx], m_idfs[[idx]])
    }
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
      sel[[experiment]] <- read_intermediate_data(
        path = data_dir, 
        step = tail(get_pipeline_steps(), n = 1), 
        experiment = experiment
      )

    } else {
    
      message("Processing ", experiment)
      data_type <- experiment
      nested_identifiers <- inl$nested_identifiers_l[[data_model(data_type)]]

      run_vars <- list(
        partial_run = partial_run,
        start_from = start_from,
        exp = experiment,
        data_dir = data_dir
      )
    
      # 1st step - Create SE
      se <- run_pipeline_step(
        run_vars = run_vars,
        step = "create_SE", 
        se = se,
        step_fun = create_SE_step, 
        step_args = list(
          inl = inl, 
          exp = experiment,
          data_type = data_type, 
          readout = readout, 
          nested_identifiers = nested_identifiers,
          override_untrt_controls = override_untrt_controls
          )
      )

      # 2nd step - Normalize SE
      se <- run_pipeline_step(
        run_vars = run_vars,
        step = "normalize_SE", 
        se = se,
        step_fun = normalize_SE, 
        step_args = list(
          se = se$result,
          data_type = data_type,
          nested_identifiers = nested_identifiers,
          nested_confounders = inl$nested_confounders,
          control_mean_fxn = control_mean_fxn,
          control_assay = control_assay,
          raw_treated_assay = raw_treated_assay,
          normalized_assay = normalized_assay,
          ndigit_rounding = ndigit_rounding
        )
      )

      # 3rd step - Average SE
      se <- run_pipeline_step(
        run_vars = run_vars,
        step = "average_SE", 
        se = se,
        step_fun = average_SE, 
        step_args = list(
          se = se$result, 
          data_type = data_type,
          series_identifiers = nested_identifiers,
          override_masked = override_masked, 
          normalized_assay = normalized_assay, 
          averaged_assay = averaged_assay
        ), 
        if_paste_warnings = TRUE
      )
    
      # 4th step - Fit SE
      if (data_type == "matrix") {
        step_args <- list(
          se = se$result,
          data_type = data_type,
          series_identifiers = nested_identifiers,
          averaged_assay = averaged_assay
        )
        step_fun <- fit_SE.combinations
      } else {
        step_args <- list(
          se = se$result,
          data_type = data_type,
          averaged_assay = averaged_assay,
          metrics_assay = metrics_assay, 
          n_point_cutoff = n_point_cutoff
        )
        step_fun <- fit_SE
      }
      se <- run_pipeline_step(
        run_vars = run_vars,
        step = "fit_SE", 
        se = se,
        step_fun = step_fun, 
        step_args = step_args, 
        if_read_intermediate_data = FALSE
      )
    
      paste_warnings(se$warnings)
      sel[[experiment]] <- se$result
      se <- list()
    }
   }

  gDRutils::reset_env_identifiers()
  MultiAssayExperiment::MultiAssayExperiment(experiments = sel)
}
