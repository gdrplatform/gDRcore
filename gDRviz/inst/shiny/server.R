server <- function(input, output, session) {

  # load config data (default and optionally custom ones)
  CONFIG <- gDRcomponents::get_config(dump = TRUE)

  # set proper gdr_username
  gDRcomponents::set_gdr_username(request = session$request)

  S_LOG_PREFIX <- "GV::SESSION"
  SS_LOG_PREFIX <- "GV::SESSION_STARTED"
  SE_LOG_PREFIX <- "GV::SESSION_ENDED"
  logExtra::mark_timestamp(S_LOG_PREFIX)
  logExtra::info(SS_LOG_PREFIX, sprintf("user:%s", gDRcomponents::get_gdr_username()))

  supported_exp_names <- gDRutils::get_supported_experiments()
  sa_exp_name <- gDRutils::get_supported_experiments("sa")
  combo_exp_name <- gDRutils::get_supported_experiments("combo")

  #-##### spinner
  gDRcomponents::show_spinner()
  #-#####


  #-##### global reactiveValues
  # Define reactiveValues to communicate across the modules

  rv <- reactiveValues(
    session = session,
    pidfs = reactiveVal(NULL),
    nonunique_additional_perturbation = reactiveVal(NULL),
    data_list = reactive(NULL),
    triggerViz = reactiveVal(NULL),
    response_metrics = reactive(NULL),
    response_data = reactive(NULL)
  )

  output$restart <- renderUI({
    a("RESTART", href = paste0("javascript:history.go(0)"))
  })

  #-#### front page
  gDRcomponents::frontPageSERVER("frontpage", parent_session = session)

  #-####

  #-#### logs
  callModule(gDRcomponents::moduleDisplayLogs, "displayLogs", filepath = LOG_FILE)
  #-####

  #-##### admin mode
  # to prevent password guessing, we introduce 3 attempts lines per session
  remaining_attempts <- reactiveVal(3)

  output$login_admin <- renderUI({
    if (remaining_attempts() > 0) {
      shinyWidgets::actionBttn(
        "admin_login",
        style = "material-circle",
        color = "default",
        icon = icon("lock")
      )
    } else {
      shinyWidgets::actionBttn(
        "admin_login",
        style = "material-circle",
        color = "danger",
        icon = icon("lock")
      )

    }
  })

  observeEvent(input$admin_login, {

    if (remaining_attempts() > 0) {
      showModal(
        modalDialog(
          title = "Admin area - login required",
          "Enter password:",
          fluidRow(
            column(3, passwordInput("password", label = NULL)),
            column(3, actionButton("password_button", "OK"))
          ),
          shinyjs::hidden(
            span(id = "failed_login_span",
                 textOutput("failed_login"), style = "color:red"
            )),
          easyClose = TRUE
        )
      )
    } else {
      showModal(
        modalDialog(
          title = "Admin area - access restricted",
          easyClose = TRUE
        )
      )
    }
  })

  output$failed_login <- renderText(sprintf("remaining attempts: %i", remaining_attempts()))

  observeEvent(input$password_button, {
    if (input$password == Sys.getenv("ADMIN_PASSWORD", "gdrp4sw0rd")) {
      warning(remaining_attempts())
      initAdminMode(input, CONFIG$plugins)
      showNotification("Successful login.", type = "message")
      removeModal()
    } else {
      warning(remaining_attempts())
      shinyjs::show("failed_login_span")
      remaining_attempts(remaining_attempts() - 1)
    }
    if (remaining_attempts() == 0) {
      showNotification("Login failed", type = "warning")
      removeModal()
    }
  })
  #-#####

  #-##### additional dataset (present/absent) with logic based on envs defined in global.R
  # TODO create a module that support multiple tabs/project-views
  #   and/or upgrade while switching to elasticsearch-based solution
  observe({
    rv$new_tab_filtration <- NEW_TAB_FILTRATION
    if (any(c(NEW_TAB_DESC, NEW_TAB_VALUE, NEW_TAB_FILTRATION) == "")) {
      shinyjs::addClass(selector = sprintf("#tabs > li > a[data-value='%s']", NEW_TAB_VALUE), class = "hidden")
    } else {
      shinyjs::removeClass(selector = sprintf("#tabs > li > a[data-value='%s']", NEW_TAB_VALUE), class = "hidden")
    }
  })
  #-#####

  # # launch via url
  rv$selected_url_params <- eventReactive(session$clientData$url_search, {

    url <- session$clientData$url_search

    url_params <- lapply(parseQueryString(url), function(x) {
      stringr::str_replace_all(
        x,
        c("\"" = "")
      )
    })

    if (!is.null(url_params$in_app) && url_params$in_app == "gDRsearch2") {
      return(url_params)
    }

    req(nchar(url) > 0)

    id_params <- get_url_params_as_list(
      url,
      id = get_url_param(),
      regex = TRUE)

    if (is.null(id_params$ds_id)) {
      # Assumption: valid url must have at least "ds_id"
      showNotification("The given URL is invalid.", type = "warning")
    }

    req(id_params$ds_id)
    verified_project <-
      check_project_availability(get_project_tables(use_local_cache = TRUE),
                                 id_params$ds_id,
                                 id_params$meta_id)

    if (is.null(unlist(verified_project))) {
      showNotification("Selected project is not available.", type = "warning")
      out <- NULL
    } else {
      out <- id_params # all selected param
      out$ds_id <- verified_project$ds_id # update with verified ds_id
      out$meta_id <- verified_project$meta_id # update with verified meta_id
    }

    out
  })

  observeEvent(rv$selected_url_params(), {
    req(rv$selected_url_params())

    if (all(!is.null(rv$selected_url_params()$in_app),
            rv$selected_url_params()$in_app == "gDRsearch2")) {
      shinyjs::show("in_plugins", asis = TRUE)
      shinyjs::delay(500, updateTabItems(session, "main_tabs", "gDRsearch"))
    } else if (all(!is.null(unlist(rv$selected_url_params())))) {
      shinyjs::show("in_plugins", asis = TRUE)
      shinyjs::delay(500, updateTabItems(session, "main_tabs", "gDRexplorer"))
    }
  })

  #-##### too-big dataset
  # TODO: rewrite the logic for this part
  # the current logic for enabling/disabling some modules based on the dataset size
  # is overly complex and outdated (prepared for single SE, not valid for MAE with potentially three experiments)
  issue_with_mae <- reactiveValues(
    too_big_mae = FALSE
  )
  observeEvent(rv$error, {
    issue_with_mae$too_big_mae <- TRUE
    updateTabItems(session, "main_tabs", "gDRsearch")
  })

  observeEvent(input$main_tabs, {
    if (input$main_tabs %in% c("gDRin", "gDRsearch", "gDRsearch2")) {
      gDRutils::reset_env_identifiers()
    }
  }, ignoreInit = TRUE)

  #-#####

  #-##### load dataset from any data source/(aka tab )
  observeEvent(input$main_tabs, {
    if (input$main_tabs == "frontpage") {
      shinyjs::hide(id = "miniSummary")
    } else {
      shinyjs::show(id = "miniSummary")
    }

    #TODO: we can add modal to tell user about possibility to lost all progress after switch input plugin
    if (input$main_tabs == "gDRsearch") {
      rv$mae_source <- "gDRsearch"
    } else if (input$main_tabs == "gDRin") {
      rv$mae_source <- "gDRin"
    } else if (input$main_tabs == "gDRsearch2") {
      rv$mae_source <- "gDRsearch2"
    }
  })


  rv$id_data <- reactive({
    req(rv$pidfs())
    id <- gDRutils::prettify_flat_metrics(gDRutils::get_header("id"), human_readable = TRUE)
    unlist(c(rv$pidfs()[c("drug", "drug2", "cellline")], id))
  })

  .ds_raw <- reactive({
    s_dtype <- c("MultiAssayExperiment")

    my_ds <- if (is.null(rv$data_src)) {
      NULL
    } else if (rv$data_src == "gDRsearch") {
      req(rv$mae_search())
    } else if (rv$data_src == "gDRin") {
      rv$mae
    }
    req(my_ds)
    if (inherits(my_ds, "MultiAssayExperiment")) {
      gDRutils::reset_env_identifiers()
      mae_idfs <- gDRutils::get_MAE_identifiers(my_ds)
      invisible(gDRutils::update_env_idfs_from_mae(mae_idfs))
      my_ds
    } else {
      stop(sprintf("Unsupported dataset type. Currently supported: '%s'", toString(s_dtype)))
    }
  })
  #-#####


  #-##### load summary module
  # TODO: refactor for MAE dataset (show data for all experiments not only the first one)
  callModule(gDRcomponents::miniSummary, id = "miniSummary", ds = .ds)

  #-##### validate and transform metadata of the dataset
  # TODO: get rid of 'drug_name' temporary fix (once dataset are updated)
  # TODO: wrap the code into functions/reactives
  # - validate MAE
  # - disambiguate cell-line/drugs data
  # - substitute 'drugname' to 'drug_name' column in metadataa
  # - check if 'too-big-mae' condition applies
  .ds <- reactive({
    futile.logger::flog.trace("Main App: \t disambiguating", name = "trace.logger")

    mae <- req(.ds_raw())
    # checking if MAE identifiers are the same env
    mae_idfs <- gDRutils::get_MAE_identifiers(mae)
    gDRutils::update_env_idfs_from_mae(mae_idfs)

    rv$pidfs <- reactiveVal(gDRutils::get_prettified_identifiers(simplify = TRUE))

    # validate MAE
    v <- tryCatch(gDRutils::validate_MAE(mae),
                  error = function(e) {
                    e
                  })
    # show popup if validation failed
    if (!is.null(unlist(unname(v)))) {
      futile.logger::flog.trace(sprintf("Main App: \t 'validate_mae' failed with '%s'", v$message),
                                name = "trace.logger")
      # inform the user that given dataset in not gDRviz compatible and can't be used
      shinyalert::shinyalert(
        sprintf(
          "Selected dataset is not compatible with gDRviz. Please select another one."
        ),
        confirmButtonCol = "#3c8dbc"
      )
      req(FALSE)
    }

    for (i in seq_along(mae)) {
      # renaming drugnameX identifiers to drug_nameX
      # this is a temporary fix until the identifier in the incoming data and testData is updated
      if (any(grepl("drugname", names(S4Vectors::metadata(mae[[i]])$identifiers)))) {
        names(S4Vectors::metadata(mae[[i]])$identifiers) <-
          gsub("drugname", "drug_name", names(metadata(mae[[i]])$identifiers))
      }

      # refine colData
      SummarizedExperiment::colData(mae[[i]]) <-
        gDRutils::refine_coldata(SummarizedExperiment::colData(mae[[i]]), mae[[i]])
      # refine rowData
      SummarizedExperiment::rowData(mae[[i]]) <-
        gDRutils::refine_rowdata(SummarizedExperiment::rowData(mae[[i]]), mae[[i]])
    }

    search_simple_max_cell_lines <- gDRcomponents::check_env_for_value(NULL, "search_simple_max_cell_lines")
    search_simple_max_drugs <- gDRcomponents::check_env_for_value(NULL, "search_simple_max_drugs")

    for (i in seq_along(mae)) {
      if (nrow(SummarizedExperiment::colData(mae[[i]])) < as.numeric(search_simple_max_cell_lines) /
          20 &&
          nrow(SummarizedExperiment::rowData(mae[[i]])) < as.numeric(search_simple_max_drugs) /
          20) {
        shinyjs::js$enableTab("shiny-tab-grid")
      } else {
        shinyjs::js$disableTab("shiny-tab-grid")
      }
    }

    issue_with_mae$too_big_mae  <- FALSE
    return(gDRutils::set_unique_identifiers(mae))
  })
  #-#####

  #-##### fit source
  rv$fit_source <- reactive({
    futile.logger::flog.trace("Main App: extracting fit sources", name = "trace.logger")
    data <- req(assay_metrics_initial())
    fit_sources <-
      lapply(data, function(x) {
        unique(na.omit(x$fit_source))
      })
    lapply(fit_sources, function(x) {
      checkmate::assert_true(length(x) > 0)
    })

    fs_v <- intersect(names(fit_sources), names(.ds()))
    if (length(fs_v)) {
      fit_sources[fs_v]
    } else {
      NULL
    }

  })

  #-##### transform dataset assays
  # -- for all assays (including combo-matrix):
  # - convert BumpyMatrice(s) to data.table(s)
  # -- for non-combo-matrix assays:
  # - dirty hack for ec50/c50 assay data
  # - flatten with normalization type/fit source
  # - prettify metrics
  # - merge assay-data with drug/cell_line data
  # - drop excessive columns (drug2-related if only untreated data found)
  # - run capVals (on metrics only)
  # TODO: split into helper function (see especialy repetitive code for non-matrix assays)
  # TODO: wrap the code into functions/reactives/shiny modules

  assay_metrics_initial <- reactive({
    futile.logger::flog.trace("Main App: extracting initial Metrics data", name = "trace.logger")
    ds <- req(.ds())
    shinybusy::show_modal_spinner(spin = "orbit", text = "Loading the project")
    gDRutils::MAEpply(
      ds,
      gDRutils::convert_se_assay_to_custom_dt,
      assay_name = "Metrics",
      output_table = "Metrics_initial"
    )
  })


  assay_metrics_raw <- reactive({
    futile.logger::flog.trace("Main App: extracting Metrics data", name = "trace.logger")
    req(assay_metrics_initial())
    ds <- req(.ds())

    data <-
      gDRutils::MAEpply(
        ds,
        gDRutils::convert_se_assay_to_custom_dt,
        assay_name = "Metrics",
        output_table = "Metrics_raw"
      )
    gDRcomponents::remove_spinner()
    return(data)
  })

  assay_metrics <- reactive({
    futile.logger::flog.trace("Main App: \t capping values", name = "trace.logger")
    assay_metrics <- assay_metrics_raw()
    lapply(assay_metrics, gDRutils::capVals)
  })

  assay_normalized <- reactive({
    futile.logger::flog.trace("Main App: extracting Normalized data", name = "trace.logger")
    ds <- req(.ds())

    gDRutils::MAEpply(
      ds,
      gDRutils::convert_se_assay_to_custom_dt,
      assay_name = "Normalized"
    )
  })

  assay_averaged <- reactive({
    futile.logger::flog.trace("Main App: extracting Averaged data", name = "trace.logger")
    ds <- req(.ds())

    futile.logger::flog.trace("Main App: \t extracting", name = "trace.logger")

    gDRutils::MAEpply(
      ds,
      gDRutils::convert_se_assay_to_custom_dt,
      assay_name = "Averaged"
    )
  })

  combo_dts <- reactive({
    ds <- req(.ds())
    combo_dt <- NULL
    for (i in seq_along(ds)) {
      if (isTRUE(gDRutils::is_combo_data(ds[[i]]))) {
        combo_dt <- gDRutils::convert_combo_data_to_dt(ds[[i]])
      }
    }
    combo_dt
  })

  # create object that stores all assays
  assay_objects <- reactive({
    req(.ds())
    rv$triggerViz(FALSE)
    futile.logger::flog.trace("Main App: building list of assays", name = "trace.logger")
    newObjects <- list()
    for (i in seq_along(.ds())) {
      newObjects[[i]] <- list(
        metrics_raw = assay_metrics_raw()[[i]],
        metrics = assay_metrics()[[i]],
        normalized = assay_normalized()[[i]],
        averaged = assay_averaged()[[i]]
      )
      if (!is.null(combo_dts()) && names(.ds())[i] == combo_exp_name) {
        futile.logger::flog.trace("Main App: appending combo data", name = "trace.logger")
        newObjects[[i]] <- c(newObjects[[i]], combo_dts())
      }
    }
    names(newObjects) <- names(.ds())
    newObjects
  })

  rv$assay_object_manage_data <- reactive({

    futile.logger::flog.trace("Main App: preparing assay object for managing data", name = "trace.logger")
    rv$triggerViz(FALSE)
    req(rv$data_src)

    # TODO: extend the logic while working on the new combo model: GDR-2244
    if (rv$data_src == "gDRsearch2") {
      data_l <- list()
      # gDRsearch2 returns ready to use response_metrics & response_data
      # but these values can be also NULL (if user selects combination data only)
      # this is a new edge-case not possible with data from gDRin/gDRsearch
      if (is.reactive(rv$response_data2)) {
        data_l[[sa_exp_name]] <- list(
          metrics = rv$response_metrics2(),
          averaged = rv$response_data2()
        )
      }

      # Show combination tab only when combo exp was created, and have any data
      if (is.reactive(rv$combo2) && NROW(rv$combo2()$excess) != 0) {
        data_l[[combo_exp_name]] <- rv$combo2()
        # get rid of the "_combo" postfixes in list names
        names(data_l[[combo_exp_name]]) <- gsub("_combo", "", names(data_l[[combo_exp_name]]))
      }
      data_l
    } else {
      # gDRin and gDRsearch return the same object
      req(assay_objects())
      assay_objects()
    }

  })

  #-##### disable grid module if too many drug vs cell-line combinations
  observeEvent(input$elementExist, {

    assay_metrics <- assay_objects()[["metrics"]]
    n_cell_lines <- length(unique(assay_metrics$clid))
    n_drugs <- length(unique(assay_metrics$Gnumber))

    search_simple_max_cell_lines <- as.numeric(gDRcomponents::check_env_for_value(NULL, "search_simple_max_cell_lines"))
    search_simple_max_drugs <- as.numeric(gDRcomponents::check_env_for_value(NULL, "search_simple_max_drugs"))


    if (n_cell_lines < search_simple_max_cell_lines / 20 &&
        n_drugs < search_simple_max_drugs / 20) {
      shinyjs::js$enableTab("shiny-tab-grid")
    } else {
      shinyjs::js$disableTab("shiny-tab-grid")
    }

    # handle tooltips below
    toltip_text <-
      sprintf(
        "Dose Response Overview is available for max %i drugs and %i cell lines",
        as.numeric(search_simple_max_drugs) / 20,
        as.numeric(search_simple_max_cell_lines) / 20
      )
    shinyjs::runjs(sprintf('$("#gDRvizTabs > ul > li:nth-child(6) > a").attr("title_dis", "%s");', toltip_text))

    shinyjs::runjs(sprintf('$("#gDRvizTabs > ul > li:nth-child(8) > a").attr("title_dis", "%s");',
                           "Drug Combo Treatments is available for combo data only"))
  })
  #-#####

  observeEvent(c(rv$data_list,
                 rv$gDRsearch2_triger), ignoreInit = FALSE, ignoreNULL = FALSE, {

                   plugins_table <- data.table::rbindlist(CONFIG$plugins, fill = TRUE, idcol = TRUE)
                   plugins_table$.id <- gDRcomponents:::plugin_make_id(parent_id = NULL, plugins_table$.id)

                   combo_plugins <- plugins_table[plugins_table$type == "output - combo"]$.id
                   sa_plugins <- plugins_table[plugins_table$type == "output - sa"]$.id


                   # user's just started using the app no data selected yet
                   if ((is.null(rv$data_list()$response_data) ||
                        is.null(rv$data_list()$response_metrics)) &&
                       is.null(rv$data_list()$excess)) {
                     # hiding/showing is analogous in every condition
                     # collapse tabs
                     shinyjs::hide("out_plugins")
                     shinyjs::hide("manage_plugins")
                     # hide tabs
                     shinyjs::js$hideSidebarTab("shiny-tab-out_plugins")
                     shinyjs::js$hideSidebarTab("shiny-tab-manage_plugins")
                     # hide subtab for each single-agent and combo plugin
                     lapply(combo_plugins, shinyjs::hide)
                     lapply(sa_plugins, shinyjs::hide)
                   } else {
                     # single-agent data is ready-to-use
                     if (!is.null(rv$data_list()$response_data) &&
                         !is.null(rv$data_list()$response_metrics)) {
                       shinyjs::js$showSidebarTab("shiny-tab-out_plugins")
                       shinyjs::show("out_plugins")
                       shinyjs::js$showSidebarTab("shiny-tab-manage_plugins")
                       shinyjs::show("manage_plugins")
                       lapply(sa_plugins, shinyjs::show)
                     }
                     # combination data is ready-to-use
                     if (!is.null(rv$data_list()$excess)) {
                       shinyjs::js$showSidebarTab("shiny-tab-out_plugins")
                       shinyjs::show("out_plugins")
                       shinyjs::js$showSidebarTab("shiny-tab-manage_plugins")
                       shinyjs::show("manage_plugins")
                       lapply(combo_plugins, shinyjs::show)
                     }
                   }
                 })


  observeEvent(rv$assay_object_manage_data(), {
    req(rv$assay_object_manage_data())
    req(rv$id_data())
    id_data <- rv$id_data()

    # is single-agent data available?
    if (!is.null(rv$assay_object_manage_data()[[sa_exp_name]][["metrics"]])) {
      metrics <- rv$assay_object_manage_data()[[sa_exp_name]][["metrics"]]
      averaged <- rv$assay_object_manage_data()[[sa_exp_name]][["averaged"]]
      id_data_intersect <- Reduce(intersect, list(id_data, names(metrics), names(averaged)))

      metrics[, (id_data_intersect) := NULL]
      averaged[, (id_data_intersect) := NULL]

      rv$response_metrics <- reactiveVal(metrics)
      rv$response_data <- reactiveVal(averaged)
    } else {
      averaged <- NULL
      metrics <- NULL
    }

    # Condition for combo data. Firstly we check if we have any data source (that is whether
    # any data was generated) or when data is not from gDRsearch2 if there is experiment `combination`.
    # If not combo is NULL. Otherwise for gDRsearch2 we check rv$combo2 object for combo data,
    # and for other modules we extract combo data from SE
    combo <- if (is.null(rv$data_src)) {
      NULL
    } else if (rv$data_src == "gDRsearch2") {
      if (isTruthy(rv$combo2)) {
        req(rv$combo2())
      } else {
        NULL
      }
    } else {
      # gDRin and gDRsearch
      if (combo_exp_name %in% names(.ds())) {
        data <- gDRutils::convert_combo_data_to_dt(.ds()[[combo_exp_name]])
        out <- lapply(data, function(x) {
          x[, (id_data_intersect) := NULL]
        })
        out
      } else {
        NULL
      }
    }

    rv$combo_object <- reactive(combo)
    futile.logger::flog.trace("Main App: \t preparing data list", name = "trace.logger")
    rv$data_list <- reactive(c(list(response_data = averaged,
                                    response_metrics = metrics),
                               combo))
  })


  observeEvent(rv$data_list, ignoreInit = TRUE, {
    req(rv$data_list())
    futile.logger::flog.trace("Main App: \t checking additional metadata", name = "trace.logger")
    additional_metadata <- gDRutils::get_additional_variables(rv$data_list(), prettified = TRUE)
    futile.logger::flog.trace(sprintf("Main App: \t found additional metadata: %s",
                                      additional_metadata), name = "trace.logger")
    rv$nonunique_additional_perturbation <- reactiveVal(additional_metadata)
  })

  observeEvent(rv$nonunique_additional_perturbation(), ignoreInit = TRUE, {
    additional_metadata <- req(rv$nonunique_additional_perturbation())
    futile.logger::flog.trace("Main App: \t managing additional metadata", name = "trace.logger")
    showModal(modalDialog(
      shinyFeedback::useShinyFeedback(),
      withTags({
        lapply(rv$nonunique_additional_perturbation(), function(var) {
          button_id <- sprintf("button_%s", var)
          default_choice <- if (var == gDRutils::prettify_flat_metrics("source_id", human_readable = TRUE)) {
            "average"
          } else {
            NULL
          }
          radioButtons(
            button_id,
            sprintf("How to proceed with %s?", var),
            choiceNames = c(
              "combine with Drug",
              "combine with Cell Line",
              "average"
            ),
            choiceValues = c("toDrug", "toCellLine", "average"),
            selected = default_choice,
            inline = TRUE
          )
        })
      }),
      title = "Duplicates or additional perturbations have been identified",
      size = "l",
      easyClose = FALSE,
      footer = tagList(
        actionButton("merge_button", "OK")
      )
    ))

    observeEvent(input[["merge_button"]], ignoreInit = TRUE, {
      removeModal()
      futile.logger::flog.trace("Main App: \t processing additional metadata", name = "trace.logger")
      shinybusy::show_modal_spinner(spin = "orbit", text = "Processing additional perturbations")
      old_data <- isolate(rv$data_list())
      for (var in additional_metadata) {
        # get input ids
        button_id <- sprintf("button_%s", var)
        # capture the current state of the data objects
        # run modification; S3 dispatch decides the action
        for (data_name in names(old_data)) {
          futile.logger::flog.trace(sprintf("Main App: \t\t  processing %s as an additional variable in %s data",
                                            var, data_name),
                                    name = "trace.logger")
          if (!is.null(old_data[[data_name]])) {
            data <- gDRutils::addClass(old_data[[data_name]], newClass = var)
            if (!var %in% names(data)) {
              next
            }
            new_data <- gDRutils::modifyData(data,
                                             option = input[[button_id]],
                                             keep = NULL)
            # mark new the list of assays as modified
            old_data[[data_name]] <- gDRutils::addClass(new_data, "modifiedAssays")
          }
        }
      }

      old_combo <- rv$combo_object()
      for (data_name in names(old_data)) {
        if (data_name %in% c("response_data", "response_metrics") && !is.null(old_data[[data_name]])) {
          rv[[data_name]] <- reactiveVal(old_data[[data_name]])
        } else {
          old_combo[[data_name]] <- old_data[[data_name]]
        }
      }
      shinybusy::remove_modal_spinner()
      rv$combo_object <- reactive(old_combo)
    })
  })

  # jump to MetricClustering plugin if selected data is mixed single-agent/combination or single-agent only
  # jump to DrugComboTreatment plugin if selected data is combination only
  observeEvent(
    c(rv$response_metrics(), rv$response_data(), rv$combo_object),
    ignoreInit = TRUE,
    {
      # either single-agent or combination data is available
      req((
        isTruthy(rv$response_data()) &&
          isTruthy(rv$response_metrics())
      ) || isTruthy(rv$combo_object()))

      req(rv$pidfs())
      futile.logger::flog.trace("Main App: \t checking uniqueness of vars", name = "trace.logger")

      if (isTruthy(rv$response_data()) &&
          isTruthy(rv$response_metrics())) {
        if (is.null(gDRutils::get_additional_variables(
          list(rv$response_metrics(), rv$response_data()),
          prettified = TRUE
        ))) {
          rv$triggerViz(TRUE)
          shinyjs::delay(800,
                         updateTabItems(rv$session, "main_tabs", "MetricClustering"))
        }
      } else if (isTruthy(rv$combo_object())) {
        rv$triggerViz(TRUE)
        shinyjs::delay(800,
                       updateTabItems(rv$session, "main_tabs", "DrugComboTreatment"))
      }
    })

  observeEvent(c(rv$mae_source, rv$data_src), {

    # control visibility of data source badges
    ls_in <- c(
      gDRin = "gDRin",
      gDRsearch = "gDRexplorer",
      gDRsearch2 = "gDRsearch"
    )
    if (is.null(rv$data_src)) {
      # starting application
      shinyjs::js$showSidebarTab("shiny-tab-in_plugins")
      shinyjs::show("in_plugins")
      shinyjs::delay(
        100,
        for (i in ls_in) {
          shinyjs::addClass(selector = sprintf("#%s > a > small", i), class = "hidden")
        }
      )
    } else {
      # change of data source
      for (i in ls_in) {
        shinyjs::addClass(selector = sprintf("#%s > a > small", i), class = "hidden")
      }
      tab_to_show <- ls_in[rv$data_src]
      shinyjs::removeClass(selector = sprintf("#%s > a > small", tab_to_show), class = "hidden")
    }

  }, ignoreNULL = FALSE)

  #-##### render visualisation tabs


  #-##### plugins
  gDRcomponents::init_plugins(
    plugins_config = CONFIG$plugins,
    session = session,
    rv = rv
  )

  #-##### spinner removal
  shinyjs::delay(1500, gDRcomponents::remove_spinner()) # remove it when done
  #-#####

  session$onSessionEnded(
    function() {
      logExtra::save_timestamp_mark(
        tag = S_LOG_PREFIX,
        ds = "",
        extra = data.frame(user = gDRcomponents::get_gdr_username())
      )
      logExtra::info(SE_LOG_PREFIX, sprintf("user:%s", gDRcomponents::get_gdr_username()))
    })

}
