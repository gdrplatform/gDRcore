#' merge_data
#'
#' Merge all the input data into a single dataframe
#'
#' @param manifest a data frame with a manifest info
#' @param treatments a data frame with a treaatments info
#' @param data a data frame with a raw data info
#'
#' @return a dataframe with merged data
#' @export

merge_data <- function(manifest, treatments, data) {
  # Assertions:
  stopifnot(inherits(manifest, "data.frame"))
  stopifnot(inherits(treatments, "data.frame"))
  stopifnot(inherits(data, "data.frame"))

  futile.logger::flog.info("Merging data")

  # first unify capitalization in the headers of treatments with manifest
  duplicated_col <-
    setdiff(colnames(treatments)[toupper(colnames(treatments)) %in%
                                   toupper(colnames(manifest))],
            colnames(treatments)[colnames(treatments) %in% colnames(manifest)])
  for (m_col in duplicated_col) {
    colnames(treatments)[colnames(treatments) == m_col] <-
      colnames(manifest)[toupper(m_col) == toupper(colnames(manifest))]
    futile.logger::flog.trace("Header %s in templates corrected to match case with manifest", m_col)
  }
  # merge manifest and treatment files first
  df_metadata <- merge(manifest, treatments, by = "Template")
  futile.logger::flog.info("Merging the metadata (manifest and treatment files")

  # sort out duplicate metadata columns
  duplicated_col <-
    setdiff(intersect(colnames(manifest), colnames(treatments)), "Template")
  for (m_col in duplicated_col) {
    df_metadata[, m_col] <-
      df_metadata[, paste0(m_col, ".y")] # parse template values
    missing_idx <-
      is.na(df_metadata[, m_col]) | df_metadata[, m_col] %in% c("", "-")
    # add manifest values when missing in template
    df_metadata[missing_idx, m_col] <-
      df_metadata[missing_idx, paste0(m_col, ".x")]
    # check for conflicts
    double_idx <- !(is.na(df_metadata[, paste0(m_col, ".x")]) |
                      df_metadata[, paste0(m_col, ".x")] %in% c("", "-")) &
      !(is.na(df_metadata[, paste0(m_col, ".y")]) |
          df_metadata[, paste0(m_col, ".y")] %in% c("", "-"))
    if (any(double_idx) &&
        any(df_metadata[, paste0(m_col, ".x")] != df_metadata[, paste0(m_col, ".y")], na.rm = TRUE)) {
      futile.logger::flog.warn("Merge data: metadata field %s found in both the manifest
                               and some templates with inconsistent values;
                               values in template supersede the ones in the manifest", m_col)
    }
    df_metadata[, paste0(m_col, ".x")] <- NULL
    df_metadata[, paste0(m_col, ".y")] <- NULL
  }

  # check for the expected columns
  expected_headers <- gDRutils::get_identifier("cellline")
  headersOK <- expected_headers %in% colnames(df_metadata)
  if (any(!headersOK)) {
    stop(sprintf(
      "df_metadata does not contains all expected headers: %s required",
      paste(expected_headers[!(expected_headers %in% col_df)], collpase = " ; ")
    ))
  }


  # remove wells not labeled
  df_metadata_trimmed <-
    df_metadata[!is.na(df_metadata[, gDRutils::get_identifier("drug")]),]
  futile.logger::flog.warn("%i wells discarded for lack of annotation, %i data point selected",
                           dim(df_metadata_trimmed)[1],
                           sum(is.na(df_metadata[, gDRutils::get_identifier("drug")])))

  # clean up the metadata
  cleanedup_metadata <-
    cleanup_metadata(df_metadata_trimmed)
  stopifnot(dim(cleanedup_metadata)[1] == dim(df_metadata_trimmed)[1]) # should not happen

  df_merged <- merge(cleanedup_metadata, data, by = c("Barcode",
                                                      gDRutils::get_identifier("WellPosition")))
  if (dim(df_merged)[1] != dim(data)[1]) {
    # need to identify issue and output relevant warning
    futile.logger::flog.warn("merge_data: Not all results have been matched with treatments;
                             merged table is smaller than data table")
  }
  if (dim(df_merged)[1] != dim(df_metadata)[1]) {
    # need to identify issue and print relevant warning
    futile.logger::flog.warn("merge_data: Not all treatments have been matched with results;
                             merged table is smaller than metadata table")
  }

  # remove wells not labeled
  df_raw_data <-
    df_merged[!is.na(df_merged[, gDRutils::get_identifier("drug")]), ]
  futile.logger::flog.warn("%i well loaded, %i discarded for lack of annotation, %i data point selected",
      dim(data)[1],
      sum(is.na(df_merged[, gDRutils::get_identifier("drug")])),
      dim(df_raw_data)[1]
    )

  # reorder the columns
  df_raw_data <- Order_result_df(df_raw_data)

  return(df_raw_data)
}

#' normalize_SE
#'
#' Normalize raw DR data
#'
#' @param df_raw_data a dataframe with raw data
#' @param selected_keys a vector of keys that should be included in the normalization (NULL by default)
#' @param key_values a list of values for keys that should be included in the normalization (NULL by default)
#' @param discard_keys a vector of keys that should be discarded (NULL by default)
#' @param control_mean_fct the function used for averaging controls (trimmed arithmetic mean with trim = .25 by default)
#' @param nDigits_rounding an integer with number of digits of rounding (4 by default)
#'
#' @return a SummarizedExperiment object with normalized assays
#' @export
#'
normalize_SE <- function(df_raw_data,
                  selected_keys = NULL,
                  key_values = NULL,
                  discard_keys = NULL,
                  control_mean_fct = function(x) mean(x, trim = .25), # used for averaging controls
                  nDigits_rounding = 4 # rounding of normalized response values
                ) {
    # Assertions
    stopifnot(inherits(df_raw_data, "data.frame"))
    checkmate::assert_vector(selected_keys, null.ok = TRUE)
    checkmate::assert_list(key_values, null.ok = TRUE)
    checkmate::assert_vector(discard_keys, null.ok = TRUE)
    checkmate::assert_function(control_mean_fct, null.ok = TRUE)
    checkmate::assert_number(nDigits_rounding)

    # average technical replicates and assign the right controls to each treated well
    Keys <- identify_keys(df_raw_data)
    if (!is.null(key_values)) {
      checkmate::assert_true(all(names(key_values) %in% unlist(Keys)))
    } 
    Keys$discard_keys <- discard_keys
    if (!is.null(selected_keys)) {
        Keys[names(selected_keys)] <- selected_keys[names(selected_keys)]
    }
    if (!is.null(discard_keys)) {
      Keys$DoseResp <- setdiff(Keys$DoseResp, discard_keys)
    }

    # adding 'masked = F' if missing from df_raw_data
    if ( !(gDRutils::get_identifier('masked_tag') %in% colnames(df_raw_data))) {
      df_raw_data[,gDRutils::get_identifier('masked_tag')] = FALSE
    }

    # remove background value to readout (at least 1e-10 to avoid artefactual normalized values)
    df_raw_data$CorrectedReadout = pmax(df_raw_data$ReadoutValue -
                    df_raw_data$BackgroundValue, 1e-10)
    # creates the DataFrameMatrix and fill with the treated/untreated data
    normSE <- gDR::createSE(df_raw_data, data_type = "treated", discard_keys = discard_keys)
    SummarizedExperiment::assayNames(normSE) = "Normalized"
    ctrlSE <- gDR::createSE(df_raw_data, data_type = "untreated", discard_keys = discard_keys)

    # enforced key values for end points (override selected_keys) --> for rows of the SE
    Keys$untrt_Endpoint <- setdiff(Keys$untrt_Endpoint, names(key_values))
    row_endpoint_value_filter <- array(TRUE, nrow(ctrlSE))
    if (!is.null(key_values) & length(key_values) > 0) {
        for (i in which(names(key_values) %in% names(SummarizedExperiment::rowData(ctrlSE)))) {
            if (is.numeric(key_values[i])) {
                row_endpoint_value_filter <- row_endpoint_value_filter &
                    (SummarizedExperiment::rowData(ctrlSE)[, names(key_values)[i] ] == key_values[i] &
                            !is.na(SummarizedExperiment::rowData(ctrlSE)[, names(key_values)[i] ]))
            } else {
                row_endpoint_value_filter <- row_endpoint_value_filter &
                    (SummarizedExperiment::rowData(ctrlSE)[, names(key_values)[i] ] %in% key_values[i])
            }}}

    # perform the mapping for normalization
    # first the rows
    # matching the reference end point without any treatment
    row_maps_end <- mapSE(normSE, ctrlSE, row_endpoint_value_filter, Keys)
    names(row_maps_end) <- rownames(normSE)

    # matching the reference end point with the same co-treatment (all the same but conc=0/Gnumber="vehicle")
    row_maps_cotrt <- lapply(rownames(normSE), function(x) {
      ref_metadata_idx = intersect(Keys$ref_Endpoint,
                            names(SummarizedExperiment::rowData(ctrlSE)))
      names(ref_metadata_idx) = ref_metadata_idx

      rownames(ctrlSE)[which(apply(as.matrix(  IRanges::LogicalList(
            lapply(ref_metadata_idx, function(y) # matching the metadata
              SummarizedExperiment::rowData(ctrlSE)[,y] ==
                  (SummarizedExperiment::rowData(normSE)[x, y])
            ))),
            2, all))]
        })
    names(row_maps_cotrt) <- rownames(normSE)

    # reassess the cases without a match to find equivalent drug and concentration (only 2 drugs)
    # test if one can use one of the treatment as a reference
    if ('Gnumber_2' %in% colnames(SummarizedExperiment::rowData(normSE))) {
      for (rnames in names(row_maps_cotrt)[sapply(row_maps_cotrt, length)==0]) {
        # loop through the rows without co-treatment matched
        ref_metadata_idx = setdiff(intersect(Keys$ref_Endpoint,
                              names(SummarizedExperiment::rowData(ctrlSE))),
                            c('Gnumber_2', "DrugName_2", 'Concentration_2'))
        names(ref_metadata_idx) = ref_metadata_idx

        ref_match = apply(as.matrix(  c(IRanges::LogicalList(
          lapply(ref_metadata_idx, function(y) # matching the metadata
                unlist(SummarizedExperiment::rowData(normSE)[,y, drop=F] ==
                    (SummarizedExperiment::rowData(normSE)[rnames, y, drop=F]))
              )), # matching the drugs with mapping from Gnumber to Gnumber_2
              list( Gnumber = SummarizedExperiment::rowData(normSE)$Gnumber ==
                SummarizedExperiment::rowData(normSE)[rnames,'Gnumber_2']))),
              2, all)
        if (any(ref_match)) row_maps_cotrt[rnames] = rownames(normSE)[which(ref_match)]
      }
    }

    # matching the reference at time 0 (if available)
    row_maps_T0 <- mapSE(normSE, ctrlSE, row_endpoint_value_filter, Keys, T0 = TRUE)
    names(row_maps_T0) <- rownames(normSE)

    # mapping for columns; 1 to 1 unless overridden by key_values
    col_maps <- array(colnames(ctrlSE), dimnames = list(colnames(normSE)))
    if (any(names(key_values) %in% names(SummarizedExperiment::colData(normSE)))) {
        col_maps[] <- colnames(ctrlSE)[
                which(key_values[names(key_values) %in% names(SummarizedExperiment::colData(normSE))] ==
                    SummarizedExperiment::colData(ctrlSE)[, names(SummarizedExperiment::colData(ctrlSE)) %in% names(key_values)])]
    }

    # creates the DataFrameMatrix for controls
    SummarizedExperiment::assay(normSE, "Controls") <- matrix(lapply(1:prod(dim(normSE)), function(x) S4Vectors::DataFrame()),
            nrow = nrow(normSE), ncol = ncol(normSE))

    # temporary optimization (use 'normSE_n' and 'normSE_c' to avoid using 'assay<-` in for loops)
    # TODO: refactor this part of code once we switch to DataFrameMatrix class
    ctrl_original = SummarizedExperiment::assay( aapply(ctrlSE, function(x) x[!x$masked,]) )
    # need to keep original data for the case in which reference is such that Gnumber == Gnumber_2
    normSE_n <- normSE_original <- SummarizedExperiment::assay(normSE, "Normalized")
    normSE_c <- SummarizedExperiment::assay(normSE, "Controls")


    for (i in rownames(normSE_n)) {
        for (j in colnames(normSE_n)) {
            if (nrow(normSE_original[[i, j]]) == 0) next # skip if no data
            # get all the control endpoint data
            df_end <- do.call(rbind,
                    lapply(row_maps_end[[i]], function(x) ctrl_original[[x, col_maps[j]]]))
            if (nrow(df_end) == 0) {
              futile.logger::flog.warn(
                  "Missing control data.
                  Treatment Id: %s
                  Cell_line Id: %s",
                  i,
                  j
                )
              next
            }
            df_end <- df_end[, c("CorrectedReadout",
                    intersect(Keys$untrt_Endpoint, colnames(df_end))), drop = F]
            colnames(df_end)[1] <- "UntrtReadout"
            if (ncol(df_end)>1) {
              df_end <- aggregate(df_end[, 1, drop = FALSE],
                by = as.list(df_end[, -1, drop = FALSE]),
                function(x) control_mean_fct(x))
            } else {
              df_end = DataFrame(UntrtReadout = control_mean_fct(df_end$UntrtReadout))
            }
            # reference co-treatment is not always present
            if (i %in% names(row_maps_cotrt) && length(row_maps_cotrt[[i]])>0) {
              if (all(row_maps_cotrt[[i]] %in% rownames(ctrlSE))) {
                # get all the co-treatment reference endpoint data
                df_ref <- do.call(rbind,
                        lapply(row_maps_cotrt[[i]], function(x) ctrl_original[[x, col_maps[j]]]))
                df_ref <- df_ref[, c("CorrectedReadout",
                        intersect(Keys$ref_Endpoint, colnames(df_ref))), drop = F]
                colnames(df_ref)[1] <- "RefReadout"
                if (ncol(df_ref)>1) {
                  df_ref <- aggregate(df_ref[, 1, drop = FALSE],
                    by = as.list(df_ref[, -1, drop = FALSE]),
                    function(x) control_mean_fct(x))
                } else {
                  df_ref = DataFrame(RefReadout = control_mean_fct(df_ref$RefReadout))
                }

                # check if all control have matching co-treated wells are on the same plate
                if (all(df_end$Barcode %in% df_ref$Barcode) && all(df_ref$Barcode %in% df_end$Barcode)) {
                  df_end <- merge(df_end, df_ref, # merge based on Barcode
                    by = intersect(colnames(df_end), c('Barcode', Keys$discard_keys)))
                } else {
                  futile.logger::flog.warn(
                      "Control data for the drug are propagated to other plates with co-drug controls.
                      Treatment Id: %s
                      Cell_line Id: %s",
                      i,
                      j
                    )
                  # merge (knowing that there will be NA)
                  df_end <- merge(df_end, df_ref, by = "Barcode", all = TRUE)
                  # propagate average values to the other plates
                  mean_UntrtReadout <- mean(df_end$UntrtReadout, na.rm = TRUE)
                  mean_RefReadout <- mean(df_end$RefReadout, na.rm = TRUE)
                  df_end$UntrtReadout[is.na(df_end$UntrtReadout)] <- mean_UntrtReadout
                  df_end$RefReadout[is.na(df_end$RefReadout)] <- mean_RefReadout
                }

                #gladkia: assert for control data
                if (nrow(df_end) == 0) {
                  stop(sprintf("Control dataframe failed.
                      Treatment Id: '%s'
                      Cell_line Id: %s",
                               i,
                               j))
                }
              } else if (all(row_maps_cotrt[[i]] %in% rownames(normSE))) {
                # case of the reference being with Gnumber == Gnumber_2
                ref_conc = SummarizedExperiment::rowData(normSE)[i,'Concentration_2']
                df_ref <- do.call(rbind,
                        lapply(row_maps_cotrt[[i]], function(x) {
                          if (any(normSE_original[[x, col_maps[j]]]$Concentration == ref_conc)) {
                            # the reference value with same concentration is found
                            normSE_original[[x, col_maps[j]]][
                              normSE_original[[x, col_maps[j]]]$Concentration == ref_conc,]
                          } else {
                            # the reference with proper concentration will be inferred --> need fits
                            ref_drc = normSE_original[[x, col_maps[j]]]
                            drc_fit = drc::drm(
                              CorrectedReadout ~ Concentration,
                              data = ref_drc,
                              fct = drc::LL.4(), # para = c(Hill, x_inf, x0, c50)
                              start = c(2, min(ref_drc$CorrectedReadout),
                                          max(ref_drc$CorrectedReadout),
                                          median(ref_drc$Concentration)),
                              lowerl = c(1e-5, min(ref_drc$CorrectedReadout)*.8, # wide range
                                          min(ref_drc$CorrectedReadout)*.9,
                                          min(ref_drc$Concentration)/1e3),
                              upperl =  c(12, max(ref_drc$CorrectedReadout)*1.1,
                                          max(ref_drc$CorrectedReadout)*1.2,
                                          min(ref_drc$Concentration)*1e3)
                            )
                            df_ref = data.frame(Concentration = ref_conc,
                                  CorrectedReadout = predict(drc_fit,
                                          data.frame(Concentration = ref_conc)))
                          }
                          }))
                #TODO: this piece of code is the same as above; it can be simplified after df_ref is created.
                df_ref <- df_ref[, c("CorrectedReadout",
                        intersect(Keys$ref_Endpoint, colnames(df_ref))), drop = F]
                colnames(df_ref)[1] <- "RefReadout"
                if (ncol(df_ref)>1) {
                  df_ref <- aggregate(df_ref[, 1, drop = FALSE],
                    by = as.list(df_ref[, -1, drop = FALSE]),
                    function(x) control_mean_fct(x))
                } else {
                  df_ref = DataFrame(RefReadout = control_mean_fct(df_ref$RefReadout))
                }

                # check if all control have matching co-treated wells are on the same plate
                if (all(df_end$Barcode %in% df_ref$Barcode) && all(df_ref$Barcode %in% df_end$Barcode)) {
                  df_end <- merge(df_end, df_ref,
                    by = intersect(colnames(df_end), c('Barcode', Keys$discard_keys)))
                } else {
                  futile.logger::flog.warn(
                      "Control data for the drug are propagated to other plates with co-drug controls.
                      Treatment Id: %s
                      Cell_line Id: %s",
                      i,
                      j
                    )
                  # propagate average values to the other plates
                  df_end <- merge(df_end, df_ref, by = "Barcode", all = TRUE)
                  mean_UntrtReadout <- mean(df_end$UntrtReadout, na.rm = TRUE)
                  mean_RefReadout <- mean(df_end$RefReadout, na.rm = TRUE)
                  df_end$UntrtReadout[is.na(df_end$UntrtReadout)] <- mean_UntrtReadout
                  df_end$RefReadout[is.na(df_end$RefReadout)] <- mean_RefReadout
                }
              } else {
                stop(sprintf("Reference failed.
                    Treatment Id: '%s'
                    Cell_line Id: %s",
                             i,
                             j))
              }
            } else if (i %in% names(row_maps_cotrt) && length(row_maps_cotrt[[i]])==0) {
              futile.logger::flog.warn("No reference condition found for
                  Treatment Id: '%s'
                  Cell_line Id: %s",
                           i,
                           j)
              df_end$RefReadout <- NA
            } else {
              df_end$RefReadout <- df_end$UntrtReadout
            }

            if (length(row_maps_T0[[i]]) > 0) {
              df_0 <- do.call(rbind,
                      lapply(row_maps_T0[[i]], function(x) ctrl_original[[x, col_maps[j]]]))
              df_0 <- df_0[, c("CorrectedReadout", intersect(Keys$Day0, colnames(df_0)))]
              colnames(df_0)[1] <- "Day0Readout"
              df_0 <- aggregate(df_0[, 1, drop = FALSE], by = as.list(df_0[, -1, drop = FALSE]),
                  function(x) control_mean_fct(x))

              if (!is.null(Keys$discard_keys) && all(Keys$discard_keys %in% colnames(df_0))) {
                df_ctrl <- merge(df_0[, setdiff(colnames(df_0), "Barcode")], df_end, all.y = TRUE, by = Keys$discard_keys)
              } else {
                df_ctrl <- merge(df_0[, setdiff(colnames(df_0), "Barcode")], df_end, all.y = TRUE)
                colnames(df_ctrl)[1] <- "Day0Readout"
              }
            } else {
              df_ctrl = df_end
              df_ctrl$Day0Readout = NA
            }

            # calculating the normalized response value for the control
            df_ctrl$RefRelativeViability <- round(df_ctrl$RefReadout/df_ctrl$UntrtReadout,
                nDigits_rounding)
            df_ctrl$RefGRvalue <- round(2 ** (
                    log2(df_ctrl$RefReadout / df_ctrl$Day0Readout) /
                    log2(df_ctrl$UntrtReadout / df_ctrl$Day0Readout) ), nDigits_rounding) - 1
            df_ctrl$DivisionTime <- round(
                    SummarizedExperiment::rowData(normSE)[i,gDRutils::get_identifier("duration")] /
                        log2(df_ctrl$UntrtReadout / df_ctrl$Day0Readout), nDigits_rounding)


            #gladkia: assert for merged study/control data
            ctrl_bcodes <- sort(unique(df_ctrl$Barcode))
            trt_bcodes <-
              sort(unique(normSE_original[[i, j]]$Barcode))
            # check if all treated values have matching controls on the same plate
            if (!all(trt_bcodes %in% ctrl_bcodes)) {
              # if not, propagate to all plates
              futile.logger::flog.warn(
                  "Control data are averaged and propagated to treatment plates.
                      Treatment Id: %s (plates %s)
                      Control plates: %s",
                  i, paste(trt_bcodes, collapse = ", "),
                  paste(ctrl_bcodes, collapse = ", ")
                )
              dplyr::bind_rows(df_ctrl, cbind(data.frame(Barcode = setdiff(trt_bcodes, ctrl_bcodes)),
                        t(colMeans(df_ctrl[, setdiff(colnames(df_ctrl), "Barcode")]))))
            }

            # works with by = character(0) but changes the order of rows
            df_merged <- merge(
              data.frame(normSE_original[[i, j]]),
                    data.frame(df_ctrl),
                    by = intersect(colnames(df_ctrl), c('Barcode', Keys$discard_keys)),
                    all.x = T)

            # merge the data with the controls assuring that the order of the records is preseved
            # removing this line because failed when   by = character(0)
            # df_merged <- dplyr::left_join(
            #         data.frame(normSE_original[[i, j]]),
            #         data.frame(df_ctrl),
            #         by = intersect(colnames(df_ctrl), c('Barcode', Keys$discard_keys)))

            # calculate the normalized values
            df_merged$RelativeViability <- round(df_merged$CorrectedReadout / df_merged$UntrtReadout, nDigits_rounding)

            df_merged$GRvalue = round(2 ** (
              log2(df_merged$CorrectedReadout / df_merged$Day0Readout) /
                log2(df_merged$UntrtReadout / df_merged$Day0Readout)
            ), nDigits_rounding) - 1

            # use the reference doubling Time (ReferenceDivisionTime) for GRvalue if day 0 missing
            if ( any(is.na(df_merged$Day0Readout)) ) {

                if ( !(gDRutils::get_header('add_clid')[3] %in% colnames(SummarizedExperiment::colData(normSE))) ||
                  is.na(SummarizedExperiment::colData(normSE)[j, gDRutils::get_header('add_clid')[3]]) ) {
                    # missing division time
                    futile.logger::flog.warn(paste(
                      "No day 0 information and no reference doubling time for cell line", SummarizedExperiment::colData(normSE)[j,gDRutils::get_header('add_clid')[1]],
                      "--> GR values are NA"))
                } else if (SummarizedExperiment::colData(normSE)[j, gDRutils::get_header('add_clid')[3]] >
                    1.5 * SummarizedExperiment::rowData(normSE)[i, gDRutils::get_identifier("duration")]) {
                      # check if experiment is long enough relative to division time
                      futile.logger::flog.warn(paste( "Reference doubling time for cell line",
                      SummarizedExperiment::colData(normSE)[j,gDRutils::get_header('add_clid')[1]], "is",
                      SummarizedExperiment::colData(normSE)[j, gDRutils::get_header('add_clid')[3]],
                        "which is too long for GR calculation compared to assay duration (",
                      SummarizedExperiment::rowData(normSE)[i, gDRutils::get_identifier("duration")],
                        "--> GR values are NA"))
                 } else {
                   # division time is correct for calculation
                  refDivisionTime = as.numeric(SummarizedExperiment::colData(normSE)[j, gDRutils::get_header('add_clid')[3]])

                  futile.logger::flog.warn(paste(
                    "Missing day 0 information --> calculate GR value based on reference doubling time for", SummarizedExperiment::colData(normSE)[j,gDRutils::get_header('add_clid')[1]]))

                  df_merged$GRvalue <-
                  round(2 ^ (1 + (
                    log2(pmin(1.25, # capping to avoid artefacts
                              df_merged[, "RelativeViability"])) /
                      (SummarizedExperiment::rowData(normSE)[i, gDRutils::get_identifier("duration")] / refDivisionTime)
                  )), nDigits_rounding) - 1

                  df_ctrl$RefGRvalue <-
                  round(2 ^ (1 + (
                    log2(pmin(1.25, # capping to avoid artefacts
                              df_ctrl[, "RefRelativeViability"])) /
                      (SummarizedExperiment::rowData(normSE)[i, gDRutils::get_identifier("duration")] / refDivisionTime)
                  )), nDigits_rounding) - 1
                }
            }

            # more robust assignment in case the order of df_merged has changed
            normSE_n[[i, j]] = merge(normSE_n[[i, j]],
                df_merged[, c(colnames(normSE_n[[i, j]]), 'GRvalue', 'RelativeViability')],
                by = colnames(normSE_n[[i, j]]))
            normSE_c[[i,j]] <- DataFrame(df_ctrl)
        }
    }
    metadata(normSE) <- c(metadata(normSE),
            list(df_raw_data = df_raw_data,
                Keys = Keys,
                row_maps = list(end = row_maps_end,
                                cotrt = row_maps_cotrt,
                                T0 = row_maps_T0)
                ))

    SummarizedExperiment::assay(normSE, 'Normalized') <- normSE_n
    SummarizedExperiment::assay(normSE, 'Controls') <- normSE_c

    return(normSE)
}

#' average_SE
#'
#' Avereage normalized SummarizedExperiment of DR data
#'
#' @param normSE a SummarizedExperiment with normalized DR data
#' @param TrtKeys a vector of keys used for averaging (NULL by default)
#'
#' @return a SummarizedExperiment with additional assay with averaged DR data
#' @export
#'

average_SE <- function(normSE, TrtKeys = NULL) {

  # Assertions:
  checkmate::assert_class(normSE, "SummarizedExperiment")
  checkmate::assert_vector(TrtKeys, null.ok = TRUE)


    avgSE <- normSE
    if (is.null(TrtKeys)) {
        if ("Keys" %in% names(metadata(normSE))) {
          TrtKeys <- metadata(normSE)$Keys$Trt
          TrtKeys <- setdiff(TrtKeys, metadata(normSE)$Keys$discard_keys)
        } else {
          TrtKeys <- identify_keys(normSE)$Trt
        }
    }
    metadata(normSE)$Keys$Trt <- TrtKeys

    SummarizedExperiment::assay(avgSE, "Averaged") <- SummarizedExperiment::assay(avgSE, "Normalized")
    avgSE <- aapply(avgSE, function(x) {
        if (nrow(x) > 1) {
            subKeys <- intersect(TrtKeys, colnames(x))
            if (all(x$masked)) {
              df_ = as.data.frame(matrix(0,0,length(subKeys)+5))
              colnames(df_) = c(subKeys,
                    c("GRvalue", "RelativeViability","CorrectedReadout"),
                    paste0("std_", c("GRvalue", "RelativeViability")))
              return(df_)
            }
            df_av <- aggregate(x[ !x$masked ,
                                  c("GRvalue", "RelativeViability","CorrectedReadout")],
                            by = as.list(x[ !x$masked , subKeys, drop = FALSE]),
                            FUN = function(y) mean(y, na.rm = TRUE))
            df_std <- aggregate(x[!x$masked, c("GRvalue", "RelativeViability")],
                                by = as.list(x[ !x$masked, subKeys, drop = FALSE]),
                                FUN = function(x) sd(x, na.rm = TRUE))
            colnames(df_std)[colnames(df_std) %in% c("GRvalue", "RelativeViability")] =
                paste0("std_",
                    colnames(df_std)[colnames(df_std) %in% c("GRvalue", "RelativeViability")])
            return( merge(df_av, df_std, by = subKeys) )
        } else return(x)
    }, "Averaged")

    SummarizedExperiment::assay(avgSE, "Avg_Controls") <- SummarizedExperiment::assay(avgSE, "Controls")
    avgSE <- aapply(avgSE, function(x) {
        if (nrow(x) > 1) {
            subKeys <- intersect(TrtKeys, colnames(x))
            df_av <- DataFrame(lapply(x[, c("Day0Readout", "UntrtReadout",
                    "RefGRvalue", "RefRelativeViability",
                    "RefReadout", "DivisionTime")], FUN = function(y) mean(y, na.rm = TRUE)))
            return( df_av )
        } else return(x)
    }, "Avg_Controls")

    return(avgSE)
}

#' metrics_SE
#'
#' Calculate metrics for DR data
#'
#' @param avgSE a SummarizedExperiment with averaged and normalized assays
#' @param studyConcThresh a numeric with study concentration threshold (4 by default)
#'
#' @return a SummarizedExperiment with additional assay with metrics
#' @export
#'

metrics_SE = function(avgSE, studyConcThresh = 4) {

    # Assertions:
    checkmate::assert_class(avgSE, "SummarizedExperiment")
    checkmate::assert_number(studyConcThresh)

    stopifnot(is.numeric(studyConcThresh))
    # this is not used as we enforce the same conditions as the input SE; not collapsing allowed
    # if (is.null(DoseRespKeys)) {
    #     if ("Keys" %in% names(metadata(avgSE))) DoseResp = metadata(avgSE)$Keys$DoseResp
    #     else DoseRespKeys = identify_keys(avgSE)$DoseResp
    # } else {
    #     metadata(avgSE)$Keys$DoseResp = DoseRespKeys
    # }

    metricsSE <- avgSE
    SummarizedExperiment::assay(metricsSE, "Metrics") <- SummarizedExperiment::assay(metricsSE, "Averaged")

    # temporary optimization (use 'normSE_n' and 'normSE_c' to avoid using 'assay<-` in a foor loops)
    # TODO: refactor this part of code once we switch to DataFrameMatrix class
    mSE_m <- SummarizedExperiment::assay(metricsSE, "Metrics")
    a_SE = SummarizedExperiment::assay(metricsSE, "Averaged")
    aCtrl_SE = SummarizedExperiment::assay(metricsSE, "Avg_Controls")
    for (i in rownames(metricsSE)) {
        for (j in colnames(metricsSE)) {
            df_ <- a_SE[[i, j]]
            if (!is.null(df_) && length(df_) > 0 && nrow(df_) > 0) { # studyConcThresh is embeded in RVGRfits
                mSE_m[[i, j]] <- DataFrame(gDRutils::RVGRfits(df_,
                    e_0 = aCtrl_SE[[i, j]]$RefRelativeViability,
                    GR_0 = aCtrl_SE[[i, j]]$RefGRvalue,
                    n_point_cutoff = studyConcThresh))
            } else {
                out <- DataFrame(matrix(NA, 0, length(gDRutils::get_header("response_metrics"))+2))
                colnames(out) <- c(gDRutils::get_header("response_metrics"), "maxlog10Concentration", "N_conc")
                mSE_m[[i, j]] <- out
            }
        }
    }
    SummarizedExperiment::assay(metricsSE, "Metrics") <- mSE_m
    return(metricsSE)
}

#' identify_keys
#'
#' Identify keys in the DR data represented by dataframe or SummarizedExperiment or MultiAssayExperiment objects
#'
#' @param df_se_mae a dataframe or SummarizedExperiment or MultiassayExperiment with keys
#'
#' @return a list of keys
#' @export
#'

identify_keys <- function(df_se_mae) {

  # Assertions:
  stopifnot(inherits(df_se_mae, c("data.frame", "MultiAssayExperiment", "SummarizedExperiment")))


    if (any(class(df_se_mae) %in% c("MultiAssayExperiment", "SummarizedExperiment"))) {
        if ("MultiAssayExperiment" %in% class(df_se_mae)) {
            # if MAE, convert to SE based on the treated SE (could be optimized)
            df_se_mae <- df_se_mae[["treated"]]
            se_untrt <-  df_se_mae[["untreated"]]
        } else se_untrt <- NULL
        all_keys <- unique(c(
            colnames(SummarizedExperiment::rowData(df_se_mae)),
            colnames(SummarizedExperiment::colData(df_se_mae)),
            unlist(lapply(SummarizedExperiment::assay(df_se_mae), colnames))))
    } else { # case of a data frame
        all_keys <- colnames(df_se_mae)
    }

    keys <- list(Trt = setdiff(all_keys, "Barcode"),
            DoseResp = setdiff(all_keys,  "Barcode"),
            ref_Endpoint = setdiff(all_keys, c("Concentration",
                                            gDRutils::get_identifier("drug"),
                                            gDRutils::get_identifier("drugname"))),
            untrt_Endpoint = all_keys[ c(-agrep("Concentration", all_keys),
                                            -agrep(gDRutils::get_identifier("drug"), all_keys),
                                            -agrep(gDRutils::get_identifier("drugname"), all_keys))])
    keys[["Day0"]] <- setdiff(keys[["untrt_Endpoint"]], gDRutils::get_identifier("duration"))
    keys <- lapply(keys, function(x) setdiff(x, c(gDRutils::get_header("raw_data"),
        gDRutils::get_header("normalized_results"), "Template", gDRutils::get_identifier("WellPosition"), gDRutils::get_header("averaged_results"),
            gDRutils::get_header("metrics_results"), "ReferenceDivisionTime"
    )))
    keys <- lapply(keys, sort)

    # check if all values of a key is NA
    for (k in keys[["untrt_Endpoint"]]) {

        if ("SummarizedExperiment" %in% class(df_se_mae)) {
            # check the metadata fields for NA
            if (k %in% colnames(SummarizedExperiment::rowData(df_se_mae))) df_ <- SummarizedExperiment::rowData(df_se_mae)
            else if (k %in% colnames(SummarizedExperiment::colData(df_se_mae))) df_ <- SummarizedExperiment::colData(df_se_mae)
            else next # not a metadata

            if (all(is.na(df_[,k]))) keys <- lapply(keys, function(x) setdiff(x, k))

            if (!is.null(se_untrt) && k %in% colnames(SummarizedExperiment::rowData(se_untrt))) {
                df_ <- SummarizedExperiment::rowData(se_untrt)
                if (all(is.na(df_[df_[,gDRutils::get_identifier("duration")] == 0, k]))) {
                    keys[["Day0"]] <- setdiff(keys[["Day0"]], k)
                }
            }
        } else { # case of a data frame
            if (all(is.na(df_se_mae[, k]))) {
                keys <- lapply(keys, function(x) setdiff(x, k))
            }
            if (all(is.na(df_se_mae[df_se_mae[,gDRutils::get_identifier("duration")] == 0, k]))) {
                keys[["Day0"]] <- setdiff(keys[["Day0"]], k)
            }
        }
    }
  return(keys)
}

#' cleanup_metadata
#'
#' Cleanup a dataframe with metadata
#'
#' @param df_metadata a dataframe with metadata
#'
#' @return a dataframe with cleaned metadata
#'
#' @export

cleanup_metadata <- function(df_metadata) {

  # Assertions:
  stopifnot(inherits(df_metadata, "data.frame"))

  # clean up numberic fields
  df_metadata[, gDRutils::get_identifier("duration")] <-
    round(as.numeric(df_metadata[, gDRutils::get_identifier("duration")]), 6)
  # identify potential numeric fields and replace NA by 0 - convert strings in factors
  for (c in setdiff(1:dim(df_metadata)[2], c(
    agrep(gDRutils::get_identifier("drug"), colnames(df_metadata)),
    agrep("Concentration", colnames(df_metadata)),
    grep(paste(
      c(
        gDRutils::get_identifier("cellline"),
        gDRutils::get_header("manifest"),
        gDRutils::get_identifier("WellPosition")
      ),
      collapse = "|"
    ), colnames(df_metadata))
  ))) {
    vals <- unique(df_metadata[, c])

    if (is.character(vals)) {
      num_vals <- as.numeric(vals)
      if (sum(is.na(num_vals)) > 2 || all(is.na(num_vals))) {
        df_metadata[, c] <- factor(df_metadata[, c])
        futile.logger::flog.warn("Metadata field %s converted to factors",
                colnames(df_metadata)[c])
      } else {
        is.na(df_metadata[, c]) <- 0
        df_metadata[, c] <- as.numeric(df_metadata[, c])
        futile.logger::flog.warn("Metadata field %s converted to numeric values",
                colnames(df_metadata)[c])
      }
    }
  }
    # TODO: specific to GNE database --> need to be replaced by a function
    df_metadata <- add_CellLine_annotation(df_metadata)

    # check that Gnumber_* are in the format 'G####' and add common name (or Vehicle or Untreated)

    for (i in agrep(gDRutils::get_identifier("drug"), colnames(df_metadata))) { # correct case issues
        for (w in gDRutils::get_identifier("untreated_tag")) {
            df_metadata[grep(w, df_metadata[,i], ignore.case = T),i] <- w
        }
    }
    # -----------------------

    df_metadata <- add_Drug_annotation(df_metadata)

    # clean up concentration fields
    for (i in agrep("Concentration", colnames(df_metadata))) {
        trt_n <- ifelse(regexpr("_\\d", colnames(df_metadata)[i]) > 0,
                            substr(colnames(df_metadata)[i], 15, 20), 1)
        DrugID_col <- ifelse(trt_n == 1, gDRutils::get_identifier("drug"), paste0(gDRutils::get_identifier("drug"), "_", trt_n))
        df_metadata[df_metadata[,DrugID_col] %in% gDRutils::get_identifier("untreated_tag"), i] <- 0 # set all untreated to 0

        DrugID_0 <- setdiff(unique(df_metadata[ df_metadata[,i] == 0, DrugID_col]), gDRutils::get_identifier("untreated_tag"))
        DrugID_0 <- DrugID_0[!is.na(DrugID_0)]
        if (length(DrugID_0) > 0) {
          futile.logger::flog.warn("Some concentration for %s are 0: %s",
                                   DrugID_col,
                                   paste(DrugID_0, collapse = " ; "))

        }
        df_metadata[,i] <- round(as.numeric(df_metadata[, i]), 6) # avoid mismatch due to string truncation
    }
    df_metadata[, i] <-
      round(as.numeric(df_metadata[, i]), 6) # avoid mismatch due to string truncation

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

Order_result_df <- function (df_) {

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
        gDRutils::get_header("add_clid")[1],
        gDRutils::get_identifier("duration"),
        gDRutils::get_identifier("drugname"),
        "Concentration",
        paste0(c(
          paste0(gDRutils::get_identifier("drugname"), "_"), "Concentration_"
        ),
        sort(c(2:10, 2:10))),
        setdiff(colnames(df_), c(
          gDRutils::get_header("ordered_1"), gDRutils::get_header("ordered_2")
        ))
      ),
      cols
    )

  df_ <- df_[do.call(order, df_[, row_order_col]), cols]

  return(df_)
}


#' add_CellLine_annotation
#'
#' add cellline annotation to a dataframe with metadata
#'
#' @param df_metadata a dataframe with metadata
#' @param fill_DB_wiith_unknown a logical indicating whether DB should be filled with unknown cell lines
#' 
#' @return a dataframe with metadata with annotated cell lines
#' @export

add_CellLine_annotation <- function(df_metadata,
                                    fill_DB_wiith_unknown = FALSE) {
  
    # Assertions:
    stopifnot(inherits(df_metadata, "data.frame"))
    checkmate::assert_logical(fill_DB_wiith_unknown)
    
    DB_cellid_header <- "cell_line_identifier"
    DB_cell_annotate <- c("cell_line_name", "primary_tissue", "doubling_time")
    # corresponds to columns gDRutils::get_header("add_clid"): name, tissue, doubling time
    
    # the logic of adding celline annotation for df_metadata is based on the function get_cell_lines from the gDRwrapper
    # we added additional parameter 'fill_DB_wiith_unknown' that allows to fill the DB with clid info for these cell lines
    # that are not present in the DB.
    # Other fields are set as "UNKNOWN". If the fill_DB_wiith_unknown is set as FALSE we add unkonown cell lines
    # only to the tibble.
    # This approach will be corrected once we will implement final solution for adding cell lines.
    validateCLs <- gDRwrapper::validate_cell_lines(unique(df_metadata[,gDRutils::get_identifier("cellline")]))
    if(!validateCLs){
      missingTblCellLines <- tibble::tibble(parental_identifier = "UNKNOWN",
                                            cell_line_name = "UNKNOWN",
                                            cell_line_identifier = unique(df_metadata[,gDRutils::get_identifier("cellline")]),
                                            doubling_time = "UNKNOWN",
                                            primary_tissue = "UNKNOWN",
                                            subtype = "UNKNOWN")
      
      if(fill_DB_wiith_unknown){
        addMissingCellLines <- gDRwrapper::add_drugs(missingTblCellLines)
      }
    }
    CLs_info <- tryCatch( {
        CLs_info <- gDRwrapper::get_cell_lines()
        CLs_info <- CLs_info[CLs_info$cell_line_identifier %in% unique(df_metadata[,gDRutils::get_identifier("cellline")]),]
        CLs_info <- CLs_info[,c(DB_cellid_header,DB_cell_annotate)]
        CLs_info
    }, error = function(e) {
      futile.logger::flog.error("Failed to load cell line info from DB: %s", e)
        data.frame()
    })

    if (nrow(CLs_info) == 0) return(df_metadata)

    colnames(CLs_info) <- c(gDRutils::get_identifier("cellline"), gDRutils::get_header("add_clid"))
    CLIDs <- unique(df_metadata[,gDRutils::get_identifier("cellline")])
    bad_CL <- !(CLIDs %in% (CLs_info %>% dplyr::pull(gDRutils::get_identifier("cellline"))))
    if (any(bad_CL)) {
        futile.logger::flog.warn("Cell line ID %s not found in cell line database",
                     paste(CLIDs[bad_CL], collapse = " ; "))
        temp_CLIDs = data.frame(CLIDs[bad_CL], CLIDs[bad_CL])
        temp_CLIDs[, 1+(2:length(gDRutils::get_header("add_clid")))] = NA
        colnames(temp_CLIDs) = c(gDRutils::get_identifier("cellline"),
                      gDRutils::get_header("add_clid"))
        CLs_info = rbind(CLs_info, temp_CLIDs)
        }

    futile.logger::flog.info("Merge with Cell line info")
    nrows_df <- nrow(df_metadata)
    df_metadata <- base::merge(df_metadata, CLs_info, by = gDRutils::get_identifier("cellline"), all.x = TRUE)
    stopifnot(nrows_df == nrow(df_metadata))

    return(df_metadata)

}

#' add_Drug_annotation
#'
#' add drug annotation to a dataframe with metadata
#'
#' @param df_metadata a dataframe with metadata
#' @param fill_DB_wiith_unknown a logical indicating whether DB should be filled with unknown drugs
#' 
#'
#' @return a dataframe with metadata with annotated drugs
#' @export


add_Drug_annotation <- function(df_metadata,
                                fill_DB_wiith_unknown = FALSE) {
  
        # Assertions:
        stopifnot(inherits(df_metadata, "data.frame"))
        checkmate::assert_logical(fill_DB_wiith_unknown)
  
        nrows_df <- nrow(df_metadata)

        DB_drug_identifier <- "gnumber"
        # the logic of adding drug annotation for df_metadata is based on the function get_drugs from the gDRwrapper
        # we added additional parameter 'fill_DB_wiith_unknown' that allows to fill the DB with drug_name and gnumber, for these drugs,
        # that are not present in the DB
        # Other fields are set as "UNKNOWN". If the fill_DB_wiith_unknown is set as FALSE we add unkonown cell lines
        # only to the tibble.
        # This approach will be corrected once we will implement final solution for adding cell lines.

        drugsTreated <- unique(df_metadata[, gDRutils::get_identifier("drug")])
        
        drugsTreated <- drugsTreated[!drugsTreated%in% gDRutils::get_identifier("untreated_tag")]
        validateDrugs <- gDRwrapper::validate_drugs(drugsTreated)
        if(!validateDrugs){
          missingTblDrugs <- tibble::tibble(drug_name = drugsTreated,
                                            drug_moa = "UNKNOWN",
                                            gnumber = drugsTreated)
          if(fill_DB_wiith_unknown){
            addMissingDrugs <- gDRwrapper::add_drugs(missingTblDrugs)
          }
          
        }
        Drug_info <- tryCatch({
          # TODO: refactor this part of code once we switch to DataFrameMatrix class
          gDrugs <- gDRwrapper::get_drugs()[, c(DB_drug_identifier, "drug_name")]
          #gDrugs[, 1] <- gsub("\\..*", "", gDrugs$gnumber) # remove batch number from DB_drug_identifier
          gDrugs
        }, error = function(e) {
          futile.logger::flog.error("Failed to load drug info from DB: %s", e)
            data.frame()
        })

        if (nrow(Drug_info) == 0) {
            df_metadata[, gDRutils::get_identifier("drugname")] = df_metadata[, gDRutils::get_identifier("drug")]
            return(df_metadata)
        }
        
        # -----------------------

        colnames(Drug_info) <- c("drug", "DrugName")
        Drug_info <-
          rbind(data.frame(
            drug = gDRutils::get_identifier("untreated_tag"),
            DrugName = gDRutils::get_identifier("untreated_tag")
          ),
          Drug_info)
        Drug_info <- dplyr::distinct(Drug_info, drug, .keep_all = TRUE)
        DrIDs <- unique(unlist(df_metadata[,agrep(gDRutils::get_identifier("drug"), colnames(df_metadata))]))
        if(any(!drugsTreated %in% Drug_info$drug)){
          Drug_info <- rbind(Drug_info, data.table::setnames(missingTblDrugs[!drugsTreated %in% Drug_info$drug, c(3,1)], names(Drug_info)))
        }
        bad_DrID <- !(DrIDs %in% Drug_info$drug) & !is.na(DrIDs)
        if (any(bad_DrID)) {
            # G number, but not registered
            ok_DrID <- attr(regexpr("^G\\d*",DrIDs), "match.length")==9
            if (any(ok_DrID)) {
              futile.logger::flog.warn("cleanup_metadata: Drug %s  not found in gCSI database; use G# as DrugName",
                                       paste(DrIDs[ok_DrID & bad_DrID], collapse = " ; "))
              Drug_info <-
                rbind(Drug_info, data.frame(drug = DrIDs[ok_DrID & bad_DrID],
                                            DrugName = DrIDs[ok_DrID & bad_DrID]))
            } else {
              futile.logger::flog.error("Drug %s not in the correct format for database",
                  paste(DrIDs[!ok_DrID], collapse = ' ; '))
            }
        }
        colnames(Drug_info)[2] <- gDRutils::get_identifier("drugname")
        futile.logger::flog.info("Merge with Drug_info for Drug 1")
        df_metadata <- base::merge(df_metadata, Drug_info, by.x = gDRutils::get_identifier("drug"), by.y = "drug", all.x = TRUE)
        # add info for columns Gnumber_*
        for (i in grep(paste0(gDRutils::get_identifier("drug"),"_\\d"), colnames(df_metadata))) {
            df_metadata[ is.na(df_metadata[,i]), i] = gDRutils::get_identifier("untreated_tag")[1] # set missing values to Untreated
            Drug_info_ <- Drug_info
            colnames(Drug_info_)[2] <- paste0(colnames(Drug_info_)[2], substr(colnames(df_metadata)[i], 8, 12))
            futile.logger::flog.info("Merge with Drug_info for %s", i)
            df_metadata <- merge(df_metadata, Drug_info_, by.x = i, by.y = "drug", all.x = TRUE)
        }
        df_metadata[, colnames(df_metadata)[grepl(gDRutils::get_identifier("drugname"), colnames(df_metadata))]] =
          droplevels(df_metadata[, colnames(df_metadata)[grepl(gDRutils::get_identifier("drugname"), colnames(df_metadata))]])

    stopifnot(nrows_df == nrow(df_metadata))

    return(df_metadata)

}


#' mapSE
#'
#' Perfmorm mapping for normalization
#'
#' @param normSE a SummarizedExperiment with normalization assaay
#' @param ctrlSE a SummarizedExperiment object with information for controls
#' @param row_endpoint_value_filter an array with key values for end points
#' @param Keys a list of all identified keys
#' @param T0 a logical indicating if the mapping should be performer for Time=0 (FALSE by default)
#'
#' @return a list of mapping
#' @export

mapSE <- function(normSE, ctrlSE, row_endpoint_value_filter, Keys, T0 = FALSE){
    # Assertions:
    checkmate::assert_class(normSE, "SummarizedExperiment")
    checkmate::assert_class(ctrlSE, "SummarizedExperiment")
    checkmate::assert_array(row_endpoint_value_filter)
    checkmate::assert_list(Keys)
    checkmate::assert_logical(T0)

    mappingFactor <- ifelse(T0, "Day0", "untrt_Endpoint")

    keyValuesList <- list(key_values = row_endpoint_value_filter)
    if(T0){
      keyValuesList <- list(T0 = SummarizedExperiment::rowData(ctrlSE)[, gDRutils::get_identifier("duration")] == 0)
    }

    matchFactor <- ifelse(T0, "T0", gDRutils::get_identifier("duration"))

    lapply(rownames(normSE), function(x) {
    # define matix with matching metadata
      ctrl_metadata_idx = intersect(Keys[[mappingFactor]],
                                    names(SummarizedExperiment::rowData(ctrlSE)))
      names(ctrl_metadata_idx) = ctrl_metadata_idx
      match_mx <-
        IRanges::LogicalList(c(
          lapply(ctrl_metadata_idx, function(y) # matching the metadata
            SummarizedExperiment::rowData(ctrlSE)[,y] ==
              SummarizedExperiment::rowData(normSE)[x,y]),
          c(keyValuesList, list(conc = apply(cbind(array(0, nrow(ctrlSE)), # padding to avoid empty df;
                                  SummarizedExperiment::rowData(ctrlSE)[, agrep("Concentration",
                                                                                colnames(SummarizedExperiment::rowData(ctrlSE))), drop = FALSE]), 1,
                            function(x)
                              all(x == 0))

          ))))
      match_idx <- which(apply(as.matrix(match_mx), 2, all)) # test matching conditions
      if (length(match_idx) == 0) {
        # if not exact match, try to find best match (as many metadata fields as possible)
        futile.logger::flog.warn("Missing untreated controls %s for: %s",
            ifelse(T0, '(T=0)', '(endpoint)'), x)
        idx <-
          apply(as.matrix(match_mx), 2, function(y)
            sum(y, na.rm = TRUE)) *
          match_mx[[matchFactor]]
        if (any(idx > 0)) {
          match_idx <- which.max(idx)
          futile.logger::flog.warn("Found partial match:",
                                   rownames(ctrlSE)[match_idx])
        } else { # failed to find any potential match
          futile.logger::flog.warn("No partial match found")
        }
      }
      return(rownames(ctrlSE)[match_idx])
    })
}



#' Add codrug group
#'
#' @param SE 
#'
#' @return
#' @export
#'
#' @examples
add_codrug_group = function(SE) {

  r_data = SummarizedExperiment::rowData(SE)
  if (!('Gnumber_2' %in% colnames(r_data))) return(SE)

  # find the pairs of drugs with relevant metadata
  drug_ids = paste0(gDRutils::get_identifier()$drug, c('', '_2'))
  other_metadata = setdiff(colnames(r_data), c('Concentration_2', drug_ids,
                    paste0(gDRutils::get_identifier()$drugname, c('', paste0('_',1:10)))))
  drug_pairs = unique(r_data[, c(drug_ids, other_metadata)])
  drug_pairs = drug_pairs[ !(drug_pairs[,drug_ids[2]] %in% gDRutils::get_identifier('untreated_tag')),]

  pair_list = vector('list', nrow(drug_pairs))
  # loop through the pairs to assess the number of individual concentration pairs
  for (idp in 1:nrow(drug_pairs)) {
    row_idx = r_data[,drug_ids[1]] %in% unlist(drug_pairs[idp, drug_ids]) &
            r_data[,drug_ids[2]] %in% c(unlist(drug_pairs[idp, drug_ids]),
                gDRutils::get_identifier('untreated_tag')) &
            apply(as.matrix(
                IRanges::LogicalList(c(
                  lapply(other_metadata, function(y) # matching the metadata
                    r_data[,y] == drug_pairs[idp,y])
                  ))), 2, all)

    # reverse engineer the type of combination experiment
    flat_data = gDRutils::assay_to_df(SE[row_idx, ], 'Averaged')
    flat_data = flat_data[flat_data$Concentration_2 > 0,]
    conc_1 = table(flat_data$Concentration)
    conc_2 = table(flat_data$Concentration_2)
    n_conc_pairs = nrow(unique(flat_data[,c('Concentration', 'Concentration_2')]))
    conc_ratio = table(round(log10(flat_data$Concentration / flat_data$Concentration_2),2))
    conc_ratio = conc_ratio[!names(conc_ratio) %in% c('Inf', '-Inf')]

    condition = paste(paste(other_metadata, unlist(drug_pairs[idp,other_metadata]), sep = '='),
                  collapse=' ')
    if (length(conc_ratio) <= 2) {
      type = 'co-dilution'
      print(sprintf('Found %s combination with %s and %s: ratio of %.2f, %i concentrations (%s)',
          type, drug_pairs[idp,1], drug_pairs[idp,2], 10**as.numeric(names(conc_ratio)),
            length(conc_2), condition))
    } else if (n_conc_pairs == length(conc_1)*length(conc_2) & length(conc_2) >= 4) {
      type = 'matrix'
      print(sprintf('Found %s combination with %s and %s: %i x %i concentrations (%s)',
          type, drug_pairs[idp,1], drug_pairs[idp,2], length(conc_1), length(conc_2), condition))
    } else if (length(conc_2)<4) {
      type = 'fixed'
      print(sprintf('Found %s combination of %s with %s at %.3g uM (%s)',
          type, drug_pairs[idp,1], drug_pairs[idp,2], as.numeric(names(conc_2)), condition))
    } else {
      type = 'other'
      print(sprintf('Found %s combination with %s and %s: %i concentration pairs (%s)',
        type, drug_pairs[idp,1], drug_pairs[idp,2], n_conc_pairs, condition))
    }

    pair_list[[idp]] = list(Gnumbers = unlist(drug_pairs[idp,]),
                          rows = rownames(r_data)[row_idx],
                          type = type)
  }

  metadata(SE)$drug_combinations = pair_list
  return(SE)
}

      # for (iCL in 1:ncol(SE)) {
      #     flat_data = assay_to_df(SE[row_idx, iCL], 'Averaged')
      #     tail(unique(flat_data[, c('Gnumber', 'Concentration', 'Gnumber_2', 'Concentration_2')]),30)
      # }
