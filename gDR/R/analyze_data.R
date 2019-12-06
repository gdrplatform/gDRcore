
#' @import gneDB
#' @import reshape2
#' @import dplyr

#########################################
### TODO:
#########################################


#' @export
Overall_function = function(manifest_file, template_file, results_file,
                output_files, selected_keys = NULL, key_values = NULL,
                instrument='EnVision') {
    # output_files should contain file names for :
    #   log_file, QC_file, raw_result, process_results, metrics_results

    log_str = 'Report from gDR pipeline'

    lData = load_data(manifest_file, template_file, results_file,
                log_str, instrument)
    df_raw_data = merge_data(lData$manifest, lData$treatments, lData$data, log_str)

    #returning two SEs with 'df_raw_data' assay and creating MAE from it
    untreatedSE <-
      gDR::createSE(df_raw_data, data_type = "untreated")
    treatedSE <-
      gDR::createSE(df_raw_data, data_type = "treated")
    rawMAE <- gDR::createMAE_raw(untreatedSE, treatedSE)

    # output_QC_byPlate(df_raw_data, output_files['QC_file']) # TODO: check column/row bias

    # TODO: this function could be overload to have an MAE as input; ok for now
    Keys = identify_keys(df_raw_data) # may be manually changed
    if (!is.null(selected_keys)) {
        Keys[names(selected_keys)] = selected_keys[names(selected_keys)]
    }

    # # adding df_normalized assay
    # df_normalized_assay <- get_normalized_data_asay(mae, log_str, Keys, key_values)
    # gDR::addAssayToMAE(mae,
    #                    assay = df_normalized_assay,
    #                    assay_name = "df_normalized",
    #                    exp_name = "treated")
    # #
    # adding averaged
    # df_averaged_assay <- get_averaged_data_asay(treatedSE, untreatedSE, Keys$Trt)
    # gDR::addAssayToMAE(mae,
    #                    assay = df_averaged_assay,
    #                    assay_name = "df_averaged",
    #                    exp_name = "treated")

    # adding metrics
    # df_metrics_assay <- get_metrics_data_asay(treatedSE, untreatedSE, Keys$DoseResp)
    # gDR::addAssayToMAE(mae,
    #                    assay = df_metrics_assay,
    #                    assay_name = "df_metrics",
    #                    exp_name = "treated")
    # return(mae)

    df_normalized = normalize_data(df_raw_data, log_str, Keys, key_values)
    df_averaged = average_replicates(df_normalized, Keys$Trt)

    df_metrics = calculate_DRmetrics(df_averaged, Keys$DoseResp)

    log_file <- file(output_files['log_file'], open = "wt")
    writeLines(log_str, log_file)
    close(log_file)


    return(list(raw=df_raw_data,
            normalized=df_normalized,
            averaged=df_averaged,
            metrics=df_metrics))
}


#' @export
merge_data = function(manifest, treatments, data, log_str) {

    log_str = c(log_str, '', 'merge_data')

    # first unify capitalization in the headers of treatments with manifest
    duplicated_col = setdiff(colnames(treatments)[ toupper(colnames(treatments)) %in%
                                                    toupper(colnames(manifest)) ],
                            colnames(treatments)[ colnames(treatments) %in% colnames(manifest) ])
    for (m_col in duplicated_col) {
        colnames(treatments)[ colnames(treatments) == m_col] =
                colnames(manifest)[ toupper(m_col) == toupper(colnames(manifest))]
        print(paste('Header', m_col, "in templates corrected to match case with manifest"))
    }
    # merge manifest and treatment files first
    df_metadata = merge(manifest, treatments, by = 'Template')
    print('Merging the metadata files:')
    print(head(df_metadata))

    # sort out duplicate metadata columns
    duplicated_col = setdiff(intersect(colnames(manifest), colnames(treatments)), 'Template')
    for (m_col in duplicated_col) {
        df_metadata[,m_col] = df_metadata[,paste0(m_col,'.y')] # parse template values
        missing_idx = is.na(df_metadata[,m_col]) | df_metadata[,m_col] %in% c('', '-')
        # add manifest values when missing in template
        df_metadata[missing_idx,m_col] = df_metadata[missing_idx,paste0(m_col,'.x')]
        # check for conflicts
        double_idx = !(is.na(df_metadata[,paste0(m_col,'.x')]) |
                            df_metadata[,paste0(m_col,'.x')] %in% c('', '-')) &
                     !(is.na(df_metadata[,paste0(m_col,'.y')]) |
                            df_metadata[,paste0(m_col,'.y')] %in% c('', '-'))
        if (any(double_idx) &&
            any(df_metadata[,paste0(m_col,'.x')] != df_metadata[,paste0(m_col,'.y')], na.rm=T)) {
                WarnMsg = paste('Metadata field', m_col,
                    'found in both the manifest and some templates with inconsistent values;',
                    'values in template supersede the ones in the manifest')
                log_str = c(log_str, 'Warning in merge_data:')
                log_str = c(log_str, WarnMsg)
                warning(WarnMsg)
            }
        df_metadata[,paste0(m_col,'.x')] = NULL
        df_metadata[,paste0(m_col,'.y')] = NULL
    }

    # check for the expected columns
    expected_headers = get_identifier('cellline')
    headersOK = expected_headers %in% colnames(df_metadata)
    if (any(!headersOK)) {
        ErrorMsg = paste('df_metadata',
            'does not contains all expected headers: ',
            paste(expected_headers[ !(expected_headers %in% col_df) ], collpase = ' ; '),
            ' required')
        log_str = c(log_str, 'Error in merge_data:')
        log_str = c(log_str, ErrorMsg)
        stop(ErrorMsg)
    }


    # remove wells not labeled
    df_metadata_trimmed = df_metadata[ !is.na(df_metadata[,get_identifier('drug')]), ]
    WarnMsg = sprintf('%i wells discarded for lack of annotation, %i data point selected',
                dim(df_metadata_trimmed)[1],
                sum(is.na(df_metadata[,get_identifier('drug')])))

    # clean up the metadata
    print(colnames(df_metadata_trimmed))
    print(dim(df_metadata_trimmed))
    cleanedup_metadata = cleanup_metadata(df_metadata_trimmed, log_str)
    print(colnames(cleanedup_metadata))
    print(dim(cleanedup_metadata))
    stopifnot( dim(cleanedup_metadata)[1] == dim(df_metadata_trimmed)[1] ) # should not happen

    df_merged = merge(cleanedup_metadata, data, by = c('Barcode',
                                get_identifier('WellPosition')))
    if (dim(df_merged)[1] != dim(data)[1]) {# need to identify issue and output relevant warning
        WarnMsg = 'Not all results have been matched with treatments; merged table is smaller than data table'
        log_str = c(log_str, 'Warning in load_merge_data:')
        log_str = c(log_str, WarnMsg)
        warning(WarnMsg)
    }
    if (dim(df_merged)[1]!=dim(df_metadata)[1]) {#need to identify issue and print relevant warning
        WarnMsg = 'Not all treatments have been matched with results; merged table is smaller than metadata table'
        log_str = c(log_str, 'Warning in load_merge_data:')
        log_str = c(log_str, WarnMsg)
        warning(WarnMsg)
    }

    # remove wells not labeled
    df_raw_data = df_merged[ !is.na(df_merged[,get_identifier('drug')]), ]
    WarnMsg = sprintf('%i well loaded, %i discarded for lack of annotation, %i data point selected',
                dim(data)[1], sum(is.na(df_merged[,get_identifier('drug')])), dim(df_raw_data)[1])
    print(WarnMsg)
    log_str = c(log_str, WarnMsg)

    # reorder the columns
    df_raw_data = Order_result_df(df_raw_data)

    return(df_raw_data)
}


#' @export
normalize_SE = function(df_raw_data, log_str, selected_keys = NULL,
                key_values = NULL, discard_keys = NULL) {
    # average technical replicates and assign the right controls to each treated well

    Keys = identify_keys(df_raw_data)
    Keys$discard_keys = discard_keys
    if (!is.null(selected_keys)) {
        Keys[names(selected_keys)] = selected_keys[names(selected_keys)]
    }
    if (!is.null(discard_keys)) {
      Keys$DoseResp = setdiff(Keys$DoseResp, discard_keys)
    }

    # the normalized SE only contains the treated conditions
    normSE = gDR::createSE(df_raw_data, data_type = "treated", discard_keys = discard_keys)
    SummarizedExperiment::assayNames(normSE) = 'Normalized'
    ctrlSE = gDR::createSE(df_raw_data, data_type = "untreated", discard_keys = discard_keys)

    # enforced key values for end points (override selected_keys) --> for rows of the SE
    Keys$untrt_Endpoint = setdiff(Keys$untrt_Endpoint, names(key_values))
    row_endpoint_value_filter = array(TRUE, nrow(ctrlSE))
    if (!is.null(key_values) & length(key_values)>0) {
        for (i in which(names(key_values) %in% names(SummarizedExperiment::rowData(ctrlSE)))) {
            if (is.numeric(key_values[i])) {
                row_endpoint_value_filter = row_endpoint_value_filter &
                    (SummarizedExperiment::rowData(ctrlSE)[, names(key_values)[i] ] == key_values[i] &
                            !is.na(SummarizedExperiment::rowData(ctrlSE)[, names(key_values)[i] ]))
            } else {
                row_endpoint_value_filter = row_endpoint_value_filter &
                    (SummarizedExperiment::rowData(ctrlSE)[ ,names(key_values)[i] ] %in% key_values[i])
            }}}

    #TODO gladkia: move mapping code to the separate function
    # perform the mapping for normalization
    # first the rows
    # matching the reference end point without any treatment
    row_maps_end = lapply(rownames(normSE), function(x) {
        # define matix with matching metadata
        match_mx = c(
            (SummarizedExperiment::rowData(ctrlSE) == (SummarizedExperiment::rowData(normSE)[x,]))[
                intersect(Keys$untrt_Endpoint,names(SummarizedExperiment::rowData(ctrlSE)))],
            IRanges::LogicalList(key_values = row_endpoint_value_filter,
                conc = apply(cbind(array(0, nrow(ctrlSE)),# padding to avoid empty df
                    SummarizedExperiment::rowData(ctrlSE)[,agrep('Concentration',
                    colnames(SummarizedExperiment::rowData(ctrlSE))),drop=F]),1,
                        function(x) all(x==0))))
        match_idx = which(apply(as.matrix(match_mx), 2, all))
        if (length(match_idx)==0) {
            # if not exact match, try to find best match
            WarnMsg = paste('Missing treated contols for:', x)
            idx = apply(as.matrix(match_mx), 2, function(y) sum(y, na.rm=T)) *
                            match_mx[[get_identifier('duration')]]
            if (any(idx>0)) {
                match_idx = which.max(idx)
                WarnMsg = paste(WarnMsg,'; found partial match:',
                        rownames(ctrlSE)[match_idx])
            } else WarnMsg = paste(WarnMsg,'; no partial match found')
            warning(WarnMsg)
        }
        return(rownames(ctrlSE)[match_idx])
    })
    names(row_maps_end) = rownames(normSE)

    # matching the reference end point with the same co-treatment (all the same but conc=0/Gnumber='vehicle')
    row_maps_cotrt = lapply(rownames(normSE), function(x)
        rownames(ctrlSE)[which(apply(as.matrix(
            (SummarizedExperiment::rowData(ctrlSE) == (SummarizedExperiment::rowData(normSE)[x,]))[
                intersect(Keys$ref_Endpoint,names(SummarizedExperiment::rowData(ctrlSE)))] ),
            2, all))])
    names(row_maps_cotrt) = rownames(normSE)
   
    #keep only valid row_maps_cotrt ---> that should not be necessary
    #row_maps_cotrt <-
     # row_maps_cotrt[vapply(row_maps_cotrt, length, numeric(1)) > 0]
    
    # matching the reference at time 0 (if available)
    row_maps_T0 = lapply(rownames(normSE), function(x) {
        # define matix with matching metadata
        match_mx = c(
            (SummarizedExperiment::rowData(ctrlSE) == (SummarizedExperiment::rowData(normSE)[x,]))[
                intersect(Keys$Day0,names(SummarizedExperiment::rowData(ctrlSE)))],
            IRanges::LogicalList(#key_values = row_endpoint_value_filter,
                T0 = SummarizedExperiment::rowData(ctrlSE)[, get_identifier('duration')] == 0,
                conc = apply(cbind(array(0, nrow(ctrlSE)),# padding to avoid empty df
                    SummarizedExperiment::rowData(ctrlSE)[,agrep('Concentration',
                    colnames(SummarizedExperiment::rowData(ctrlSE))),drop=F]),1,
                        function(x) all(x==0)) ))
        match_idx = which(apply(as.matrix(match_mx), 2, all))
        if (length(match_idx)==0) {
            # if not exact match, try to find best match
            WarnMsg = paste('Missing day 0 plate for:', x)
            idx = apply(as.matrix(match_mx), 2, function(y) sum(y, na.rm=T)) *
                            match_mx[['T0']]
            if (any(idx>0)) {
                match_idx = which.max(idx)
                WarnMsg = paste(WarnMsg,'; found partial match:',
                        rownames(ctrlSE)[match_idx])
            } else WarnMsg = paste(WarnMsg,'; no partial match found')
            warning(WarnMsg)
        }
        return(rownames(ctrlSE)[match_idx])
    })
    names(row_maps_T0) = rownames(normSE)

    # mapping for columns; 1 to 1 unless overridden by key_values
    col_maps = array(colnames(ctrlSE), dimnames = list(colnames(normSE)))
    if (any(names(key_values) %in% names(SummarizedExperiment::colData(normSE)))) {
        col_maps[] = colnames(ctrlSE)[
                which(key_values[names(key_values) %in% names(SummarizedExperiment::colData(normSE))] ==
                    SummarizedExperiment::colData(ctrlSE)[, names(SummarizedExperiment::colData(ctrlSE)) %in% names(key_values)])]
    }



    # remove background value to readout
    normSE = aapply(normSE, function(x) {
        x$CorrectedReadout = pmax(x$ReadoutValue - x$BackgroundValue,1)
        return(x)},
        'Normalized')
    ctrlSE = aapply(ctrlSE, function(x) {
        x$CorrectedReadout = pmax(x$ReadoutValue - x$BackgroundValue,1)
        return(x)})

    SummarizedExperiment::assay(normSE, 'Controls') <- matrix(lapply(1:prod(dim(normSE)), function(x) DataFrame()),
            nrow = nrow(normSE), ncol = ncol(normSE))

    # run through all conditions to assign controls and normalize the data
    # TODO: optimize (could that be replaced by a lapply?? or dyplr function??)
    for (i in rownames(normSE)) {
        for (j in colnames(normSE)) {

            if (nrow(SummarizedExperiment::assay(normSE,'Normalized')[[i,j]]) == 0) next # skip if no data

            df_end = do.call(rbind,
                    lapply(row_maps_end[[i]], function(x) SummarizedExperiment::assay(ctrlSE)[[x, col_maps[j]]]))
            df_end = df_end[, c('CorrectedReadout',
                    intersect(Keys$untrt_Endpoint,colnames(df_end)))]
            colnames(df_end)[1] = 'UntrtReadout'
            df_end = aggregate(df_end[,1,drop=F], by = as.list(df_end[,-1,drop=F]),
                function(x) mean(x, trim= .25))

            # not always present
            if (i %in% names(row_maps_cotrt)) {
                df_ref = do.call(rbind,
                        lapply(row_maps_cotrt[[i]], function(x) SummarizedExperiment::assay(ctrlSE)[[x, col_maps[j]]]))
                df_ref = df_ref[, c('CorrectedReadout',
                        intersect(Keys$ref_Endpoint,colnames(df_ref)))]
                colnames(df_ref)[1] = 'RefReadout'
                df_ref = aggregate(df_ref[,1,drop=F], by = as.list(df_ref[,-1,drop=F]),
                    function(x) mean(x, trim= .25))

                # check if control and co-treated wells are on the same plate
                if (all(df_end$Barcode %in% df_ref$Barcode) && all(df_ref$Barcode %in% df_end$Barcode)) {
                  df_end = merge(df_end, df_ref, by=c('Barcode', Keys$discard_keys))
                } else {
                  wrnMsg1 <-
                    sprintf(
                      "Control data for the drug are propagated to other plates with co-drug controls.
                      Treatment Id: %s
                      Cell_line Id: %s",
                      i,
                      j
                    )
                  warning(wrnMsg1)
                  # propagate average values to the other plates
                  df_end = merge(df_end, df_ref, by='Barcode', all=T)
                  mean_UntrtReadout = mean(df_end$UntrtReadout, na.rm=T)
                  mean_RefReadout = mean(df_end$RefReadout, na.rm=T)
                  df_end$UntrtReadout[is.na(df_end$UntrtReadout)] = mean_UntrtReadout
                  df_end$RefReadout[is.na(df_end$RefReadout)] = mean_RefReadout
                }
                
                #gladkia: assert for control data
                if (nrow(df_end) == 0) {
                  errMsg1 <-
                    sprintf(
                      "Control dataframe failed.
                      Treatment Id: '%s'
                      Cell_line Id: %s",
                      i,
                      j
                    )
                  stop(errMsg1)
                }
            } else {
                df_end$RefReadout = df_end$UntrtReadout
            }

            df_0 = do.call(rbind,
                    lapply(row_maps_T0[[i]], function(x) SummarizedExperiment::assay(ctrlSE)[[x, col_maps[j]]]))
            df_0 = df_0[, c('CorrectedReadout', intersect(Keys$Day0,colnames(df_0)))]
            colnames(df_0)[1] = 'Day0Readout'
            df_0 = aggregate(df_0[,1,drop=F], by = as.list(df_0[,-1,drop=F]),
                function(x) mean(x, trim= .25))

            if (!is.null(Keys$discard_keys) && all(Keys$discard_keys %in% colnames(df_0))) {
              df_ctrl = merge(df_0[, setdiff(colnames(df_0), 'Barcode')], df_end, all.y = T, by = Keys$discard_keys)
            } else {
              df_ctrl = merge(df_0[, setdiff(colnames(df_0), 'Barcode')], df_end, all.y = T)
              colnames(df_ctrl)[1] = 'Day0Readout'
            }

            df_ctrl$RefRelativeViability = round(df_ctrl$RefReadout/df_ctrl$UntrtReadout,4)

            df_ctrl$RefGRvalue = round(2 ** (
                    log2(df_ctrl$RefReadout / df_ctrl$Day0Readout) /
                    log2(df_ctrl$UntrtReadout / df_ctrl$Day0Readout) ), 4) - 1

            df_ctrl$DivisionTime = round(
                    SummarizedExperiment::rowData(normSE)[i,get_identifier('duration')] /
                        log2(df_ctrl$UntrtReadout / df_ctrl$Day0Readout) , 4)
            SummarizedExperiment::assay(normSE, 'Controls')[[i,j]] = DataFrame(df_ctrl)

            #gladkia: assert for merged study/control data
            ctrl_bcodes <- sort(unique(df_ctrl$Barcode))
            trt_bcodes <-
              sort(unique(SummarizedExperiment::assay(normSE, 'Normalized')[[i, j]]$Barcode))
            if (!all(trt_bcodes %in% ctrl_bcodes)) {
              wrnMsg1 <-
                sprintf(
                  "Control data are averaged and propagated to treatment plates.
                      Treatment Id: %s (plates %s)
                      Control plates: %s",
                  i, paste(trt_bcodes, collapse = ', '),
                  paste(ctrl_bcodes, collapse = ', ')
                )
              warning(wrnMsg1)
              bind_rows(df_ctrl, cbind(data.frame(Barcode = setdiff(trt_bcodes, ctrl_bcodes)),
                        t(colMeans(df_ctrl[, setdiff(colnames(df_ctrl), 'Barcode')]))))
            }
            
            # merge the data with the controls assuring that the order of the records is preseved
            df_merged = dplyr::left_join(data.frame(SummarizedExperiment::assay(normSE, 'Normalized')[[i,j]]),
                    data.frame(df_ctrl), by = c('Barcode', Keys$discard_keys))

            # calculate the normalized values
            
            SummarizedExperiment::assay(normSE, 'Normalized')[[i,j]]$RelativeViability =
                round(df_merged$CorrectedReadout/df_merged$UntrtReadout,4)
            
            SummarizedExperiment::assay(normSE, 'Normalized')[[i,j]]$GRvalue = round(2 ** (
                    log2(df_merged$CorrectedReadout / df_merged$Day0Readout) /
                    log2(df_merged$UntrtReadout / df_merged$Day0Readout) ), 4) - 1
           
        }
    }
    metadata(normSE) = c(metadata(normSE),
            list(df_raw_data = df_raw_data,
                Keys = Keys,
                row_maps = list(end = row_maps_end,
                                cotrt = row_maps_cotrt,
                                T0 = row_maps_T0)
                ))

    return(normSE)
}

#' @export
normalize_data = function(df_raw_data, log_str, selected_keys = NULL,
                key_values = NULL) {
    # average technical replicates and assign the right controls to each treated well

    # remove unused columns but keep barcodes to normalize by plate
    df_normalized = df_raw_data[, setdiff(colnames(df_raw_data),
                        c('Template', get_identifier('WellPosition')) ) ]

    # Identify keys for assigning the controls
    Keys = identify_keys(df_normalized)
    if (!is.null(selected_keys)) {
        Keys[names(selected_keys)] = selected_keys[names(selected_keys)]
    }

    df_normalized$CorrectedReadout = pmax(df_normalized$ReadoutValue -
                                            df_normalized$BackgroundValue,1)

    # enforced key values for end points (override selected_keys)
    Keys$untrt_Endpoint = setdiff(Keys$untrt_Endpoint, names(key_values))
    endpoint_value_filter = array(TRUE, dim(df_raw_data)[1])
    if (!is.null(key_values) & length(key_values)>0) {
        for (i in 1:length(key_values)) {
            if (is.numeric(key_values[i])) {
                endpoint_value_filter = endpoint_value_filter &
                            (df_normalized[, names(key_values)[i] ] == key_values[i] &
                                !is.na(df_normalized[, names(key_values)[i] ]))
            } else {
                endpoint_value_filter = endpoint_value_filter &
                            (df_normalized[ ,names(key_values)[i] ] %in% key_values[i])
            }}}
    # get the untreated controls at endpoint and perform interquartile mean
    df_end_untrt = df_normalized[df_normalized[,get_identifier('duration')]>0 & endpoint_value_filter &
        apply(df_normalized[,agrep('Concentration', colnames(df_normalized)),drop=F]==0,1,all),]
    df_end_mean = aggregate(df_end_untrt[,'CorrectedReadout'],
                    by = as.list(df_end_untrt[,Keys$untrt_Endpoint]), function(x) mean(x, trim= .25))
    colnames(df_end_mean)[dim(df_end_mean)[2]] = 'UntrtReadout'


    # get the untreated controls at Day 0 and perform interquartile mean
    df_day0 = df_normalized[df_normalized[,get_identifier('duration')]==0 &
        apply(df_normalized[,agrep('Concentration', colnames(df_normalized)),drop=F]==0,1,all), ]
    df_day0_mean = aggregate(df_day0[,'CorrectedReadout'],
                by = as.list(df_day0[,Keys$Day0]), function(x) mean(x, trim= .25))
    colnames(df_day0_mean)[dim(df_day0_mean)[2]] = 'Day0Readout'

    df_controls = merge(df_end_mean, df_day0_mean[, setdiff(colnames(df_day0_mean),
                    c(get_identifier('duration'), 'Barcode'))], all.x = T)
    if (length(setdiff(Keys$untrt_Endpoint, Keys$Day0))>0) {
        WarnMsg = paste('Not all control conditions found on the day 0 plate,',
            'dispatching values for field: ',
            paste(setdiff(Keys$untrt_Endpoint, Keys$Day0), collapse = ' ; '))
        log_str = c(log_str, 'Warning in normalize_data:')
        log_str = c(log_str, WarnMsg)
        warning(WarnMsg)
    }
    # identify missing values in the Day0 that needs to be matched (usually for co-treatments)
    df_controls_NA = which(is.na(df_controls$Day0Readout))

    if (length(df_controls_NA)>0) {
        dispatched = NULL
        for (i in df_controls_NA){
            matches = t(apply(df_day0_mean[, setdiff(Keys$Day0,
                                                c(get_identifier('duration'), 'Barcode'))], 1,
                function(x) df_controls[i, setdiff(Keys$Day0, c(get_identifier('duration'), 'Barcode')),
                                drop=F] == x))
            colnames(matches) = setdiff(Keys$Day0, c(get_identifier('duration'), 'Barcode'))
            # try to find a good match for the day 0 (enforce same cell line)
            idx = rowSums(matches) * matches[, get_identifier('cellline')]
            if (all(idx==0)) {next}
            match_idx = which.max(idx)
            mismatch = df_day0_mean[match_idx, setdiff(Keys$Day0,
                                                c(get_identifier('duration'), 'Barcode'))] !=
                        df_controls[i, setdiff(Keys$Day0, c(get_identifier('duration'), 'Barcode'))]
            dispatched = c(dispatched, colnames(mismatch)[mismatch])
            df_controls[i, 'Day0Readout'] = df_day0_mean[match_idx, 'Day0Readout']
        }
        WarnMsg = paste('Not all control conditions found on the day 0 plate,')
        WarnMsg = ifelse(length(dispatched)>0,
                    paste(WarnMsg,'dispatching values for mismatches in field: ',
                    paste(unique(dispatched), collapse = ' ; ')),
                    paste(WarnMsg,'some Day0 are not being matched'))
        log_str = c(log_str, 'Warning in normalize_data:')
        log_str = c(log_str, WarnMsg)
        warning(WarnMsg)
    }

    df_to_norm = df_normalized[df_normalized[,get_identifier('duration')]>0 &
        (apply(df_normalized[,agrep('Concentration', colnames(df_normalized)),drop=F]!=0,1,any) |
                !endpoint_value_filter),]

    df_to_norm_conditions = unique(df_to_norm[, intersect(colnames(df_to_norm),
                                colnames(df_controls))])

    # if missing barcodes --> dispatch for similar conditions
    if (!all(df_to_norm_conditions$Barcode %in% df_controls$Barcode)) {
        WarnMsg = paste('Not all control conditions found at the end of treatment,')
        WarnMsg = paste(WarnMsg,'dispatching values for plates: ',
                    paste(setdiff(df_to_norm_conditions$Barcode, df_controls$Barcode),
                        collapse = ' ; '))
        log_str = c(log_str, 'Warning in normalize_data:')
        log_str = c(log_str, WarnMsg)
        warning(WarnMsg)

        df_ctrl_mean = aggregate(df_controls[,c('UntrtReadout', 'Day0Readout')],
                    by = as.list(subset(df_controls,
                        select = -c(UntrtReadout,Day0Readout,Barcode))), mean)
        df_controls = rbind(df_controls, merge(df_ctrl_mean,
                df_to_norm_conditions[!(df_to_norm_conditions$Barcode %in%
                        df_controls$Barcode),]))
    }

    df_normalized = merge(df_to_norm, df_controls)

    df_normalized$RelativeViability = round(df_normalized$CorrectedReadout/
                                            df_normalized$UntrtReadout,4)
    df_normalized$GRvalue = round(2 ** (
            log2(df_normalized$CorrectedReadout / df_normalized$Day0Readout) /
                log2(df_normalized$UntrtReadout / df_normalized$Day0Readout) ), 4) - 1

    df_normalized$DivisionTime = round( df_normalized[,get_identifier('duration')] /
                    log2(df_normalized$UntrtReadout / df_normalized$Day0Readout) , 4)


    if (any(is.na(df_normalized$Day0Readout))) {
        # need to use the reference doubling Time if day 0 missing
        InferedIdx = is.na(df_normalized$Day0Readout)
        filtered = df_normalized$ReferenceDivisionTime > (df_normalized[,get_identifier('duration')]*2) |
            is.na(df_normalized$ReferenceDivisionTime)
        WarnMsg = paste('Missing day 0 information --> calculate GR value based on reference',
            'doubling time')
        WarnMsg = ifelse(!any(filtered & InferedIdx), WarnMsg,
            paste(WarnMsg, '; filtering', sum(filtered & InferedIdx),
                'conditions because of too short assay:',
                paste(unique(df_normalized$CellLineName[filtered & InferedIdx]), collpase=' ; ')))
        log_str = c(log_str, 'Warning in normalize_data:')
        log_str = c(log_str, WarnMsg)
        warning(WarnMsg)
        InferedIdx = !filtered & InferedIdx
        # calculate GR values using formula from https://www.nature.com/articles/nbt.3882
        df_normalized$GRvalue[InferedIdx] = round(2 ^ ( 1 + (log2(pmin(1.25,
                df_normalized[InferedIdx,'RelativeViability'])) /
                (df_normalized[,get_identifier('duration')][InferedIdx] /
                        df_normalized$ReferenceDivisionTime[InferedIdx]) ) ),4) - 1
    }

    df_normalized = cbind(df_normalized[, 1:(which(colnames(df_normalized)=='ReadoutValue')-1)],
        df_normalized[, c('GRvalue', 'RelativeViability', 'DivisionTime')],
        df_normalized[, which(colnames(df_normalized)=='ReadoutValue'):(dim(df_normalized)[2]-3)])
    df_normalized = Order_result_df(df_normalized)
    print('df normalized:')
    print(head(df_normalized))
    return(df_normalized)
}



#' @export
average_SE = function(normSE, TrtKeys = NULL) {
  
    avgSE = normSE
    if (is.null(TrtKeys)) {
        if ('Keys' %in% names(metadata(normSE))) {
          TrtKeys = metadata(normSE)$Keys$Trt
          TrtKeys = setdiff(TrtKeys, metadata(normSE)$Keys$discard_keys)
        } else TrtKeys = identify_keys(normSE)$Trt
    }
    metadata(normSE)$Keys$Trt = TrtKeys
      
    SummarizedExperiment::assay(avgSE, 'Averaged') = SummarizedExperiment::assay(avgSE, 'Normalized')
    avgSE = aapply(avgSE, function(x) {
        if (nrow(x) > 1) {
            subKeys = intersect(TrtKeys, colnames(x))
            df_av = aggregate(x[, c('GRvalue', 'RelativeViability',"CorrectedReadout")],
                            by = as.list(x[,subKeys,drop=F]), FUN = function(y) mean(y, na.rm=T))
            df_std = aggregate(x[, c('GRvalue', 'RelativeViability')],
                                by = as.list(x[,subKeys,drop=F]), FUN = function(x) sd(x, na.rm=T))
            colnames(df_std)[colnames(df_std) %in% c('GRvalue', 'RelativeViability')] =
                paste0('std_',
                    colnames(df_std)[colnames(df_std) %in% c('GRvalue', 'RelativeViability')])
            return( merge(df_av, df_std, by = subKeys) )
        } else return(x)
    }, 'Averaged')

    SummarizedExperiment::assay(avgSE, 'Avg_Controls') = SummarizedExperiment::assay(avgSE, 'Controls')
    avgSE = aapply(avgSE, function(x) {
        if (nrow(x) > 1) {
            subKeys = intersect(TrtKeys, colnames(x))
            df_av = DataFrame(lapply(x[, c('Day0Readout', 'UntrtReadout',
                    'RefGRvalue', 'RefRelativeViability',
                    'RefReadout', 'DivisionTime')], FUN = function(y) mean(y, na.rm=T)))
            return( df_av )
        } else return(x)
    }, 'Avg_Controls')
    
    return(avgSE)
}

#' @export
average_replicates = function(df_normalized, TrtKeys = NULL) {
    if (is.null(TrtKeys)) { TrtKeys = identify_keys(df_normalized)$Trt }

    df_averaged = aggregate(df_normalized[, c('GRvalue', 'RelativeViability', "CorrectedReadout",
                    "UntrtReadout", "Day0Readout", "DivisionTime", "ReferenceDivisionTime")],
                    by = as.list(df_normalized[,TrtKeys]), FUN = function(x) mean(x, na.rm=T))
    df_std = aggregate(df_normalized[, c('GRvalue', 'RelativeViability')],
                        by = as.list(df_normalized[,TrtKeys]), FUN = function(x) sd(x, na.rm=T))
    colnames(df_std)[colnames(df_std) %in% c('GRvalue', 'RelativeViability')] = paste0('std_',
            colnames(df_std)[colnames(df_std) %in% c('GRvalue', 'RelativeViability')])
    df_averaged = merge(df_averaged, df_std, by = TrtKeys)

    #reorganize column order:
    df_averaged = Order_result_df(df_averaged)

    print('df averaged:')
    print(head(df_averaged))
    return(df_averaged)
}


#' @export
metrics_SE = function(avgSE, studyConcThresh = 4) {

    # this is not used as we enforce the same conditions as the input SE; not collapsing allowed
    # if (is.null(DoseRespKeys)) {
    #     if ('Keys' %in% names(metadata(avgSE))) DoseResp = metadata(avgSE)$Keys$DoseResp
    #     else DoseRespKeys = identify_keys(avgSE)$DoseResp
    # } else {
    #     metadata(avgSE)$Keys$DoseResp = DoseRespKeys
    # }

    metricsSE = avgSE
    SummarizedExperiment::assay(metricsSE, 'Metrics') = SummarizedExperiment::assay(metricsSE, 'Averaged')

    for (i in rownames(metricsSE)) {
        for (j in colnames(metricsSE)) {
            df_ = SummarizedExperiment::assay(metricsSE, 'Averaged')[[i,j]]
            if (length(unique(df_$Concentration)) >= studyConcThresh) {
                SummarizedExperiment::assay(metricsSE, 'Metrics')[[i,j]] = DataFrame(ICGRfits(df_,
                    e_0 = SummarizedExperiment::assay(metricsSE, 'Avg_Controls')[[i,j]]$RefRelativeViability,
                    GR_0 = SummarizedExperiment::assay(metricsSE, 'Avg_Controls')[[i,j]]$RefGRvalue))
            } else if (nrow(df_) == 0) {
                out = DataFrame(matrix(NA, 0, length(get_header('response_metrics'))+2))
                colnames(out) = c(get_header('response_metrics'), 'maxlog10Concentration', 'N_conc')
                SummarizedExperiment::assay(metricsSE, 'Metrics')[[i,j]] = out
            } else {
                out = DataFrame(matrix(NA, 2, length(get_header('response_metrics'))))
                colnames(out) = get_header('response_metrics')
                out$maxlog10Concentration = max(log10(df_$Concentration))
                out$N_conc = length(unique(df_$Concentration))
                SummarizedExperiment::assay(metricsSE, 'Metrics')[[i,j]] = out
            }
        }

    }
    return(metricsSE)
}



#' @importFrom dplyr arrange_at group_by_at left_join summarise
#' @export
calculate_DRmetrics <-
  function(df_averaged,
           DoseRespKeys = NULL,
           studyConcThresh = 4) {

   df_a = df_averaged
   colnames(df_a)[ colnames(df_a) == get_identifier('drugname') ] = 'DrugName'

   if (is.null(DoseRespKeys)) {
       DoseRespKeys = identify_keys(df_a)$DoseResp
   } else {
       DoseRespKeys [ DoseRespKeys == get_identifier('drugname') ] = 'DrugName'
   }
   DoseRespKeys = setdiff(DoseRespKeys, 'Concentration')
   DoseRespKeys = c(DoseRespKeys, 'DivisionTime')
   DoseRespKeys = intersect(DoseRespKeys, colnames(df_a))

   df_a$log10Concentration = log10(df_a$Concentration)

   metrics = names(ICGRlogisticFit(c(-7, -6, -5, -4), c(1, .9, .8, .7), c(1, .9, .8, .7)))
    # dummy call to get variable names

    #define set of key for merging control and study data
    mergeKeys <- setdiff(DoseRespKeys, c(get_identifier('drug'), 'DrugName'))

    #get avereage GRvalue ('GR_0') for control data
    controlSets <-
      df_a %>% filter(DrugName %in% get_identifier('untreated_tag')) %>%
      dplyr::group_by_at(mergeKeys) %>%
        dplyr::summarise(GR_0 = mean(GRvalue), e_0 = mean(RelativeViability))

    #get study data
    studySets <-
      df_a %>% filter(!DrugName %in% get_identifier('untreated_tag'))

    #join study and control data
    # i.e. get  reference (average control) GRvalue ('GR_0') for study data
    fSets <-
      dplyr::left_join(studySets, controlSets, by = mergeKeys)
    # for study sets with no reference GRvalue, assing GRValue0 to 1
    fSets[is.na(fSets$GR_0), "GR_0"] <- 1
    fSets[is.na(fSets$e_0), "e_0"] <- 1

    #group study data by 'DoseRespKeys'
    gSets <- fSets %>% dplyr::group_by_at(DoseRespKeys) %>% group_split()

    # filter to have at least 4 records with non-NA RelativeViability
    gSets = gSets[ lapply(gSets, function(x) sum(!is.na(x$RelativeViability))) >= studyConcThresh ]

    print(paste(
      'Metadata variables for dose response curves:',
      paste(setdiff(
        DoseRespKeys, c(get_identifier('drug'), get_identifier('cellline'), paste(get_identifier('drug'),'_', 2:10))
      ),
      collapse = ' '),
      '(',
      length(gSets),
      'groups )'
    ))


    #iterate over study groups
    resL <- lapply(1:length(gSets), function(x) {
      # the 'DoseRespKeys' columns in given grup are identical for each entry
      # let's get the first record then
      repCols <- as.vector(gSets[[x]][1, DoseRespKeys])
      #get selected columns ('metrics') from GRlogisticFit output
      # (if at least 4 records with non-NA RelativeViability)
      if (sum(!is.na(gSets[[x]]$RelativeViability)) >= studyConcThresh) {
        grLogCols <-
          ICGRlogisticFit(gSets[[x]]$log10Concentration,
                        gSets[[x]]$RelativeViability,
                        gSets[[x]]$GRvalue,
                        e_0 = gSets[[x]]$e_0[1],
                        GR_0 = gSets[[x]]$GR_0[1])[metrics]
      } else {
        grLogCols <- rep(NA, length(metrics))
        names(grLogCols) = metrics
        grLogCols$N_conc = sum(!is.na(gSets[[x]]$RelativeViability))
      }
      cbind(repCols, t(grLogCols))
    })

    #return final data.frame
    resDf <- do.call(rbind, resL)
    resDf = resDf [resDf$N_conc >= studyConcThresh,]
    resDf = resDf %>% dplyr::arrange_at(DoseRespKeys)
    colnames(resDf)[ colnames(resDf) == 'DrugName'] = get_identifier('drugname')
    resDf = Order_result_df(resDf)
    return(resDf)
  }

#' @export
identify_keys = function(df_se_mae) {

    if (any(class(df_se_mae) %in% c('MultiAssayExperiment', 'SummarizedExperiment'))) {
        if ('MultiAssayExperiment' %in% class(df_se_mae)) {
            # if MAE, convert to SE based on the treated SE (could be optimized)
            df_se_mae = df_se_mae[['treated']]
            se_untrt =  df_se_mae[['untreated']]
        } else se_untrt = NULL
        all_keys = unique(c(
            colnames(SummarizedExperiment::rowData(df_se_mae)),
            colnames(SummarizedExperiment::colData(df_se_mae)),
            unlist(lapply(SummarizedExperiment::assay(df_se_mae), colnames))))
    } else { # case of a data frame
        all_keys = colnames(df_se_mae)
    }

    keys = list(Trt = setdiff(all_keys, "Barcode"),
            DoseResp = setdiff(all_keys,  'Barcode'),
            ref_Endpoint = setdiff(all_keys, c('Concentration',
                                            get_identifier('drug'),
                                            get_identifier('drugname'))),
            untrt_Endpoint = all_keys[ c(-agrep('Concentration', all_keys),
                                            -agrep(get_identifier('drug'), all_keys),
                                            -agrep(get_identifier('drugname'), all_keys))])
    keys[['Day0']] = setdiff(keys[['untrt_Endpoint']], get_identifier('duration'))
    keys = lapply(keys, function(x) setdiff(x, c(get_header('raw_data'),
        get_header('normalized_results'), "Template", get_identifier('WellPosition'), get_header('averaged_results'),
            get_header('metrics_results'), "ReferenceDivisionTime"
    )))
    keys = lapply(keys, sort)

    # check if all values of a key is NA
    for (k in keys[['untrt_Endpoint']]) {

        if ('SummarizedExperiment' %in% class(df_se_mae)) {
            # check the metadata fields for NA
            if (k %in% colnames(SummarizedExperiment::rowData(df_se_mae))) df_ = SummarizedExperiment::rowData(df_se_mae)
            else if (k %in% colnames(SummarizedExperiment::colData(df_se_mae))) df_ = SummarizedExperiment::colData(df_se_mae)
            else next # not a metadata

            if (all(is.na(df_[,k]))) keys = lapply(keys, function(x) setdiff(x, k))

            if (!is.null(se_untrt) && k %in% colnames(SummarizedExperiment::rowData(se_untrt))) {
                df_ = SummarizedExperiment::rowData(se_untrt)
                if (all(is.na(df_[df_[,get_identifier('duration')]==0,k]))) {
                    keys[['Day0']] = setdiff(keys[['Day0']], k)
                }
            }
        } else { # case of a data frame
            if (all(is.na(df_se_mae[,k]))) {
                keys = lapply(keys, function(x) setdiff(x, k))
            }
            if (all(is.na(df_se_mae[df_se_mae[,get_identifier('duration')]==0,k]))) {
                keys[['Day0']] = setdiff(keys[['Day0']], k)
            }
        }
    }
    return(keys)
}


#' @export
#' @importFrom gneDB annotateCLIDs
#' @importFrom gCellGenomics getDrugs
cleanup_metadata = function(df_metadata, log_str) {
    log_str = c(log_str, '    cleanup_metadata')
    # clean up numberic fields
    df_metadata[,get_identifier('duration')] = round(as.numeric(df_metadata[,get_identifier('duration')]),6)
    # identify potential numeric fields and replace NA by 0 - convert strings in factors
    for (c in setdiff(1:dim(df_metadata)[2], c( agrep(get_identifier('drug'), colnames(df_metadata)),
            agrep('Concentration', colnames(df_metadata)),
            grep(paste(c(get_identifier('cellline'), get_header('manifest'), get_identifier('WellPosition')),
                collapse = '|'), colnames(df_metadata)) ))) {
        vals = unique(df_metadata[,c])

        if (is.character(vals)) {
            num_vals = as.numeric(vals)
            if (sum(is.na(num_vals))>2 || all(is.na(num_vals))) {
                df_metadata[,c] = factor(df_metadata[,c])
                WarnMsg = paste('Metadata field ', colnames(df_metadata)[c],
                                    ' converted to factors')
                log_str = c(log_str, 'Warning in cleanup_metadata:')
                log_str = c(log_str, WarnMsg)
                warning(WarnMsg)
            } else {
                is.na(df_metadata[,c]) = 0
                df_metadata[,c] = as.numeric(df_metadata[,c])
                WarnMsg = paste('Metadata field ', colnames(df_metadata)[c],
                                    ' converted to numeric values')
                log_str = c(log_str, 'Warning in cleanup_metadata:')
                log_str = c(log_str, WarnMsg)
                warning(WarnMsg)
            }
        }
    }

    # TODO: specific to GNE database --> need to be replaced by a function
    df_metadata = gDR::add_CellLine_annotation(df_metadata)

    print(dim(df_metadata))

    # check that Gnumber_* are in the format 'G####' and add common name (or Vehicle or Untreated)

    for (i in agrep(get_identifier('drug'), colnames(df_metadata))) { # correct case issues
        for (w in get_identifier('untreated_tag')) {
            df_metadata[grep(w, df_metadata[,i], ignore.case = T),i] = w
        }
    }
    # -----------------------

    df_metadata = add_Drug_annotation(df_metadata)

    # clean up concentration fields
    for (i in agrep('Concentration', colnames(df_metadata))) {
        trt_n = ifelse(regexpr('_\\d', colnames(df_metadata)[i])>0,
                            substr(colnames(df_metadata)[i], 15,20), 1)
        DrugID_col = ifelse(trt_n == 1, get_identifier('drug'), paste0(get_identifier('drug'), '_', trt_n))
        df_metadata[df_metadata[,DrugID_col] %in% get_identifier('untreated_tag'),i] = 0 # set all untreated to 0

        DrugID_0 = setdiff(unique(df_metadata[ df_metadata[,i] == 0, DrugID_col]), get_identifier('untreated_tag'))
        DrugID_0 = DrugID_0[!is.na(DrugID_0)]
        if (length(DrugID_0)>0) {
            WarnMsg = paste('Some concentration for ', DrugID_col,
                            ' are 0: ', paste(DrugID_0, collapse = ' ; '))
            log_str = c(log_str, 'Warning in cleanup_metadata:')
            log_str = c(log_str, WarnMsg)
            warning(WarnMsg)
        }
        df_metadata[,i] = round(as.numeric(df_metadata[,i]),6) # avoid mismatch due to string truncation
    }

    return(df_metadata)
}


#' @export
Order_result_df = function (df_) {
    cols = c(get_header('ordered_1'), setdiff(colnames(df_),
        c(get_header('ordered_1'), get_header('ordered_2'))), get_header('ordered_2'))
    cols = intersect(cols, colnames(df_))

    row_order_col = intersect( c(get_header('add_clid')[1], get_identifier('duration'),
            get_identifier('drugname'), 'Concentration',
            paste0(c(paste0(get_identifier('drugname'),'_'),'Concentration_'),
                        sort(c(2:10,2:10))),
            setdiff(colnames(df_), c(get_header('ordered_1'), get_header('ordered_2'))) ),
            cols)

    df_ = df_[do.call(order, df_[,row_order_col]), cols]

    return(df_)
}


#
add_CellLine_annotation = function(df_metadata) {

    DB_cellid_header = 'clid'
    DB_cell_annotate = c('celllinename', 'primarytissue', 'doublingtime')
    # corresponds to columns get_header('add_clid'): name, tissue, doubling time
    CLs_info = tryCatch( {
        CLs_info = gneDB::annotateCLIDs(unique(df_metadata[,get_identifier('cellline')]))
        CLs_info = CLs_info[,c(DB_cellid_header,DB_cell_annotate)]
        CLs_info
    }, error = function(e) {
        print('failed to load cell line info from DB')
        print(e)
        data.frame()
    })

    if (nrow(CLs_info)==0) return(df_metadata)

    colnames(CLs_info) = c(get_identifier('cellline'), get_header('add_clid'))
    CLIDs = unique(df_metadata[,get_identifier('cellline')])
    bad_CL = !(CLIDs %in% CLs_info[,get_identifier('cellline')])
    if (any(bad_CL)) {
        ErrorMsg = paste('Cell line ID ', paste(CLIDs[bad_CL], collapse = ' ; '),
            ' not found in cell line database')
        log_str = c(log_str, 'Error in cleanup_metadata:')
        log_str = c(log_str, ErrorMsg)

        stop(ErrorMsg)
        }

    print('merge with Cell line info')
    nrows_df = nrow(df_metadata)
    df_metadata = merge(df_metadata, CLs_info, by.x = get_identifier('cellline'),
                by.y = DB_cellid_header, all.x = T)
    stopifnot(nrows_df == nrow(df_metadata))

    return(df_metadata)

}


add_Drug_annotation = function(df_metadata) {
        nrows_df = nrow(df_metadata)

        DB_drug_identifier = 'drug'
        Drug_info = tryCatch( {
                gDrugs = gCellGenomics::getDrugs()[,c(DB_drug_identifier, 'gcsi_drug_name')]
                gDrugs[,1] = substr(gDrugs[,1], 1, 9) # remove batch number from DB_drug_identifier
                gDrugs
        }, error = function(e) {
            print('failed to load drug info from DB')
            print(e)
            data.frame()
        })

        if (nrow(Drug_info) == 0) {
            df_metadata[, get_identifier('drugname')] = df_metadata[,get_identifier('drug')]
            return(df_metadata)
        }
        # -----------------------

        colnames(Drug_info) = c('drug', 'DrugName')
        Drug_info = rbind(data.frame(drug=get_identifier('untreated_tag'), DrugName=get_identifier('untreated_tag')), Drug_info)
        Drug_info = unique(Drug_info)
        DrIDs = unique(unlist(df_metadata[,agrep(get_identifier('drug'), colnames(df_metadata))]))
        bad_DrID = !(DrIDs %in% Drug_info$drug) & !is.na(DrIDs)
        if (any(bad_DrID)) {
            # G number, but not registered
            ok_DrID = attr(regexpr('^G\\d*',DrIDs), 'match.length')==9
            if (any(ok_DrID)) {
                WarnMsg = paste('Drug ', paste(DrIDs[ok_DrID & bad_DrID], collapse = ' ; '),
                    ' not found in gCSI database; use G# as DrugName')
                Drug_info = rbind(Drug_info, data.frame(drug=DrIDs[ok_DrID & bad_DrID],
                        DrugName=DrIDs[ok_DrID & bad_DrID]))
                log_str = c(log_str, 'Warning in cleanup_metadata:')
                log_str = c(log_str, WarnMsg)
                warning(WarnMsg)
            } else {
                ErrorMsg = paste('Drug ', paste(DrIDs[!ok_DrID], collapse = ' ; '),
                    ' not found in gCSI database')
                log_str = c(log_str, 'Error in cleanup_metadata:')
                log_str = c(log_str, ErrorMsg)

                stop(ErrorMsg)
            }
        }
        colnames(Drug_info)[2] = get_identifier('drugname')
        print('merge with Drug_info for Drug 1')
        df_metadata = merge(df_metadata, Drug_info, by.x=get_identifier('drug'), by.y='drug', all.x = T)
        print(dim(df_metadata))
        # add info for columns Gnumber_*
        for (i in grep(paste0(get_identifier('drug'),'_\\d'), colnames(df_metadata))) {
            df_metadata[ is.na(df_metadata[,i]), i] = get_identifier('untreated_tag')[1] # set missing values to Untreated
            Drug_info_ = Drug_info
            colnames(Drug_info_)[2] = paste0(colnames(Drug_info_)[2], substr(colnames(df_metadata)[i], 8, 12))
            print(paste('merge with Drug_info for ', i))
            df_metadata = merge(df_metadata, Drug_info_, by.x=i, by.y='drug', all.x = T)
            print(dim(df_metadata))
        }
        df_metadata[, colnames(df_metadata)[grepl(get_identifier('drugname'), colnames(df_metadata))] ] =
            droplevels(df_metadata[,colnames(df_metadata)[grepl(get_identifier('drugname'), colnames(df_metadata))]])

    stopifnot(nrows_df == nrow(df_metadata))

    return(df_metadata)

}
