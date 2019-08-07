
#' @import gneDB
#' @import reshape2
#' @import dplyr

#########################################
### TODO:
### Need to put all reserved header names (CLID, DrugName, ...) as global variables
### Need to put all output header names (GRvalue, RelativeViability, ...) as global variables
### see identify_keys for all names to protect
#########################################



#' @export
Overall_function = function(manifest_file, template_file, results_file,
                output_files, selected_keys = NULL, key_values = NULL) {
    # output_files should contain file names for :
    #   log_file, QC_file, raw_result, process_results, metrics_results

    log_str = 'Report from gDR pipeline'

    raw_data = load_data(manifest_file, template_file, results_file, log_str)
    df_merged_data = merge_data(raw_data$manifest, raw_data$treatments, raw_data$data, log_str)

    # output_QC_byPlate(df_raw_data, output_files['QC_file']) # TODO: check column/row bias

    Keys = identify_keys(df_raw_data) # may be manually changed
    if (!is.null(selected_keys)) {
        Keys[names(selected_keys)] = selected_keys[names(selected_keys)]
    }

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
    expected_headers = cellline_identifier
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

    # clean up the metadata
    print(colnames(df_metadata))
    print(dim(df_metadata))
    cleanedup_metadata = cleanup_metadata(df_metadata, log_str)
    print(colnames(cleanedup_metadata))
    print(dim(cleanedup_metadata))
    stopifnot( dim(cleanedup_metadata)[1] == dim(df_metadata)[1] ) # should not happen

    df_merged = merge(cleanedup_metadata, data, by = c('Barcode', 'WellRow', 'WellColumn'))
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
    df_raw_data = df_merged[ !is.na(df_merged[,drug_identifier]), ]
    WarnMsg = sprintf('%i well loaded, %i discarded for lack of annotation, %i data point selected',
                dim(data)[1], sum(is.na(df_merged[,drug_identifier])), dim(df_raw_data)[1])
    print(WarnMsg)
    log_str = c(log_str, WarnMsg)

    # reorder the columns
    cols = c(add_fields_clid[1:2], drugname_identifier, 'Concentration',
        intersect(paste0(c(paste0(add_fields_clid,'_'),'Concentration_'),
                        sort(c(2:10,2:10))),colnames(df_raw_data)))
    cols = c(cols, setdiff(colnames(df_metadata), c('Barcode', 'Template', "WellRow", "WellColumn",
                        cols, drug_identifier, paste0(drug_identifier, '_', 2:10))))
    cols = c(cols, "ReadoutValue", "BackgroundValue")
    cols = c(cols, setdiff(colnames(df_raw_data), cols))
    df_raw_data = df_raw_data[, cols ]
    df_raw_data = df_raw_data[order(df_raw_data$CellLineName, df_raw_data[,drugname_identifier],
                            df_raw_data$Concentration, df_raw_data[,duration_identifier]), ]

    return(df_raw_data)
}



#' @export
normalize_data = function(df_raw_data, log_str, selected_keys = NULL, key_values = NULL) {
    # average technical replicates and assign the right controls to each treated well

    # remove unused columns but keep barcodes to normalize by plate
    df_normalized = subset(df_raw_data, select=-c(Template, WellRow, WellColumn))

    # Identify keys for assigning the controls
    Keys = identify_keys(df_normalized)
    if (!is.null(selected_keys)) {
        Keys[names(selected_keys)] = selected_keys[names(selected_keys)]
    }

    df_normalized$CorrectedReadout = pmax(df_normalized$ReadoutValue -
                                            df_normalized$BackgroundValue,1)

    # enforced key values for end points (override selected_keys)
    Keys$Endpoint = setdiff(Keys$Endpoint, names(key_values))
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
    df_end_untrt = df_normalized[df_normalized[,duration_identifier]>0 & endpoint_value_filter &
        apply(df_normalized[,agrep('Concentration', colnames(df_normalized)),drop=F]==0,1,all),]
    df_end_mean = aggregate(df_end_untrt[,'CorrectedReadout'],
                    by = as.list(df_end_untrt[,Keys$Endpoint]), function(x) mean(x, trim= .25))
    colnames(df_end_mean)[dim(df_end_mean)[2]] = 'UntrtReadout'


    # get the untreated controls at Day 0 and perform interquartile mean
    df_day0 = df_normalized[df_normalized[,duration_identifier]==0 &
        apply(df_normalized[,agrep('Concentration', colnames(df_normalized)),drop=F]==0,1,all), ]
    df_day0_mean = aggregate(df_day0[,'CorrectedReadout'],
                by = as.list(df_day0[,Keys$Day0]), function(x) mean(x, trim= .25))
    colnames(df_day0_mean)[dim(df_day0_mean)[2]] = 'Day0Readout'

    df_controls = merge(df_end_mean, df_day0_mean[, setdiff(colnames(df_day0_mean),
                    c(duration_identifier, 'Barcode'))], all.x = T)
    if (length(setdiff(Keys$Endpoint, Keys$Day0))>0) {
        WarnMsg = paste('Not all control conditions found on the day 0 plate,',
            'dispatching values for field: ',
            paste(setdiff(Keys$Endpoint, Keys$Day0), collapse = ' ; '))
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
                                                c(duration_identifier, 'Barcode'))], 1,
                function(x) df_controls[i, setdiff(Keys$Day0, c(duration_identifier, 'Barcode')),
                                drop=F] == x))
            colnames(matches) = setdiff(Keys$Day0, c(duration_identifier, 'Barcode'))
            # try to find a good match for the day 0 (enforce same cell line)
            idx = rowSums(matches) * matches[, cellline_identifier]
            if (all(idx==0)) {next}
            match_idx = which.max(idx)
            mismatch = df_day0_mean[match_idx, setdiff(Keys$Day0,
                                                c(duration_identifier, 'Barcode'))] !=
                        df_controls[i, setdiff(Keys$Day0, c(duration_identifier, 'Barcode'))]
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


    df_normalized = merge(df_normalized[df_normalized[,duration_identifier]>0 &
        (apply(df_normalized[,agrep('Concentration', colnames(df_normalized)),drop=F]!=0,1,any) |
                !endpoint_value_filter),],
                 df_controls)

    df_normalized$RelativeViability = round(df_normalized$CorrectedReadout/
                                            df_normalized$UntrtReadout,4)
    df_normalized$GRvalue = round(2 ** (
            log2(df_normalized$CorrectedReadout / df_normalized$Day0Readout) /
                log2(df_normalized$UntrtReadout / df_normalized$Day0Readout) ), 4) - 1

    df_normalized$DivisionTime = round( df_normalized[,duration_identifier] /
                    log2(df_normalized$UntrtReadout / df_normalized$Day0Readout) , 4)


    if (any(is.na(df_normalized$Day0Readout))) {
        # need to use the reference doubling Time if day 0 missing
        InferedIdx = is.na(df_normalized$Day0Readout)
        filtered = df_normalized$ReferenceDivisionTime > (df_normalized[,duration_identifier]*2) |
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
                (df_normalized[,duration_identifier][InferedIdx] /
                        df_normalized$ReferenceDivisionTime[InferedIdx]) ) ),4) - 1
    }

    df_normalized = cbind(df_normalized[, 1:(which(colnames(df_normalized)=='ReadoutValue')-1)],
        df_normalized[, c('GRvalue', 'RelativeViability', 'DivisionTime')],
        df_normalized[, which(colnames(df_normalized)=='ReadoutValue'):(dim(df_normalized)[2]-3)])
    print('df normalized:')
    print(head(df_normalized))
    return(df_normalized)
}




#' @export
average_replicates = function(df_normalized, TrtKeys = NULL) {
    if (is.null(TrtKeys)) { TrtKeys = identify_keys(df_normalized)$Trt }

    df_averaged = aggregate(df_normalized[, c('GRvalue', 'RelativeViability', "CorrectedReadout",
                    "UntrtReadout", "Day0Readout", "DivisionTime", "ReferenceDivisionTime")],
                    by = as.list(df_normalized[,TrtKeys]), FUN = function(x) mean(x, rm.na=T))
    print('df averaged:')
    print(head(df_averaged))
    return(df_averaged)
}


#' @importFrom dplyr arrange_at group_by_at left_join summarise
#' @export
calculate_DRmetrics <-
  function(df_averaged,
           DoseRespKeys = NULL,
           studyGRvalueThresh = 4) {

   df_a = df_averaged
   colnames(df_a)[ colnames(df_a) == drugname_identifier ] = 'DrugName'

   if (is.null(DoseRespKeys)) {
       DoseRespKeys = identify_keys(df_a)$DoseResp
   } else {
       DoseRespKeys [ DoseRespKeys == drugname_identifier ] = 'DrugName'
   }
   DoseRespKeys = setdiff(DoseRespKeys, 'Concentration')
   DoseRespKeys = c(DoseRespKeys, 'DivisionTime')
   DoseRespKeys = intersect(DoseRespKeys, colnames(df_a))

   df_a$log10Concentration = log10(df_a$Concentration)

   metrics = names(ICGRlogisticFit(c(-7, -6, -5, -4), c(1, .9, .8, .7), c(1, .9, .8, .7)))
    # dummy call to get variable names

    #define set of key for merging control and study data
    mergeKeys <- setdiff(DoseRespKeys, c(drug_identifier, 'DrugName'))

    #get avereage GRvalue ('GR_0') for control data
    controlSets <-
      df_a %>% filter(DrugName %in% untreated_tag) %>%
      dplyr::group_by_at(mergeKeys) %>%
        dplyr::summarise(GR_0 = mean(GRvalue), e_0 = mean(RelativeViability))

    #get study data
    studySets <-
      df_a %>% filter(!DrugName %in% untreated_tag)

    #join study and control data
    # i.e. get  reference (average control) GRvalue ('GR_0') for study data
    fSets <-
      dplyr::left_join(studySets, controlSets, by = mergeKeys)
    # for study sets with no reference GRvalue, assing GRValue0 to 1
    fSets[is.na(fSets$GR_0), "GR_0"] <- 1
    fSets[is.na(fSets$e_0), "e_0"] <- 1

    #group study data by 'DoseRespKeys'
    gSets <- fSets %>% dplyr::group_by_at(DoseRespKeys) %>% group_split()

    print(paste(
      'Metadata variables for dose response curves:',
      paste(setdiff(
        DoseRespKeys, c(drug_identifier, cellline_identifier, paste(drug_identifier,'_', 2:10))
      ),
      collapse = ' '),
      '(',
      length(gSets),
      'groups )'
    ))

    ### ################
    ### TODO:
    ### Need to implement the metric calculation for IC50/AAC
    ### ################

    #iterate over study groups
    resL <- lapply(1:length(gSets), function(x) {
      # the 'DoseRespKeys' columns in given grup are identical for each entry
      # let's get the first record then
      repCols <- as.vector(gSets[[x]][1, DoseRespKeys])
      #get selected columns ('metrics') from GRlogisticFit output
      # (if at least 4 records with non-NA GRvalue)
      if (sum(!is.na(gSets[[x]]$GRvalue)) >= studyGRvalueThresh) {
        grLogCols <-
          ICGRlogisticFit(gSets[[x]]$log10Concentration,
                        gSets[[x]]$RelativeViability,
                        gSets[[x]]$GRvalue,
                        e_0 = gSets[[x]]$e_0[1],
                        GR_0 = gSets[[x]]$GR_0[1])[metrics]
      } else {
        grLogCols <- rep(NA, length(metrics))
      }
      cbind(repCols, t(grLogCols))
    })
    #return final data.frame
    resDf <- do.call(rbind, resL)
    resDf = resDf %>% dplyr::arrange_at(DoseRespKeys)
    colnames(resDf)[ colnames(resDf) == 'DrugName'] = drugname_identifier
    return(resDf)
  }

#' @export
identify_keys = function(df) {

    keys = list(Trt = setdiff(colnames(df), "Barcode"),
            DoseResp = setdiff(colnames(df),  'Barcode'),
            Endpoint = colnames(df)[ c(-agrep('Concentration', colnames(df)),
                                            -agrep(drug_identifier, colnames(df)),
                                            -agrep(drugname_identifier, colnames(df)))])
    keys[['Day0']] = keys[['Endpoint']]
    keys = lapply(keys, function(x) setdiff(x, c(raw_data_headers,
        normalized_results_headers, "Template", "WellRow", "WellColumn", averaged_results_headers,
            metrics_results_headers, "ReferenceDivisionTime"
    )))

    # check if all values of a key is NA
    for (k in keys[['Endpoint']]) {
        if (all(is.na(df[,k]))) {keys = lapply(keys, function(x) setdiff(x, k))}
        if (all(is.na(df[df[,duration_identifier]==0,k]))) {keys[['Day0']] = setdiff(keys[['Day0']], k)}
    }
    return(keys)
}


#' @export
#' @importFrom gneDB annotateCLIDs
#' @importFrom gCellGenomics getDrugs
cleanup_metadata = function(df_metadata, log_str) {
    log_str = c(log_str, '    cleanup_metadata')
    # clean up numberic fields
    df_metadata[,duration_identifier] = round(as.numeric(df_metadata[,duration_identifier]),6)
    # identify potential numeric fields and replace NA by 0 - convert strings in factors
    for (c in setdiff(1:dim(df_metadata)[2], c( agrep(drug_identifier, colnames(df_metadata)),
            agrep('Concentration', colnames(df_metadata)),
            grep(paste(c(cellline_identifier, manifest_headers, 'WellColumn|WellRow'),
                collapse = '|'), colnames(df_metadata)) ))) {
        vals = unique(df_metadata[,c])
        # if (is.numeric(vals)) {           # removed to ensure that missing annotation is kept NA
        #     df_metadata[is.na(df_metadata[,c]),c] = 0
        # } else
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
    gCLs = gneDB::annotateCLIDs(unique(df_metadata[,cellline_identifier]))[,
        c(DB_cell_identifier, 'celllinename', 'primarytissue', 'doublingtime')]

    colnames(gCLs) = c(cellline_identifier, add_fields_clid)
    CLIDs = unique(df_metadata[,cellline_identifier])
    bad_CL = !(CLIDs %in% gCLs[,cellline_identifier])
    if (any(bad_CL)) {
        ErrorMsg = paste('Cell line ID ', paste(CLIDs[bad_CL], collapse = ' ; '),
            ' not found in gCSI database')
        log_str = c(log_str, 'Error in cleanup_metadata:')
        log_str = c(log_str, ErrorMsg)

        stop(ErrorMsg)
        }
    print('merge with gCells')
    df_metadata = merge(df_metadata, gCLs, by.x = cellline_identifier,
                by.y = DB_cell_identifier, all.x = T)
    print(dim(df_metadata))

    # check that Gnumber_* are in the format 'G####' and add common name (or Vehicle or Untreated)

    for (i in agrep(drug_identifier, colnames(df_metadata))) { # correct case issues
        for (w in untreated_tag) {
            df_metadata[grep(w, df_metadata[,i], ignore.case = T),i] = w
        }
    }
    # -----------------------

    gDrugs = gCellGenomics::getDrugs()[,c(DB_drug_identifier, 'gcsi_drug_name')]

    # -----------------------
    gDrugs[,1] = substr(gDrugs[,1], 1, 9) # remove batch number from DB_drug_identifier
    colnames(gDrugs) = c('drug', 'DrugName')
    gDrugs = rbind(data.frame(drug=untreated_tag, DrugName=untreated_tag), gDrugs)
    gDrugs = unique(gDrugs)
    Gnbrs = unique(unlist(df_metadata[,agrep(drug_identifier, colnames(df_metadata))]))
    bad_Gn = !(Gnbrs %in% gDrugs$drug) & !is.na(Gnbrs)
    if (any(bad_Gn)) {
        # G number, but not registered
        ok_Gn = attr(regexpr('^G\\d*',Gnbrs), 'match.length')==9
        if (any(ok_Gn)) {
            WarnMsg = paste('Drug ', paste(Gnbrs[ok_Gn & bad_Gn], collapse = ' ; '),
                ' not found in gCSI database; use G# as DrugName')
            gDrugs = rbind(gDrugs, data.frame(drug=Gnbrs[ok_Gn & bad_Gn],
                    DrugName=Gnbrs[ok_Gn & bad_Gn]))
            log_str = c(log_str, 'Warning in cleanup_metadata:')
            log_str = c(log_str, WarnMsg)
            warning(WarnMsg)
        } else {
            ErrorMsg = paste('Drug ', paste(Gnbrs[!ok_Gn], collapse = ' ; '),
                ' not found in gCSI database')
            log_str = c(log_str, 'Error in cleanup_metadata:')
            log_str = c(log_str, ErrorMsg)

            stop(ErrorMsg)
        }
    }
    colnames(gDrugs)[2] = drugname_identifier
    print('merge with gDrugs for Drug 1')
    df_metadata = merge(df_metadata, gDrugs, by.x=drug_identifier, by.y='drug', all.x = T)
    print(dim(df_metadata))
    # add info for columns Gnumber_*
    for (i in grep(paste0(drug_identifier,'_\\d'), colnames(df_metadata))) {
        df_metadata[ is.na(df_metadata[,i]), i] = untreated_tag[1] # set missing values to Untreated
        gDrugs_ = gDrugs
        colnames(gDrugs_)[2] = paste0(colnames(gDrugs_)[2], substr(colnames(df_metadata)[i], 8, 12))
        print(paste('merge with gDrugs for ', i))
        df_metadata = merge(df_metadata, gDrugs_, by.x=i, by.y='drug', all.x = T)
        print(dim(df_metadata))
    }
    df_metadata[, colnames(df_metadata)[grepl(drugname_identifier, colnames(df_metadata))] ] =
        droplevels(df_metadata[,colnames(df_metadata)[grepl(drugname_identifier, colnames(df_metadata))]])

    # clean up concentration fields
    for (i in agrep('Concentration', colnames(df_metadata))) {
        trt_n = ifelse(regexpr('_\\d', colnames(df_metadata)[i])>0,
                            substr(colnames(df_metadata)[i], 15,20), 1)
        DrugID_col = ifelse(trt_n == 1, drug_identifier, paste0(drug_identifier, '_', trt_n))
        df_metadata[df_metadata[,DrugID_col] %in% untreated_tag,i] = 0 # set all untreated to 0

        DrugID_0 = setdiff(unique(df_metadata[ df_metadata[,i] == 0, DrugID_col]), untreated_tag)
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
