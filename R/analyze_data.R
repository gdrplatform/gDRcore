library(gCellGenomics) # best reference for cell line and drug names?

library(devtools)
# install_git('https://stash.intranet.roche.com/stash/scm/~hafnerm6/gcsiutils.git',
    # ref = 'GRimplementation')
load_all('~/workspace/Rpackages/gcsiutils/') # local copy of gcsiutils/GRimplementation

Overall_function = function(manifest_file, template_file, results_file, output_files) {
    # output_files should contain file names for :
    #   log_file, QC_file, raw_result, process_results, metrics_results
    log_file <- file(output_files['log_file'], open = "wt")

    df_raw_data = load_merge_data(manifest_file, template_file, results_file, log_file)

    # output_QC_byPlate(df_raw_data, output_files['QC_file']) # TODO: check column/row bias

    Keys = identify_keys(df_raw_data) # may be manually changed

    df_normalized = normalize_data(df_raw_data, log_file, Keys)
    df_averaged = average_replicates(df_normalized, Keys$Trt)

    df_metrics = calculate_DRmetrics(df_averaged, Keys$DoseResp)

    return(list(raw=df_raw_data,
            normalized=df_normalized,
            averaged=df_averaged,
            metrics=df_metrics))
}

load_merge_data = function(manifest_file, template_file, results_file, log_file) {

    manifest = load_manifest(manifest_file, log_file)
    treatments = load_templates(template_file, log_file)
    data = load_results(results_file, log_file)

    # check the all template files are available
    if (!all(unique(manifest$Template[manifest$Barcode %in% data$Barcode])
                    %in% basename(template_file))) {
        ErrorMsg = paste('Some template files are missing:',
                    paste(setdiff(unique(manifest$Template[manifest$Barcode %in% data$Barcode]),
                                     basename(template_file)), collapse = ' ; '))
        writeLines('Error in load_merge_data:', log_file)
        writeLines(ErrorMsg, log_file)
        close(log_file)
        stop(ErrorMsg)
    }


    # merge manifest and treatment files first
    df_metadata = merge(manifest, treatments, by = 'Template')
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
                writeLines('Warning in load_merge_data:', log_file)
                writeLines(WarnMsg, log_file)
                warning(WarnMsg)
            }
        df_metadata[,paste0(m_col,'.x')] = NULL
        df_metadata[,paste0(m_col,'.y')] = NULL
    }

    # check for the expected columns
    expected_headers = c('CLid')
    headersOK = expected_headers %in% colnames(df_metadata)
    if (any(!headersOK)) {
        ErrorMsg = paste('df_metadata',
            'does not contains all expected headers: ',
            paste(expected_headers[ !(expected_headers %in% col_df) ], collpase = ' ; '),
            ' required')
        writeLines('Error in load_merge_data:', log_file)
        writeLines(ErrorMsg, log_file)
        close(log_file)
        stop(ErrorMsg)
    }

    # clean up the metadata
    cleanedup_metadata = cleanup_metadata(df_metadata)
    stopifnot( dim(cleanedup_metadata)[1] == dim(df_metadata)[1] ) # should not happen

    df_merged = merge(cleanedup_metadata, data, by = c('Barcode', 'WellRow', 'WellColumn'))
    if (dim(df_merged)[1] != dim(data)[1]) {# need to identify issue and output relevant warning
        WarnMsg = 'Not all results have been matched with treatments; merged table is smaller than data table'
        writeLines('Warning in load_merge_data:', log_file)
        writeLines(WarnMsg, log_file)
        warning(WarnMsg)
    }
    if (dim(df_merged)[1] != dim(df_metadata)[1]) {# need to identify issue and output relevant warning
        WarnMsg = 'Not all treatments have been matched with results; merged table is smaller than metadata table'
        writeLines('Warning in load_merge_data:', log_file)
        writeLines(WarnMsg, log_file)
        warning(WarnMsg)
    }

    # remove wells not labeled
    df_raw_data = df_merged[ !is.na(df_merged$Gnumber), ]
    WarnMsg = sprintf('%i well loaded, %i discarded for lack of annotation, %i data point selected',
                dim(data)[1], sum(is.na(df_merged$Gnumber)), dim(df_raw_data)[1])
    print(WarnMsg)
    writeLines(WarnMsg, log_file)

    # reorder the columns
    cols = c('CellLineName', 'Tissue', 'DrugName', 'Concentration',
        intersect(paste0(c('DrugName_','Concentration_'),sort(c(2:10,2:10))),colnames(df_raw_data)))
    cols = c(cols, setdiff(colnames(df_metadata), c('Barcode', 'Template', "WellRow", "WellColumn",
                        cols, 'Gnumber', paste0('Gnumber_', 2:10))))
    cols = c(cols, "ReadoutValue", "BackgroundValue")
    cols = c(cols, setdiff(colnames(df_raw_data), cols))
    df_raw_data = df_raw_data[, cols ]
    df_raw_data = df_raw_data[order(df_raw_data$CellLineName, df_raw_data$DrugName,
                            df_raw_data$Concentration, df_raw_data$Time), ]

    return(df_raw_data)
}






normalize_data = function(df_raw_data, log_file, selected_keys = NULL) {
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

    # get the untreated controls at endpoint and perform interquartile mean
    df_end_untrt = df_normalized[df_normalized$Time>0 &
        apply(df_normalized[,agrep('Concentration', colnames(df_normalized)),drop=F]==0,1,all),]
    df_end_mean = aggregate(df_end_untrt[,'CorrectedReadout'],
                    by = as.list(df_end_untrt[,Keys$Endpoint]), function(x) mean(x, trim= .25))
    colnames(df_end_mean)[dim(df_end_mean)[2]] = 'UntrtReadout'


    # get the untreated controls at Day 0 and perform interquartile mean
    df_day0 = df_normalized[df_normalized$Time==0 &
        apply(df_normalized[,agrep('Concentration', colnames(df_normalized)),drop=F]==0,1,all), ]
    df_day0_mean = aggregate(df_day0[,'CorrectedReadout'],
                by = as.list(df_day0[,Keys$Day0]), function(x) mean(x, trim= .25))
    colnames(df_day0_mean)[dim(df_day0_mean)[2]] = 'Day0Readout'

    df_controls = merge(df_end_mean, subset(df_day0_mean, select=-c(Time, Barcode)), all.x = T)
    if (length(setdiff(Keys$Endpoint, Keys$Day0))>0) {
        WarnMsg = paste('Not all control conditions found on the day 0 plate,',
            'dispatching values for field: ',
            paste(setdiff(Keys$Endpoint, Keys$Day0), collapse = ' ; '))
        writeLines('Warning in normalize_data:', log_file)
        writeLines(WarnMsg, log_file)
        warning(WarnMsg)
    }
    # identify missing values in the Day0 that needs to be matched (usually for co-treatments)
    df_controls_NA = which(is.na(df_controls$Day0Readout))

    if (length(df_controls_NA)>0) {
        dispatched = NULL
        for (i in df_controls_NA){
            matches = t(apply(df_day0_mean[, setdiff(Keys$Day0, c('Time', 'Barcode'))], 1,
                function(x) df_controls[i, setdiff(Keys$Day0, c('Time', 'Barcode')),
                                drop=F] == x))
            colnames(matches) = setdiff(Keys$Day0, c('Time', 'Barcode'))
            # try to find a good match for the day 0 (enforce same cell line)
            idx = rowSums(matches) * matches[,'CLid']
            if (all(idx==0)) {next}
            match_idx = which.max(idx)
            mismatch = df_day0_mean[match_idx, setdiff(Keys$Day0, c('Time', 'Barcode'))] !=
                        df_controls[i, setdiff(Keys$Day0, c('Time', 'Barcode'))]
            dispatched = c(dispatched, colnames(mismatch)[mismatch])
            df_controls[i, 'Day0Readout'] = df_day0_mean[match_idx, 'Day0Readout']
        }
        WarnMsg = paste('Not all control conditions found on the day 0 plate,')
        WarnMsg = ifelse(length(dispatched)>0,
                    paste(WarnMsg,'dispatching values for mismatches in field: ',
                    paste(unique(dispatched), collapse = ' ; ')),
                    paste(WarnMsg,'some Day0 are not being matched'))
        writeLines('Warning in normalize_data:', log_file)
        writeLines(WarnMsg, log_file)
        warning(WarnMsg)
    }


    df_normalized = merge(df_normalized[df_normalized$Time>0 &
        apply(df_normalized[,agrep('Concentration', colnames(df_normalized)),drop=F]!=0,1,any) ,],
                 df_controls)

    df_normalized$RelViability = round(df_normalized$CorrectedReadout/df_normalized$UntrtReadout,4)
    df_normalized$GRvalue = round(2 ** (
            log2(df_normalized$CorrectedReadout / df_normalized$Day0Readout) /
                log2(df_normalized$UntrtReadout / df_normalized$Day0Readout) ), 4) - 1

    df_normalized$DivisionTime = round( df_normalized$Time /
                    log2(df_normalized$UntrtReadout / df_normalized$Day0Readout) , 4)


    if (any(is.na(df_normalized$Day0Readout))) {
        # need to use the reference doubling Time if day 0 missing
        InferedIdx = is.na(df_normalized$Day0Readout)
        filtered = df_normalized$ReferenceDivisionTime > (df_normalized$Time*2) |
            is.na(df_normalized$ReferenceDivisionTime)
        WarnMsg = paste('Missing day 0 information --> calculate GR value based on reference',
            'doubling time')
        WarnMsg = ifelse(!any(filtered & InferedIdx), WarnMsg,
            paste(WarnMsg, '; filtering', sum(filtered & InferedIdx),
                'conditions because of too short assay:',
                paste(unique(df_normalized$CellLineName[filtered & InferedIdx]), collpase=' ; ')))
        writeLines('Warning in normalize_data:', log_file)
        writeLines(WarnMsg, log_file)
        warning(WarnMsg)
        InferedIdx = !filtered & InferedIdx
        # calculate GR values using formula from https://www.nature.com/articles/nbt.3882
        df_normalized$GRvalue[InferedIdx] = round(2 ^ ( 1 + (log2(pmin(1.25,
                df_normalized[InferedIdx,'RelViability'])) /
                (df_normalized$Time[InferedIdx] /
                        df_normalized$ReferenceDivisionTime[InferedIdx]) ) ),4) - 1
    }

    df_normalized = cbind(df_normalized[, 1:(which(colnames(df_normalized)=='ReadoutValue')-1)],
        df_normalized[, c('GRvalue', 'RelViability', 'DivisionTime')],
        df_normalized[, which(colnames(df_normalized)=='ReadoutValue'):(dim(df_normalized)[2]-3)])
}





average_replicates = function(df_normalized, TrtKeys = NULL) {
    if (is.null(TrtKeys)) { TrtKeys = identify_keys(df_normalized)$Trt }

    df_averaged = aggregate(df_normalized[, c('GRvalue', 'RelViability', "CorrectedReadout",
                    "UntrtReadout", "Day0Readout", "DivisionTime", "ReferenceDivisionTime")],
                    by = as.list(df_normalized[,TrtKeys]), FUN = function(x) mean(x, rm.na=T))
    return(df_averaged)
}


calculate_DRmetrics = function(df_averaged, DoseRespKeys = NULL, force = FALSE, cap = FALSE) {
    if (is.null(DoseRespKeys)) { DoseRespKeys = identify_keys(df_averaged)$DoseResp }
    DoseRespKeys = setdiff(DoseRespKeys, 'Concentration')

    df_GR = df_averaged # may be worthwhile to use consistent variable names at some point

    # colnames(df_GR)[colnames(df_GR) %in% 'GRvalue'] = 'GR'
    # colnames(df_GR)[colnames(df_GR) %in% 'Concentration'] = 'concentration'
    # df_metrics = calculate_GRmetrics(df_GR, meta_variables = DoseRespKeys, force=force, cap=cap)

    # copied from  calculate_GRmetrics
    df_GR$log10Concentration = log10(as.numeric(as.character(df_GR$Concentration)))
    for (v in DoseRespKeys) df_GR[, v] = as.factor(df_GR[, v])
    df_metrics = unique(df_GR[, DoseRespKeys])
    print(paste('Metadata variables for dose response curves:',
                do.call(paste,as.list(DoseRespKeys)), '(',
                dim(df_metrics)[1], 'groups)'))

    metrics = names(GRlogisticFit(c(-7,-6,-5,-4),c(1,.6,.3,.3))) # dummy call to get variable names
    df_metrics = cbind(df_metrics, as.data.frame(matrix(NA, dim(df_metrics)[1], length(metrics))))
    colnames(df_metrics) = c(DoseRespKeys, metrics)
    # loop through all conditions
    for (i in 1:dim(df_metrics)[1]) {
        sub_meta = !is.na(df_metrics[i,DoseRespKeys])
        idx = apply(df_GR[,DoseRespKeys[sub_meta]], 1,
            function(x) all(x == df_metrics[i,DoseRespKeys[sub_meta]]))
        log10concs = df_GR$log10Concentration[idx]
        GRvalues = df_GR$GRvalue[idx]
        if (sum(!is.na(GRvalues))<4) next
        df_metrics[i, metrics] = GRlogisticFit(log10concs, GRvalues)
    }


    # handling cases where DrugName = Vehicle/Untrt --> these are reference for other conditions
    df_0 = df_metrics[df_metrics$DrugName %in% c('Vehicle', 'Untreated'),DoseRespKeys]
    df_metrics = df_metrics[!(df_metrics$DrugName %in% c('Vehicle', 'Untreated')),]
    df_GR[ df_GR[,DoseRespKeys] %in% df_0 ,]

    return(df_metrics)
}



identify_keys = function(df) {
    # c(paste0('Concentration_', 2:10), paste0('Gnumber_', 2:10), paste0('DrugName_', 2:10)
    keys = list(Trt = setdiff(colnames(df), "Barcode"),
            DoseResp = setdiff(colnames(df),  'Barcode'),
            Endpoint = colnames(df)[ c(-agrep('Concentration', colnames(df)),
                                            -agrep('Gnumber', colnames(df)),
                                            -agrep('DrugName', colnames(df)))])
    keys[['Day0']] = keys[['Endpoint']]
    keys = lapply(keys, function(x) setdiff(x, c("ReadoutValue", "BackgroundValue",
            "UntrtReadout", "CorrectedReadout", "Day0Readout", "Template", "WellRow", "WellColumn",
            'GRvalue', 'RelViability', 'DivisionTime', 'ReferenceDivisionTime')))
    # check if all values of a key is NA
    for (k in keys[['Endpoint']]) {
        if (all(is.na(df[,k]))) {keys = lapply(keys, function(x) setdiff(x, k))}
        if (all(is.na(df[df$Time==0,k]))) {keys[['Day0']] = setdiff(keys[['Day0']], k)}
    }
    return(keys)
}





cleanup_metadata = function(df_metadata) {

    # clean up numberic fields
    df_metadata$Time = round(as.numeric(df_metadata$Time),6)
    # identify potential numeric fields and replace NA by 0 - convert strings in factors
    for (c in setdiff(1:dim(df_metadata)[2], c( agrep('Gnumber', colnames(df_metadata)),
            agrep('Concentration', colnames(df_metadata)),
            grep('Time|CLid|Barcode|WellRow|WellColumn|Template', colnames(df_metadata)) ))) {
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
                writeLines('Warning in cleanup_metadata:', log_file)
                writeLines(WarnMsg, log_file)
                warning(WarnMsg)
            } else {
                is.na(df_metadata[,c]) = 0
                df_metadata[,c] = as.numeric(df_metadata[,c])
                WarnMsg = paste('Metadata field ', colnames(df_metadata)[c],
                                    ' converted to numeric values')
                writeLines('Warning in cleanup_metadata:', log_file)
                writeLines(WarnMsg, log_file)
                warning(WarnMsg)
            }
        }
    }

    # check that CLid are in the format 'CL####' and add common name
    gCLs = gCellGenomics::getSamples()[,c('clid', 'celllinename', 'tissue', 'doublingtime')]
    colnames(gCLs)[2:4] = c('CellLineName', 'Tissue', 'ReferenceDivisionTime')
    CLids = unique(df_metadata$CLid)
    bad_CL = !(CLids %in% gCLs$clid)
    if (any(bad_CL)) {
        ErrorMsg = paste('Cell line ID ', paste(CLids[bad_CL], collapse = ' ; '),
            ' not found in gCSI database')
        writeLines('Error in load_all_data:', log_file)
        writeLines(ErrorMsg, log_file)
        close(log_file)
        stop(ErrorMsg)
        }
    df_metadata = merge(df_metadata, gCLs, by.x='CLid', by.y='clid', all.x = T)

    # check that Gnumber_* are in the format 'G####' and add common name (or Vehicle or Untreated)
    untrt_flag = c('Vehicle', 'Untreated')
    for (i in agrep('Gnumber', colnames(df_metadata))) { # correct case issues
        for (w in untrt_flag) {
            df_metadata[grep(w, df_metadata[,i], ignore.case = T),i] = w
        }
    }
    gDrugs = gCellGenomics::getDrugs()[,c('drug', 'gcsi_drug_name')]
    gDrugs[is.na(gDrugs$gcsi_drug_name), ]
    gDrugs$drug = substr(gDrugs$drug, 1, 9)
    colnames(gDrugs)[2] = 'DrugName'
    gDrugs = rbind(data.frame(drug=untrt_flag, DrugName=untrt_flag),
                gDrugs)
    Gnbrs = unique(unlist(df_metadata[,agrep('Gnumber', colnames(df_metadata))]))
    bad_Gn = !(Gnbrs %in% gDrugs$drug) & !is.na(Gnbrs)
    if (any(bad_Gn)) {
        # G number, but not registered
        ok_Gn = attr(regexpr('^G\\d*',Gnbrs), 'match.length')==9
        if (any(ok_Gn)) {
            WarnMsg = paste('Drug ', paste(Gnbrs[ok_Gn & bad_Gn], collapse = ' ; '),
                ' not found in gCSI database; use G# as DrugName')
            gDrugs = rbind(gDrugs, data.frame(drug=Gnbrs[ok_Gn & bad_Gn],
                    DrugName=Gnbrs[ok_Gn & bad_Gn]))
            writeLines('Warning in cleanup_metadata:', log_file)
            writeLines(WarnMsg, log_file)
            warning(WarnMsg)
        } else {
            ErrorMsg = paste('Drug ', paste(Gnbrs[!ok_Gn], collapse = ' ; '),
                ' not found in gCSI database')
            writeLines('Error in cleanup_metadata:', log_file)
            writeLines(ErrorMsg, log_file)
            close(log_file)
            stop(ErrorMsg)
        }
    }
    df_metadata = merge(df_metadata, gDrugs, by.x='Gnumber', by.y='drug', all.x = T)
    # add info for columns Gnumber_*
    for (i in grep('Gnumber_\\d', colnames(df_metadata))) {
        df_metadata[ is.na(df_metadata[,i]), i] = 'Untreated' # set missing values to Untreated
        gDrugs_ = gDrugs
        colnames(gDrugs_)[2] = paste0(colnames(gDrugs_)[2], substr(colnames(df_metadata)[i], 8, 12))
        df_metadata = merge(df_metadata, gDrugs_, by.x=i, by.y='drug', all.x = T)
    }

    # clean up concentration fields
    for (i in agrep('Concentration', colnames(df_metadata))) {
        trt_n = ifelse(regexpr('_\\d', colnames(df_metadata)[i])>0,
                            substr(colnames(df_metadata)[i], 15,20), 1)
        G_col = ifelse(trt_n == 1, 'Gnumber', paste0('Gnumber_', trt_n))
        df_metadata[df_metadata[,G_col] %in% untrt_flag,i] = 0 # set all untreated to 0

        Gnum_0 = setdiff(unique(df_metadata[ df_metadata[,i] == 0, G_col]), untrt_flag)
        Gnum_0 = Gnum_0[!is.na(Gnum_0)]
        if (length(Gnum_0)>0) {
            WarnMsg = paste('Some concentration for ', G_col,
                            ' are 0: ', paste(Gnum_0, collapse = ' ; '))
            writeLines('Warning in cleanup_metadata:', log_file)
            writeLines(WarnMsg, log_file)
            warning(WarnMsg)
        }
        df_metadata[,i] = round(as.numeric(df_metadata[,i]),6) # avoid mismatch due to string truncation
    }


    return(df_metadata)
}
