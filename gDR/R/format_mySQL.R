source("../R/analyze_data.R") # to get the function identify_keys

# DB structure is described in https://drive.google.com/open?id=1gX5ja_dSdygr2KYTqYUiENKWxu9HkOEz

#########
# when adding to the database, needs to check for already existing co-treatments and
#   split the clid and Gnumber tables
#######


#' @export
format_mySQL = function(df_averaged, df_metrics, project_id, condition_metadata_table_length,
                            treatment_metadata_table_length) {
                        # length could be replaced by connection to mySQL database

    keys = c(setdiff(identify_keys(df_averaged)$DoseResp, 'Concentration'), 'DivisionTime')

    all_response_metadata = unique(rbind(df_averaged[, keys], df_metrics[,keys]))
    all_response_metadata$project_id = project_id

    # define the conditions (all metadata but primary Gnumber)
    condition_keys = setdiff(keys, c(drugname_identifier, 'Concentration', drug_identifier,
                                    add_fields_clid,
                                    keys[grep(drugname_identifier, keys)], 'project_id'))
    condition_all_metadata = unique(all_response_metadata[, condition_keys])

    condition_all_metadata$condition_id = condition_metadata_table_length +
                                                (1:dim(condition_all_metadata)[1])

    all_response_metadata = merge(all_response_metadata, condition_all_metadata,
                        by = condition_keys)

    # get the secondary Gnumbers in their own table if any
    N_add_Drugs = condition_keys[grep(drug_identifier, condition_keys)]
    sub_condition_codrug = data.frame()
    for (d in N_add_Drugs) {
        codrug = cbind(condition_all_metadata[, c('condition_id', d,
                gsub(drug_identifier,'Concentration', d))], gsub(paste0(drug_identifier,'_'),'', d))
        colnames(codrug)[-1] = c(drug_identifier,'Concentration','Ordinality')
        sub_condition_codrug = rbind(sub_condition_codrug, codrug)
    }

    # get the other treatment conditions if any
    add_treatments = setdiff(condition_keys[ c(-grep(drug_identifier, condition_keys),
                                        -grep('Concentration', condition_keys))],
                                    c(duration_identifier, cellline_identifier, "project_id", "condition_id"))
    sub_condition_additional_treatment = data.frame()
    for (d in add_treatments) {
        cotrt = cbind(condition_all_metadata[, c('condition_id', d)], metadata_field = d)
        colnames(cotrt)[2] = 'metadata_value'
        cotrt$metadata_value = as.character(cotrt$metadata_value)
        sub_condition_additional_treatment = rbind(sub_condition_additional_treatment, cotrt)
    }

    # get the properties for a given condition
    sub_condition_metadata = condition_all_metadata[,c(cellline_identifier, duration_identifier, 'condition_id')]
    # TODO: add division time (DivisionTime)
    sub_condition_metadata = merge(sub_condition_metadata, unique(all_response_metadata[,
                            c('project_id', 'condition_id')]), by='condition_id')

    # get the different treatments for a given condition
    sub_treatment_metadata = all_response_metadata[, c(setdiff(keys,
        c('project_id', 'CellLineName', condition_keys, 'Tissue',
                keys[grep(drugname_identifier, keys)])),
                                    'condition_id')]
    sub_treatment_metadata$treatment_id = treatment_metadata_table_length +
                                                (1:dim(sub_treatment_metadata)[1])
    all_response_metadata = merge(all_response_metadata, sub_treatment_metadata)

    # annotate the data tables with the treatment indices
    df_ = merge(df_metrics, all_response_metadata[, c(keys, 'treatment_id')])
    sub_response_metrics = df_[, setdiff(colnames(df_), keys)]

    df_ = merge(df_averaged, all_response_metadata[, c(keys, 'treatment_id')])
    sub_response_mean =
        df_[, setdiff(colnames(df_), c(keys, 'DivisionTime', 'ReferenceDivisionTime'))]

    # pack everything to be added to the mySQL database
    return(list(sub_condition_metadata = sub_condition_metadata,
            sub_condition_additional_treatment = sub_condition_additional_treatment,
            sub_condition_codrug = sub_condition_codrug,
            sub_treatment_metadata = sub_treatment_metadata,
            sub_response_mean = sub_response_mean,
            sub_response_metrics = sub_response_metrics))
}


#########
# when fetching to the database, needs to reconstruct co-treatments and
#   metadata from the clid and Gnumber tables
#######

#' @export
extract_mySQL = function(project_id, condition_metadata, condition_additional_treatment,
                condition_codrug, treatment_metadata, response_mean, response_metrics) {
                    # inputs to be replaced by connection to the mySQL database

    # TODO: properly handle the subtables for co-treatments, Gnumber, and clid

    # get all conditions for a given project
    df_metadata = merge(condition_metadata[condition_metadata$project_id == project_id,],
                            treatment_metadata, by = 'condition_id')

    # reconstruct additional_treatments
    df_additional_treatment =
                    condition_additional_treatment[condition_additional_treatment$condition_id
                                            %in% df_metadata$condition_id,]
    if (dim(df_additional_treatment)[1]>0) {
        mx_additional_treatment = unique(df_additional_treatment[,'condition_id', drop=F])
        for (meta_f in as.character(unique(df_additional_treatment$metadata_field))) {
            meta_df = df_additional_treatment[df_additional_treatment$metadata_field == meta_f,
                                                c('condition_id', 'metadata_value')]
            colnames(meta_df)[2] = meta_f
            mx_additional_treatment = merge(mx_additional_treatment, meta_df)
        }
        df_metadata = merge(df_metadata, mx_additional_treatment, by = 'condition_id', all.x=T)
    }

    # reconstruct co-druging
    df_codrug = condition_codrug[condition_codrug$condition_id
                                            %in% df_metadata$condition_id,]
    if (dim(df_codrug)[1]>0) {

        mx_codrug = unique(df_codrug[,'condition_id', drop=F])
        N_codrug = unique(df_codrug$Ordinality)
        for (codrug in N_codrug) {
            codrug_df = df_codrug[df_codrug$Ordinality == codrug,
                                    c('condition_id', drug_identifier, 'Concentration')]
            colnames(codrug_df)[2:3] = paste0(c(paste0(drug_identifier,'_'), 'Concentration_'), codrug)
            mx_codrug = merge(mx_codrug, codrug_df)
        }
        df_metadata = merge(df_metadata, mx_codrug, by = 'condition_id', all.x=T)
    }

    # pack and return the result tables
    return( list(
        df_averaged = merge(df_metadata, response_mean, by='treatment_id'),
        df_metrics = merge(df_metadata, response_metrics, by='treatment_id')
    ))

}
