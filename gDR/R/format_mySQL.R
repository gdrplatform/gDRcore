source("analyze_data.R") # to get the function identify_keys
library(reshape2)

# DB structure is described in https://drive.google.com/open?id=1gX5ja_dSdygr2KYTqYUiENKWxu9HkOEz

# mimic the mySQL db
gdr_projects = data.frame()
condition_metadata = data.frame()
condition_codrug = data.frame()
condition_additional_treatment = data.frame()
treatment_metadata = data.frame()
response_metrics = data.frame()
response_mean = data.frame()

for (project in 1:4) {
    print('----')
    print(project)
    if (project == 1) {
        project_raw_files = data.frame(project_number = project,
                    file_type = c('manifest', 'template', 'template', 'results', 'results'),
                    file_uri = c('test_Kyle1/Manifest file_0077&Palbo.tsv',
                                'test_Kyle1/Template_BT474_untreated.xlsx',
                                'test_Kyle1/Template_BT474_7daytreated_KE.xlsx',
                                'test_Kyle1/cyquant_BT474_day0.xlsx',
                                'test_Kyle1/cyquant_BT474_day7.xlsx'))

        project_processed_files = data.frame(project_number = project,
                    file_type = c('normalized', 'averaged', 'metrics'),
                    file_uri = c('test_Kyle1/ref_df_normalized.tsv',
                                'test_Kyle1/ref_df_averaged.tsv',
                                'test_Kyle1/ref_df_metrics.tsv'))

        project_data = c(project_number = project,
                            username = 'kyle',
                            description = 'Kyle1',
                            approved = NULL,
                            date_processed = date())


    } else if (project == 2) {
        project_raw_files = data.frame(project_number = project,
                    file_type = c('manifest', 'template', 'template', 'results', 'results'),
                    file_uri = c('test_Kyle2/Manifest_0077vs0068.xls',
                                'test_Kyle2/Template_0077vs0068_Untreated.xls',
                                'test_Kyle2/Template_0077_0068_7daytreated.xls',
                                'test_Kyle2/Result_0077vs0068_day7.xls',
                                'test_Kyle2/Result_0077vs0068_day0.xls'))

        project_processed_files = data.frame(project_number = project,
                    file_type = c('normalized', 'averaged', 'metrics'),
                    file_uri = c('test_Kyle2/ref_df_normalized.tsv',
                                'test_Kyle2/ref_df_averaged.tsv',
                                'test_Kyle2/ref_df_metrics.tsv'))

        project_data = c(project_number = project,
                            username = 'kyle',
                            description = 'Kyle2',
                            approved = T,
                            date_processed = date())

    } else if (project == 3) {
        project_raw_files = data.frame(project_number = project,
                    file_type = c('manifest', 'template', 'template', 'results',
                                        'results', 'results', 'results'),
                    file_uri = c('test_Wei2/Manifest file_ET&Palbo.xls',
                                'test_Wei2/Template_untreated.xls',
                                'test_Wei2/Template_7daytreated.xls',
                                'test_Wei2/Result_20181024_untreated.xls',
                                'test_Wei2/Result_20181031_untreated.xls',
                                'test_Wei2/Result_20181024_7daytreated.xls',
                                'test_Wei2/Result_20181031_7daytreated.xls'))

        project_processed_files = data.frame(project_number = project,
                    file_type = c('normalized', 'averaged', 'metrics'),
                    file_uri = c('test_Wei2/ref_df_normalized.tsv',
                                'test_Wei2/ref_df_averaged.tsv',
                                'test_Wei2/ref_df_metrics.tsv'))

        project_data = c(project_number = project,
                            username = 'wei',
                            description = 'Wei2',
                            approved = T,
                            date_processed = date())

    } else if (project == 4) {
        results_file <- list.files('../inst/testdata/test_Nick', pattern = "12\\w*", full.names=T)
        project_raw_files = data.frame(project_number = project,
                    file_type = c('manifest', 'template', 'template',
                                    array('results',length(results_file))),
                    file_uri = c('test_Nick/1000700950096_manifest_final.xlsx',
                                'test_Nick/trtmt1.xlsx', 'test_Nick/untreated.xlsx',
                                gsub("../inst/testdata/", '', results_file)))

        project_processed_files = data.frame(project_number = project,
                    file_type = c('normalized', 'averaged', 'metrics'),
                    file_uri = c('test_Nick/ref_df_normalized.tsv',
                                'test_Nick/ref_df_averaged.tsv',
                                'test_Nick/ref_df_metrics.tsv'))

        project_data = c(project_number = project,
                            username = 'nick',
                            description = 'Nick',
                            approved = T,
                            date_processed = date())

    }
    # will be replace directly by the variable, not the files
    df_normalized = read.table(paste0('../inst/testdata/',project_processed_files[1,3]),
                        header = T, sep='\t')
    df_averaged = read.table(paste0('../inst/testdata/', project_processed_files[2,3]),
                        header = T, sep='\t')
    df_metrics = read.table(paste0('../inst/testdata/', project_processed_files[3,3]),
                        header = T, sep='\t')

    print(dim(df_averaged))
    print(dim(df_metrics))

    # need to update:
    #  - CLID --> clid
    #  - add DivisionTime in the df_metrics table

    keys = setdiff(identify_keys(df_averaged)$DoseResp, 'Concentration')

    all_response_metadata = unique(rbind(df_averaged[, keys], df_metrics[,keys]))
    all_response_metadata$project_number = project

    # define the conditions (all metadata but Gnumber)
    condition_keys = setdiff(keys, c("DrugName", 'Concentration', 'Gnumber', "CellLineName", "Tissue",
                                    keys[grep('DrugName', keys)], 'project_number'))
    condition_all_metadata = unique(all_response_metadata[, condition_keys])  # , 'DivisionTime')]

    condition_all_metadata$condition_number = dim(condition_metadata)[1] +
                                                (1:dim(condition_all_metadata)[1])

    all_response_metadata = merge(all_response_metadata, condition_all_metadata, by = condition_keys)

    # get the secondary Gnumbers in their own table if any
    N_add_Drugs = condition_keys[grep('Gnumber', condition_keys)]
    sub_condition_codrug = data.frame()
    for (d in N_add_Drugs) {
        codrug = condition_all_metadata[, c('condition_number', d,
                    gsub('Gnumber','Concentration', d))]
        colnames(codrug)[-1] = c('Gnumber','Concentration')
        sub_condition_codrug = rbind(sub_condition_codrug, codrug)
    }
    condition_codrug = rbind(condition_codrug, sub_condition_codrug)

    # get the other conditions if any
    add_treatments = setdiff(condition_keys[ c(-grep('Gnumber', condition_keys),
                                        -grep('Concentration', condition_keys))],
                                    c('Time', 'CLID', "project_number", "condition_number"))
    sub_condition_additional_treatment = data.frame()
    for (d in add_treatments) {
        cotrt = cbind(condition_all_metadata[, c('condition_number', d)], metadata_field = d)
        colnames(cotrt)[2] = 'metadata_value'
        cotrt$metadata_value = as.character(cotrt$metadata_value)
        sub_condition_additional_treatment = rbind(sub_condition_additional_treatment, cotrt)
    }
    condition_additional_treatment = rbind(condition_additional_treatment,
                                        sub_condition_additional_treatment)


    sub_condition_metadata = condition_all_metadata[,c('CLID', 'Time', 'condition_number')]
    # add division time
    sub_condition_metadata = merge(sub_condition_metadata, unique(all_response_metadata[,
                            c('project_number', 'condition_number')]), by='condition_number')
    condition_metadata = rbind(condition_metadata, sub_condition_metadata)


    sub_treatment_metadata = all_response_metadata[, c(setdiff(keys,
        c('project_number', 'CellLineName', condition_keys, 'Tissue', keys[grep('DrugName', keys)])),
                                    'condition_number')]
    sub_treatment_metadata$treatment_number = dim(treatment_metadata)[1] +
                                                (1:dim(sub_treatment_metadata)[1])
    all_response_metadata = merge(all_response_metadata, sub_treatment_metadata)

    treatment_metadata = rbind(treatment_metadata, sub_treatment_metadata)

    df_ = merge(df_metrics, all_response_metadata[, c(keys, 'treatment_number')])
    response_metrics = rbind(response_metrics, df_[, setdiff(colnames(df_), keys)])

    df_ = merge(df_averaged, all_response_metadata[, c(keys, 'treatment_number')])
    response_mean = rbind(response_mean,
            df_[, setdiff(colnames(df_), c(keys, 'DivisionTime', 'ReferenceDivisionTime'))] )

}

print('-------')
print('whole DB')
print(dim(treatment_metadata))
print(dim(response_mean))

# unpack (will be needed for plotting)
for (project in 1:4) {
    print('----')
    print(project)

    df_metadata = merge(condition_metadata[condition_metadata$project_number == project,],
                            treatment_metadata, by = 'condition_number')
    # reconstruct additional_treatments
    df_additional_treatment =
                    condition_additional_treatment[condition_additional_treatment$condition_number
                                            %in% df_metadata$condition_number,]
    if (dim(df_additional_treatment)[1]>0) {
        mx_additional_treatment = unique(df_additional_treatment[,'condition_number', drop=F])
        for (meta_f in as.character(unique(df_additional_treatment$metadata_field))) {
            meta_df = df_additional_treatment[df_additional_treatment$metadata_field == meta_f,
                                                c('condition_number', 'metadata_value')]
            colnames(meta_df)[2] = meta_f
            mx_additional_treatment = merge(mx_additional_treatment, meta_df)
        }
        df_metadata = merge(df_metadata, mx_additional_treatment, by = 'condition_number', all.x=T)
    }

    # reconstruct co-druging
    df_codrug = condition_codrug[condition_codrug$condition_number
                                            %in% df_metadata$condition_number,]
    if (dim(df_codrug)[1]>0) {
        df_codrug = df_codrug[order( as.character(df_codrug$Gnumber), df_codrug$Concentration),]
        df_codrug$cumul_condition =  c(1,sapply(2:dim(df_codrug)[1], function(x)
                sum( df_codrug$condition_number[x] == df_codrug$condition_number[seq(1,x-1,1)])+1))
        mx_codrug = unique(df_codrug[,'condition_number', drop=F])
        N_codrug = max(table(mx_additional_treatment))
        for (codrug in 1:N_codrug) {
            codrug_df = df_codrug[df_codrug$cumul_condition == codrug,
                                        c('condition_number', 'Gnumber', 'Concentration')]
            colnames(codrug_df)[2:3] = paste0(c('Gnumber_', 'Concentration_'), codrug+1)
            mx_codrug = merge(mx_codrug, codrug_df)
        }
        df_metadata = merge(df_metadata, mx_codrug, by = 'condition_number', all.x=T)
    }


    df_averaged_re = merge(df_metadata, response_mean, by='treatment_number')
    df_metrics_re = merge(df_metadata, response_metrics, by='treatment_number')

    print(dim(df_averaged_re))
    print(dim(df_metrics_re))
}
