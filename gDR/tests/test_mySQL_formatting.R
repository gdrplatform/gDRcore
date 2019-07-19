source("../R/format_mySQL.R") # to get the function identify_keys


# DB structure is described in https://drive.google.com/open?id=1gX5ja_dSdygr2KYTqYUiENKWxu9HkOEz

# mimic the mySQL db
gdr_projects = data.frame()
condition_metadata = data.frame()
condition_codrug = data.frame()
condition_additional_treatment = data.frame()
treatment_metadata = data.frame()
response_metrics = data.frame()
response_mean = data.frame()

# run through a bunch of projects
for (project in 6) {
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

    } else if (project == 5) {

        project_processed_files = data.frame(project_number = project,
                    file_type = c('normalized', 'averaged', 'metrics'),
                    file_uri = c('test_combo1/df_normalized.tsv',
                                'test_combo1/df_averaged.tsv',
                                'test_combo1/df_metrics.tsv'))

        project_data = c(project_number = project,
                            username = 'marc',
                            description = 'Combo1',
                            approved = T,
                            date_processed = date())

    } else if (project == 6) {

        project_processed_files = data.frame(project_number = project,
                    file_type = c('normalized', 'averaged', 'metrics'),
                    file_uri = c('test_combo_3d/df_normalized.tsv',
                                'test_combo_3d/df_averaged.tsv',
                                'test_combo_3d/df_metrics.tsv'))

        project_data = c(project_number = project,
                            username = 'marc',
                            description = 'Combo_3drugs',
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

    # format the tables into the mySQL architecture
    sub_tables = format_mySQL(df_averaged, df_metrics, project,
                dim(condition_metadata)[1], dim(treatment_metadata)[1])

    # add the tables to the mySQL database
    condition_codrug = rbind(condition_codrug, sub_tables$sub_condition_codrug)

    condition_additional_treatment = rbind(condition_additional_treatment,
                                        sub_tables$sub_condition_additional_treatment)

    condition_metadata = rbind(condition_metadata, sub_tables$sub_condition_metadata)

    treatment_metadata = rbind(treatment_metadata, sub_tables$sub_treatment_metadata)

    response_metrics = rbind(response_metrics, sub_tables$sub_response_metrics)

    response_mean = rbind(response_mean, sub_tables$sub_response_mean)
    # DONE
}

print('-------')
print('whole DB')
print(dim(treatment_metadata))
print(dim(response_mean))


for (project in unique(condition_metadata$project_number)) {
    print('----')
    print(project)

    # test unpacking (will be needed for plotting)
    response_tables = extract_mySQL(project, condition_metadata, condition_additional_treatment,
                    condition_codrug, treatment_metadata, response_mean, response_metrics)

    df_averaged_re = response_tables$df_averaged
    df_metrics_re = response_tables$df_metrics

    print(dim(df_averaged_re))
    print(dim(df_metrics_re))
}
