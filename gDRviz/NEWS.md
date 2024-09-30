## gDRviz 1.6.64 - 2024-09-30
* fix issue with managing additional perturbations

## gDRviz 1.6.63 - 2024-09-24
* fix issue with redirection

## gDRviz 1.6.62 - 2024-08-26
* fix an issue with combo modules/tabs staying active after data change

## gDRviz 1.6.61 - 2024-08-22
* update shiny test

## gDRviz 1.6.60 - 2024-08-12
* restore url access

## gDRviz 1.6.59 - 2024-08-06
* fix reactivity in response curve modules

## gDRviz 1.6.58 - 2024-08-06
* utilize `set_unique_identifiers` fun

## gDRviz 1.6.57 - 2024-08-05
* update call of `get_additional_variables` fun

## gDRviz 1.6.56 - 2024-08-01
* improve user experience after merging data

## gDRviz 1.6.55 - 2024-07-29
* refactor: server.R function
    * support averaged/metrics data for combination experiments
    * support edge-case with no single-agent data
    * improve logic for hiding/showing sa/combination plugins
    * switch directly to DrugCombo plugin if combo data only

## gDRviz 1.6.54 - 2024-07-19
* hide combo tab in manageData for missing combo data

## gDRviz 1.6.53 - 2024-07-17
* rename gDRsearch2 to gDRsearch in sidebar

## gDRviz 1.6.52 - 2024-07-15
* fix edge-case in gDRsearch2 logic (for combination data)

## gDRviz 1.6.51 - 2024-07-11
* make app log level customizable via env variable

## gDRviz 1.6.50 - 2024-06-20
* improve logging with logExtra

## gDRviz 1.6.49 - 2024-06-19
* move gDRsearch-specific code to gDRsearch

## gDRviz 1.6.48 - 2024-06-17
* update dependencies

## gDRviz 1.6.47 - 2024-06-04
* update logic for handling AWS S3 connections

## gDRviz 1.6.46 - 2024-05-31
* add admin mode button

## gDRviz 1.6.45 - 2024-05-22
* update namespace

## gDRviz 1.6.44 - 2024-05-17
* remove logic dependent on `filter_mae`

## gDRviz 1.6.43 - 2024-05-08
* move FrontPage, miniSummary, moduleDisplay and aceEditor modules to `gDRcomponents` package

## gDRviz 1.6.42 - 2024-05-07
* move deployment functions from `gDRcomponents` to `gDRviz`

## gDRviz 1.6.41 - 2024-04-30
* update namespace for `get_additional_variables` fun

## gDRviz 1.6.40 - 2024-04-29
* fix cedar test

## gDRviz 1.6.39 - 2024-04-11
* rename `get_project_table` into `get_project_tables`

## gDRviz 1.6.38 - 2024-03-06
* move pickDataset module to gDRsearch
* remove gDRsearch as dependencies
* update "Parameterized URLs" vignette

## gDRviz 1.6.37 - 2024-03-06
* update documentation links

## gDRviz 1.6.36 - 2024-02-26
* add support for combo from gDRsearch2

## gDRviz 1.6.35 - 2024-01-16
* switch from 'matrix' into 'combination'

## gDRviz 1.6.34 - 2023-11-24
* add popup for managing additional perturbations/metadata

## gDRviz 1.6.33 - 2023-11-21
* remove shinytest related code

## gDRviz 1.6.32 - 2023-11-17
* reorder tabs in gDRsearch and update shiny tests

## gDRviz 1.6.31 - 2023-11-15
* restore visualization regardless of the data source

## gDRviz 1.6.30 - 2023-10-25
* add support for data from gDRsearch2 in manageData

## gDRviz 1.6.29 - 2023-10-04
* add support for data from gDRsearch2 in manageData

## gDRviz 1.6.28 - 2023-09-27
* add support for triggering from gDRsearch2

## gDRviz 1.6.27 - 2023-09-14
* switch from `Normalized` to `Averaged` assay in dose-response curve

## gDRviz 1.6.26 - 2023-09-13
* fix for launching gDRviz via URL

## gDRviz 1.6.25 - 2023-08-16
* add support for gDRsearch2 URL query

## gDRviz 1.6.24 - 2023-08-10
* remove case sensitiveness for param values in URL

## gDRviz 1.6.23 - 2023-08-08
* introduce CD

## gDRviz 1.6.22 - 2023-08-07
* add functionality for "GO" button in URL

## gDRviz 1.6.21 - 2023-08-01
* add vignettes for parametrized URL

## gDRviz 1.6.20 - 2023-07-28
* improve logic for tests that use local cache

## gDRviz 1.6.19 - 2023-07-26
* switch from shinytest to shinytest2

## gDRviz 1.6.18 - 2023-07-26
* update admin user list

## gDRviz 1.6.17 - 2023-07-25
* fix for unreproducible behavior of shinytest

## gDRviz 1.6.16 - 2023-07-24
* update logic for launching gDRviz with specified dataset via URL

## gDRviz 1.6.15 - 2023-07-17
* add cache for projects list

## gDRviz 1.6.14 - 2023-07-10
* fix issue with NA values in ifelse statement

## gDRviz 1.6.13 - 2023-07-06
* enable launching gDRviz with specified dataset via URL

## gDRviz 1.6.12 - 2023-06-27
* fix bug in the displayLogs module

## gDRviz 1.6.11 - 2023-06-19
* simplify/standardize conversion of assays in gDRviz

## gDRviz 1.6.10 - 2023-06-15
* provide info about SAVE_FOLDER and LOG_FILE via `sv` list

## gDRviz 1.6.9 - 2023-06-15
* remove duplicated code

## gDRviz 1.6.8 - 2023-06-09
* simplify/standardize conversion of assays in gDRviz

## gDRviz 1.6.7 - 2023-05-23
* switch data.frame to data.table

## gDRviz 1.6.6 - 2023-05-13
* fix issue with login to DSDB 

## gDRviz 1.6.5 - 2023-05-06
* update dependencies (latest versions with data.table-related changes)

## gDRviz 1.6.4 - 2023-04-26
* reset identifiers when enter gDRin

## gDRviz 1.6.3 - 2023-04-25
* remove `useShinyAlerts` call

## gDRviz 1.6.2 - 2023-03-24
* update reqs for gDRstyle

## gDRviz 1.6.1 - 2023-03-10
* add support for local data storage
* add metaproject with all projects

## gDRviz 1.5.32 - 2023-03-03
* add support for wide structure in `convert_se_assay_to_dt`

## gDRviz 1.5.31 - 2023-02-22
* move functions from gDRcomponents

## gDRviz 1.5.30 - 2023-02-07
* update menu icons

## gDRviz 1.5.29 - 2023-01-31
* remove fetching drugs/cell_lines

## gDRviz 1.5.28 - 2023-01-23
* add show/hide manage tab

## gDRviz 1.5.27 - 2023-01-16
* add test for switching tab

## gDRviz 1.5.26 - 2022-01-10
* handle visualization not displaying when missing tissue metadata

## gDRviz 1.5.25 - 2023-01-10
* update package dependencies and unit tests
* switch from iheatmapr to heatmaply

## gDRviz 1.5.24 - 2022-12-21
* fix issue with misleading button color

## gDRviz 1.5.23 - 2022-12-21
* improve UI/UX for missing drugs/cell lines (submit tab; gDRin)

## gDRviz 1.5.22 - 2022-12-13
* add css style for go button

## gDRviz 1.5.21 - 2022-12-12
* add hiding gDRviz tab
* fix: hover over info naming is incorrect

## gDRviz 1.5.20 - 2022-12-09
* reset all deps to default branches

## gDRviz 1.5.19 - 2022-12-08
* hide experiment picker when only one experiment in MAE

## gDRviz 1.5.18 - 2022-12-08
* unify global variables between `gDRviz` and `gDRin`

## gDRviz 1.5.17 - 2022-12-07
* fix issue with heatmap not displaying when changing metrics

## gDRviz 1.5.16 - 2022-12-07
* plugins are available for everyone

## gDRviz 1.5.15 - 2022-12-05
* add support for instance-specific configs

## gDRviz 1.5.14 - 2022-12-01
* update dependencies

## gDRviz 1.5.13 - 2022-11-30
* add `manageData` plugin

## gDRviz 1.5.12 - 2022-11-30
* fix issue with namespacing

## gDRviz 1.5.11 - 2022-11-29
* improve navigation from front page

## gDRviz 1.5.10 - 2022-11-29
* extend notification system (popup or email)

## gDRviz 1.5.9 - 2022-11-25
* improve rendering `Visualizate data`

## gDRviz 1.5.8 - 2022-11-24
* rename front page

## gDRviz 1.5.7 - 2022-11-23
* add support for two types of output plugins

## gDRviz 1.5.6 - 2022-11-10
* all visualization modules are converted to plugins

## gDRviz 1.5.5 - 2022-10-28
* update solution for output plugins

## gDRviz 1.5.4 - 2022-10-25
* update gDR info links

## gDRviz 1.5.3 - 2022-10-24
* add initial version of front-page

## gDRviz 1.5.2 - 2022-10-24
* move some of gDRsearch staff to gDRsearch

## gDRviz 1.5.1 - 2022-10-21
* update version of `gDRcomponents`

## gDRviz 1.5.0 - 2022-10-20
* release bump

## gDRviz 1.4.12 - 2022-10-20
* fix: admin mode plugins toggle

## gDRviz 1.4.11 - 2022-10-18
* fix issue with multi-value identifiers

## gDRviz 1.4.10 - 2022-10-18
* feat: plugins available only in admin mode

## gDRviz 1.4.9 - 2022-10-18
* update version of `gDRcore`

## gDRviz 1.4.8 - 2022-10-17
* update version of `gDRimport`

## gDRviz 1.4.7 - 2022-10-11
* fix: downgrade to plotly (== 4.9.4.1)

## gDRviz 1.4.6 - 2022-10-10
* update version of `gDRcore`

## gDRviz 1.4.5 - 2022-10-03
* refactor logic for selected_metric across modules

## gDRviz 1.4.4 - 2022-10-03
* add helper js functions

## gDRviz 1.4.3 - 2022-09-26
* implement plugins

## gDRviz 1.4.2 - 2022-09-22
* data table download will server all pages

## gDRviz 1.4.1 - 2022-09-22
* use latest version of gDRtestData
* add gDRsearch to Docker image

## gDRviz 1.4.0 - 2022-09-12
* make new release

## gDRviz 1.3.11 - 2022-09-07
* fix issue with data source

## gDRviz 1.3.10 - 2022-08-29
* improve way of fetching project/multiproject data

## gDRviz 1.3.9 - 2022-08-23
* add 'restart' button and show 'logout' conditionally (shinyproxy only)

## gDRviz 1.3.8 - 2022-08-23
* move admin mode to gDRcomponents

## gDRviz 1.3.7 - 2022-08-18
* add Admin module

## gDRviz 1.3.6 - 2022-08-10
* add Admin module

## gDRviz 1.3.5 - 2022-08-09
* change the logic for displaying project name in the mini summary

## gDRviz 1.3.4 - 2022-08-01
* improve logic for setting cache dir with DSDB datasets

## gDRviz 1.3.3 - 2022-07-21
* add support for download of processed data in CSV and Excel format

## gDRviz 1.3.2 - 2022-06-15
* add support for missing data on aws

## gDRviz 1.3.1 - 2022-06-13
* allow larger number of drugs (3000) to be selected
* autosize clustering plot
* decrease width of dose-response curve in clustering module

## gDRviz 1.3.0 - 2022-06-02
* release

## gDRviz 1.2.45 - 2022-05-27
* change the way custom identifiers are handled

## gDRviz 1.2.44 - 2022-05-26
* use filtered data for mini summary

## gDRviz 1.2.43 - 2022-05-26
* use custom identifiers as env identifiers if necessary

## gDRviz 1.2.42 - 2022-05-25
* remove Drug2 from option as merging for combo data in manageData
* update the plots and menus to show the merged data options

## gDRviz 1.2.41 - 2022-05-16
* move filter_mae module to gDRcomponents

## gDRviz 1.2.40 - 2022-05-10
* enable cache on PVC
* DSDB as default data source

## gDRviz 1.2.39 - 2022-05-10
* update tests in helpers-drawMetricsDistribution.R 

## gDRviz 1.2.38 - 2022-05-05
* update unit tests and add support for polymapped identifiers

## gDRviz 1.2.37 - 2022-04-28
* add support for DSDB

## gDRviz 1.2.36 - 2022-04-27
* add PVC for deployment to `test`

## gDRviz 1.2.35 - 2022-04-27
* update logic for selecting per-group assay

## gDRviz 1.2.34 - 2022-04-22
* get rid of `gDRwrapper`

## gDRviz 1.2.33 - 2022-04-19
* add passing info about selected metric in grid

## gDRviz 1.2.32 - 2022-04-11
* add support for merging `Drugs` with additional metadata

## gDRviz 1.2.31 - 2022-04-04
* switch to a primary variable with more than two values
* add support for replicates/duration as additional variables

## gDRviz 1.2.30 - 2022-03-24
* update dependencies

## gDRviz 1.2.29 - 2022-03-21
* switch completely to MAE datasets only
    * temporarily support SE (but converted to MAE)
    * fix subsetting bug

## gDRviz 1.2.28 - 2022-03-14
* add support for additional perturbations

## gDRviz 1.2.27 - 2022-03-11
* add logout button

## gDRviz 1.2.26 - 2022-01-28
* enable proper handling of MAE datasets with different experiments

## gDRviz 1.2.25 - 2022-01-22
* refactor/clean-up server.R

## gDRviz 1.2.24 - 2022-01-21
* switch from SE to mae in server.R (I-step)

## gDRviz 1.2.23 - 2022-01-21
* fix logic for MAE/SE (combo module not properly displayed)

## gDRviz 1.2.22 - 2022-01-19
* change identifier `drugname` to `drug_name`

## gDRviz 1.2.21 - 2022-01-19
* move `helpers-manageData.R` to `gDRcomponents`

## gDRviz 1.2.20 - 2022-01-17
* new `gDRsearch` version `1.0.21`

## gDRviz 1.2.19 - 2022-01-13
* update unit tests to use MAE from `gDRtestData`

## gDRviz 1.2.18 - 2022-01-12
* temporarily turn off BiocFileCache for SE data

## gDRviz 1.2.17 - 2022-01-11
* move 'drawDrugCombo' from gDRviz to gDRcomponents

## gDRviz 1.2.16 - 2022-01-04
* update info button to include a drop down menu with multiple links

## gDRviz 1.2.15 - 2022-01-03
* 'combo_with_metrics' dataset works properly with 'combo module'

## gDRviz 1.2.14 - 2021-12-30
* fix linter issues in tests

## gDRviz 1.2.13 - 2021-12-27
* fix bug with creating nested list of experiments

## gDRviz 1.2.12 - 2021-12-22
* clean-up of the combo-related code
* combo-related helper functions move to gDRutils
* obsolete functions removed

## gDRviz 1.2.11 - 2021-12-22
* add support for MAE

## gDRviz 1.2.10 - 2021-12-12
* move selectMetric to gDRcomponents

## gDRviz 1.2.9 - 2021-12-06
* feat: click event on heatmap in `Metric Clustering` renders dose response plot

## gDRviz 1.2.8 - 2021-12-01
* fix bug in Metric Clustering, the notification "not enough data for clustering" still present after switching to good data

## gDRviz 1.2.7 - 2021-11-29
* fix bug with wrong RGB format

## gDRviz 1.2.6 - 2021-11-26
* add new tab - Manage data

## gDRviz 1.2.5 - 2021-11-24
* remove axis range edit for all plots

## gDRviz 1.2.4 - 2021-11-23
* in metrics contrast, second primary drop down, the default selection is the second option

## gDRviz 1.2.3 - 2021-11-23
* fix: daughter heatmaps don't render after switching Combo Metric and selecting same cellline and drugs

## gDRviz 1.2.2 - 2021-11-19
* feat: add isobologram ratio plot in combo module

## gDRviz 1.2.1 - 2021-11-19
* switch from browseURL to `htmltools::tagAppendAttributes`

## gDRviz 1.2.0 - 2021-11-19
* release 1.2.0

## gDRviz 1.1.45 - 2021-11-19
* use prettified assay names in (combo) mother heatmap

## gDRviz 1.1.44 - 2021-11-19
* improve support customized description in preDashboard

## gDRviz 1.1.43 - 2021-11-19
* feat: increase column width in View and Download Data table

## gDRviz 1.1.42 - 2021-11-18
* fix move button with the URL to combo presentation to the left

## gDRviz 1.1.41 - 2021-11-18
* feat add a button with the URL to combo presentation

## gDRviz 1.1.40 - 2021-11-17
* use plotly to generate pie charts

## gDRviz 1.1.39 - 2021-11-16
* fix issue when switching between combo and single agent data in `gDRsearch` hides irrelevant `gDRviz` menu item tabs

## gDRviz 1.1.38 - 2021-11-16
* trigger reset when data were selected in another module

## gDRviz 1.1.37 - 2021-11-16
* fix issue with missing pie charts

## gDRviz 1.1.36 - 2021-11-16
* fix missing sorting in Metrics Ranking

## gDRviz 1.1.35 - 2021-11-15
* fix issue when switching from drug to treatment in the summary tab

## gDRviz 1.1.34 - 2021-11-15
* set hover background in heatmap

## gDRviz 1.1.33 - 2021-11-15
* add hover over for combo metric descriptions

## gDRviz 1.1.32 - 2021-11-12
* fix improve on hover events in plotlyDC1

## gDRviz 1.1.31 - 2021-11-12
* fix error with variable mixup in `gDRcomponents`

## gDRviz 1.1.30 - 2021-11-12
* switch from `rId`, `cId` to `gDRcomponents::get_filtered_out_cols_combos()`

## gDRviz 1.1.29 - 2021-11-12
* add loading app spinner

## gDRviz 1.1.28 - 2021-11-12
* fix in combos: app crashes when no data for selected cellline and drugs

## gDRviz 1.1.27 - 2021-11-11
* create a new combo matrix tab

## gDRviz 1.1.26 - 2021-11-11
* fix edge-case when no project selected for the new metaproject
* fix: improve logic for replacing `untreated_tags` in co-drug data

## gDRviz 1.1.25 - 2021-11-11
* fix: selection update after exit dropdown in dashboard - gDRsearch

## gDRviz 1.1.24 - 2021-11-11
* filter rId and cId from combo assays

## gDRviz 1.1.23 - 2021-11-11
* remove unavailable modules completely from sidebar

## gDRviz 1.1.22 - 2021-11-10
* use the same metric through all modules

## gDRviz 1.1.21 - 2021-11-09
* fix bug GDR-1223: unchecking the last "iso levels" does not result in reactive change to daughter plots

## gDRviz 1.1.20 - 2021-11-08
* fix bug GDR-1210 where mother heatmap disappears after changing metrics

## gDRviz 1.1.19 - 2021-11-08
* make gDRviz compatible with shiny server on GRAN

## gDRviz 1.1.18 - 2021-11-05
* fix issues with matrix clustering when dimensions of data < 2L

## gDRviz 1.1.17 - 2021-11-05
* fix: observe state$cellline

## gDRviz 1.1.16 - 2021-11-04
* update color schemes for combo heatmaps
* simplify handling color schemes with 'get_combo_col_settings'

## gDRviz 1.1.15 - 2021-11-04
* all values from the heatmap available in hover over

## gDRviz 1.1.14 - 2021-11-02
* add collapsible to "assign box" in Drug Combo Treatments
* change to alphabetically sort in dropdown
* remove 'beta' flag

## gDRviz 1.1.13 - 2021-11-02
* re-render daughter heatmaps when ISO levels has change

## gDRviz 1.1.12 - 2021-11-02
* fix issue with lower limit being NA

## gDRviz 1.1.11 - 2021-10-29
* fix issue with missing concentration2 column

## gDRviz 1.1.10 - 2021-10-29
* use full space for options

## gDRviz 1.1.9 - 2021-10-29
* update admins list 

## gDRviz 1.1.8 - 2021-10-29
* improve colorscaling in plotlyDC1 and plotlyDCMother
* diagonal black line removed (plotlyDC1)
* improve labels for isobolograms (plotlyDC1)
* fix the bug in text labels (value tag in plotlyDC1)
* improve logic for drawing isobologram points (disabled by default)
* update color sets for isobolograms

## gDRviz 1.1.7 - 2021-10-28
* [Drug Combo Treatments] move selector to the top and made collapsible
* [Drug Combo Treatments] replaced modal with box below general heatmap
* [Drug Combo Treatments] app scrolls down to mini heatmaps box when user clicks on general heatmap
* [Drug Combo Treatments] increased size of mini heatmaps
* add combo data to View and Download Data module
* all modules disabled except [View and Download data] and [Drug Combo Treatments] for combo data

## gDRviz 1.1.6 - 2021-10-27
* make miniSESummary more user-friendly 

## gDRviz 1.1.5 - 2021-10-26
* change sidebar width

## gDRviz 1.1.4 - 2021-10-25
* add support for fetching project id from environmental variable

## gDRviz 1.1.3 - 2021-10-21
* add shiny modules: `drawDrugCombo_UI`, `drawDrugCombo_SERVER`

## gDRviz 1.1.2 - 2021-10-21
* add Apply button to update heatmap colors

## gDRviz 1.1.1 - 2021-10-21
* add helper functions for drawDrugCombo modules

## gDRviz 1.0.36 - 2021-10-20
* fix issue with loading examples list #251

## gDRviz 1.0.35 - 2021-10-19
* the plot was transformed into a square one in Metric Contrast

## gDRviz 1.0.34 - 2021-10-07
* fix bug with `res_cor` variable

## gDRviz 1.0.33 - 2021-10-04
* issue with non-clickable arrows in assign box tittle

## gDRviz 1.0.32 - 2021-10-01
* fix: overlap of xaxis title and plot title when unchecking Show X axis labels

## gDRviz 1.0.31 - 2021-09-24
* add namespace imports for shinyFeedback causing app crashing when trying colors in the Dose Response Overview

## gDRviz 1.0.30 - 2021-09-22
* modal for users to select colors for the heatmap in metric clustering
* the row and column annotations use fix color palettes for consistency when switching metrics

## gDRviz 1.0.29 - 2021-09-15
* switch from `get_identifier` to `get_SE_identifiers`

## gDRviz 1.0.28 - 2021-09-13
* refactor: include the `simplify` argument in get_(prettified|env)_identifiers functions

## gDRviz 1.0.27 - 2021-09-09
* fix: handle properly edge case with codrugs data (drawMetricsRanking and drawMetricsContrast modules)

## gDRviz 1.0.26 - 2021-09-07
* fix issue with dropdowns in Metric Contrast

## gDRviz 1.0.25 - 2021-09-06
* fix: make 'draw response curves' pop-up working again

## gDRviz 1.0.24 - 2021-09-02
* add "Color by" option and removed coloring from "Group by" in Metrics Distribution

## gDRviz 1.0.23 - 2021-09-01
* gDRcomponents upgraded to 1.0.10 ('Metric Clustering' improvements for data with NAs)

## gDRviz 1.0.23 - 2021-09-01
* temporary fix for 'c50/ec50' issue

## gDRviz 1.0.21 - 2021-08-27
* block gDRviz for too big SE

## gDRviz 1.0.20 - 2021-08-24
* get rid of hardcoded identifiers
* move some functions from gDRviz into gDRcomponents

## gDRviz 1.0.19 - 2021-08-20
* add collapsible to "assign box" in 3 modules

## gDRviz 1.0.18 - 2021-08-16
* add ability to set range of heatmaps

## gDRviz 1.0.18 - 2021-08-17
* limit plotly editable components to only axis and plot titles

## gDRviz 1.0.17 - 2021-08-15
* add multiple independent dropdowns in Metric Contrast

## gDRviz 1.0.16 - 2021-07-06
* switch to refactored version of `gDRcomponents`

## gDRviz 1.0.15 - 2021-06-29
* switch to `gDRtestData` version without sensitive data

## gDRviz 1.0.14 - 2021-06-23
* refactor unit tests data

## gDRviz 1.0.13 - 2021-06-16
* move helper functions and add unit tests for drawResponseCurves
* rename `module-drawResponseCurves` to `module-drawResponseOverview`
* rename `module-drawResponseCurves1` to `module-drawResponseCurves`

## gDRviz 1.0.12 - 2021-06-16
* move helper functions and add unit tests for drawResponseGrid

## gDRviz 1.0.11 - 2021-06-15
* fix error in 'Modify and View data' module

## gDRviz 1.0.10 - 2021-06-15
* fix error in clustering view

## gDRviz 1.0.9 - 2021-06-14
* move helpers and add unit test for selectGroups

## gDRviz 1.0.8 - 2021-06-11
* drop gDR dependency

## gDRviz 1.0.7 - 2021-06-10
* move helper functions and add unit tests for drawMetricsDistribution

## gDRviz 1.0.6 - 2021-06-09
* allow for switching between different metrics types
* refactor module manageData

## gDRviz 1.0.5 - 2021-06-07
* move helper functions and add unit tests for drawMetricsRanking

## gDRviz 1.0.4 - 2021-06-07
* update gDRsearch dependencies with alphabetical order of drugs/cell lines in drop-down menu

## gDRviz 1.0.3 - 2021-06-04
* remove disambiguation of rowData for merged single-agent data with duplications caused by multiple data sources

## gDRviz 1.0.2 - 2021-06-02
* improve the logic for the clustering in the MetricClustering.

## gDRviz 1.0.1 - 2021-06-01
* add support for vectors in sprintf for manageData module
* add missing namespace in logs.R

## gDRviz 1.0.0 - 2021-06-01
* release 1.0.0

## gDRviz 0.1.38 - 2021-06-01
* improve logic for namespacing

## gDRviz 0.1.37 - 2021-06-01
* fix issue in View and Download Data

## gDRviz 0.1.36 - 2021-05-31
* use reactive values for rendering warning message

## gDRviz 0.1.35 - 2021-05-31
* keep the notification in clustering module on place

## gDRviz 0.1.34 - 2021-05-27
* prepare more specific solution for combining data

## gDRviz 0.1.33 - 2021-05-27
* set secondary drug as factor
* add text in parentheses to the right of 'Choose Metric'

## gDRviz 0.1.32 - 2021-05-26
* improve error message during the clustering

## gDRviz 0.1.31 - 2021-05-26
* restore changes removed in GDR-880 [Metric Contrast]

## gDRviz 0.1.29 - 2021-05-24
* do not expect all SEs to be single-agent

## gDRviz 0.1.28 - 2021-05-24
* allow nested "data_source" column to be combined with drug or cell line

## gDRviz 0.1.27 - 2021-05-20
* change projects order

## gDRviz 0.1.26 - 2021-05-20
* add tooltips to 'Select Metric' radio buttons

## gDRviz 0.1.25 - 2021-05-20
* move helper functions and add unit tests for module-drawMetricsContrast.R

## gDRviz 0.1.24 - 2021-05-20
* logs are hidden to user

## gDRviz 0.1.23 - 2021-05-19
* add logic for validating SE objects in gDRviz

## gDRviz 0.1.22 - 2021-05-17
* remove substring on table columns
* enable filtering on table columns

## gDRviz 0.1.21 - 2021-05-17
* update logic for drawing 'correlation line' for log-transformed data (GDR-927) 

## gDRviz 0.1.20 - 2021-05-17
* fix bug with "Metrics_rownames"

## gDRviz 0.1.19 - 2021-05-17
* add custom limit of drugXcell_lines combinations

## gDRviz 0.1.18 - 2021-05-14
* update logic for legend
* add legend titles

## gDRviz 0.1.17 - 2021-05-14
* add zero line for vertical and horizontal view for "Metric Contrast"

## gDRviz 0.1.16 - 2021-05-13
* minor wording improvements

## gDRviz 0.1.15 - 2021-05-13
* round numeric data in 'View and Modify' Tab

## gDRviz 0.1.14 - 2021-05-11
* metric Ranking: remove 'rank by'

## gDRviz 0.1.13 - 2021-05-10
* add coloring option for "Metric Contrast" module

## gDRviz 0.1.12 - 2021-05-10
* do not show tabs until data is ready

## gDRviz 0.1.11 - 2021-05-10
* hide example tab

## gDRviz 0.1.10 - 2021-05-10
* show axes for log-transformed data

## gDRviz 0.1.9 - 2021-05-06
* fix bug with heatmap labels

## gDRviz 0.1.8 - 2021-05-06
* make identity line always crossing (0,0)

## gDRviz 0.1.7 - 2021-05-06
* improve logic for legend in 'Metrics Ranking'
    * do not show 'codrug concentration' in the legend in the case of single-agent data

## gDRviz 0.1.6 - 2021-05-04
* add unit tests for helpers-drawMetricsClustering.R

## gDRviz 0.1.5 - 2021-05-04
* add info that 'View all projects' tab is in beta phase
* refactor renaming of columns in assays
* fix bug with subplot in plotly

## gDRviz 0.1.4 - 2021-05-03
* dose Response Overview: major rewrite of internals, UI changes, bug fixes

## gDRviz 0.1.3 - 2021-04-30
* switch from `assay_to_dt` to `convert_se_assay_to_dt`
* add Concentration + Concentration_2 in Normalized and Averaged views
* change data in simple search

## gDRviz 0.1.2 - 2021-04-30
* make  identity line  being drawn properly for datasets with combo data 

## gDRviz 0.0.113 - 2021-04-29
* use RDS files from gDRtestData package

## gDRviz 0.0.112 - 2021-04-28
* update project order

## gDRviz 0.0.111 - 2021-04-23
* improve detection of min/max values of the y-axis

## gDRviz 0.0.110 - 2021-04-23
* 's/Generating SE/Loading the project/'

## gDRviz 0.0.109 - 2021-04-22
* move non-Shiny functions from module-drawMetricsClustering.R to helpers-drawMetricsClustering.R

## gDRviz 0.0.108 - 2021-04-21
* improve logic for legend in 'Metrics Ranking'

## gDRviz 0.0.107 - 2021-04-19
* show legend in 'Metrics Ranking' if user colors by 'Drug MOA'

## gDRviz 0.0.106 - 2021-04-16
* make development quicker by allowing a `GDRVIZ_TEST_DATA` envvar to automatically load a dataset

## gDRviz 0.0.105 - 2021-04-16
* dose response overview: reset colors after switching metric

## gDRviz 0.0.104 - 2021-04-16
* metric distribution: allow showing/hiding x axis labels

## gDRviz 0.0.103 - 2021-04-16
* add minor changes in gDRsearch

## gDRviz 0.0.102 - 2021-04-15
* metric contrast: Improved logic for drawing identity line

## gDRviz 0.0.101 - 2021-04-15
* dose response overview: fix whitespace on left column

## gDRviz 0.0.100 - 2021-04-15
* modify and view data: allow table to take up full width

## gDRviz 0.0.99 - 2021-04-14
* fix bug in 'Metric Clustering' (non-numeric variable(s) in data frame: Drug)

## gDRviz 0.0.98 - 2021-04-13
* metric clustering: fix hover when transpose

## gDRviz 0.0.97 - 2021-04-12
* refactor disambiguate functions

## gDRviz 0.0.96 - 2021-04-12
* metric contrast: make sure plot has same scale on x and y axes

## gDRviz 0.0.95 - 2021-04-12
* metric clustering: fix hover when transpose

## gDRviz 0.0.94 - 2021-04-12
* add threshold for maximum number of drugsXcell_lines combinations in single search

## gDRviz 0.0.93 - 2021-04-09
* switch from `acast` to `dcast` in module drawMetricsClustering

## gDRviz 0.0.92 - 2021-04-07
* metric ranking: allow plot to be fluid width, add option to show/hide x labels

## gDRviz 0.0.91 - 2021-04-07
* standardize use of header names as per gDRutils::get_headers()

## gDRviz 0.0.90 - 2021-04-06
* fix a bug with metric ranking

## gDRviz 0.0.89 - 2021-04-02
* metric clustering: Allow user to transpose heatmap

## gDRviz 0.0.88 - 2021-04-02
* change tabs order ('dose response overview' preceding 'dose response curves')

## gDRviz 0.0.87 - 2021-04-02
* move options to above plots to allow more room for plots

## gDRviz 0.0.86 - 2021-04-02
* change order of columns

## gDRviz 0.0.85 - 2021-03-31
* bugfix jumbled drug names in metrics clustering

## gDRviz 0.0.84 - 2021-03-30
* view Data shows unprocessed data

## gDRviz 0.0.83 - 2021-03-30
* add comments explaining NA values for EC50 metrics

## gDRviz 0.0.82 - 2021-03-30
* fix problem with global envs

## gDRviz 0.0.81 - 2021-03-29
* add fix for character data in metrics assay
* update data in simple search

## gDRviz 0.0.80 - 2021-03-26
* add conditional loading of restApiSession module (when REST API authorization enabled)

## gDRviz 0.0.78 - 2021-03-23
* bugfix assay_metrics column names by substituting data.table::setcolorder for base R equivalent 

## gDRviz 0.0.77 - 2021-03-20
* remove old dashboard
* add new version of Simple Search

## gDRviz 0.0.76 - 2021-03-18
* metric distribution: 'Group by' options did not reset when a new dataset (experiment) was chosen

## gDRviz 0.0.75 - 2021-03-18
* metric clustering: implement a workaround for bug that prevented heatmap from rendering after a plotly plot

## gDRviz 0.0.74 - 2021-03-17
* metric distribution: give more space to the plot by moving all inputs to the same row
* metric distribution: clarify that 'color by' is really 'group by' and what variable is used for the x axis
* metric distribution: disallow some options from being used as 'group by' if the resulting plot is too confusing

## gDRviz 0.0.73 - 2021-03-16
* change col order
* update BumpyMatrix and gDRutils packages
* update synthetic data

## gDRviz 0.0.72 - 2021-03-15
* fix bug where the mini se summary row would not always update at the correct times

## gDRviz 0.0.71 - 2021-03-15
* drug response curves: ability to deselect by clicking in same spot (GDR-718)

## gDRviz 0.0.70 - 2021-03-15
* remove explicit loading of {iheatmapr}

## gDRviz 0.0.69 - 2021-03-15
* add new synthetic data

## gDRviz 0.0.68 - 2021-03-09
* internal changes.

## gDRviz 0.0.67 - 2021-03-08
* remove the codilution and matrix data

## gDRviz 0.0.66 - 2021-03-03
* update README

## gDRviz 0.0.65 - 2021-03-02
* add new gDRsearch dashboard with AWS support

## gDRviz 0.0.64 - 2021-02-17
* add synthetic data

## gDRviz 0.0.63 - 2021-02-04
* add Drug Combo Treatment

## gDRviz 0.0.62 - 2021-01-27
* bug fixes and internal changes.

## gDRviz 0.0.61 - 2021-01-20
* add simpleQuerySearch from gDRsearch

## gDRviz 0.0.60 - 2021-01-18
* bug fixes.

## gDRviz 0.0.59 - 2021-01-11
* minor tweaks to UI.
* internal changes.

## gDRviz 0.0.58 - 2021-01-07
* more metric choices. Some metrics have been renamed (capitalized).
* ranking View: labels on X axis will disappear if they cannot all fit in the plot.
* ranking View: data can be grouped and colored on the same variable.
* metric axes in Distribution, Ranking and Contrast Views can be fixed or bottom-released.
* enhance annotations in Clustering View.
* some UI improvements.
* some minor updates to plots.
* internal changes.

## gDRviz 0.0.57 - 2020-12-24
* bug fix.

## gDRviz 0.0.56 - 2020-12-22
* add data management tab. Data can be modified, viewed, and saved.

## gDRviz 0.0.55 - 2020-12-18
* internal changes.

## gDRviz 0.0.54 - 2020-12-17
* Ace Editor has restricted access.
* bug fixes and minor changes to interface.

## gDRviz 0.0.53 - 2020-12-11
* new button for building preparing data.
* data is also created automatically when moving from gDRsearch to gDRviz.

## gDRviz 0.0.52 - 2020-12-08
* clustering View now uses Spearman correlation coefficients for distance calculation. Missing values are replaced with 0, meaning no correlation.
* fix issue with missing colors in Disatribution View.
* minor improvements.

## gDRviz 0.0.51 - 2020-12-07
* update gDRviz version.

## gDRviz 0.0.50 - 2020-12-04
* add gDR logo.

## gDRviz 0.0.49 - 2020-12-04
* add aceEditor to allow launching R/bash console from the app.

## gDRviz 0.0.48 - 2020-12-04
* fix issues with adding colors and button behavior in Grid View.
* curve simulation in curve view has more definition now, but plotting will be slightly slower.
* minor bug fixes.

## gDRviz 0.0.47 - 2020-12-03
* add fifth layer to Curve View that highlights max effectiveness and half-effective concentrations.
* modifie dashboard slightly.
* remove some buttons from plot mode bars.
* internal changes and bug fixes.

## gDRviz 0.0.46 - 2020-11-24
* buttons for adding color to curves in Grid View now have dynamic background colors.
* improve warnings on dropping data in Contrast, Clustering and Ranking views.
* clustering View now only drops individual NA values rather than whole rows that contain any. Rows or columns containing nothing but NAs are dropped.
* plots are editable now. Click on title or axis label to change it at will. HTML tags are accepted, e.g. `<b>`.
* infinite values of GR50, ic50 and ec50 are now high-capped at 30 uM, regardless of other experiment parameters.
* internal changes and improvements, and some bug fixes.

## gDRviz 0.0.45 - 2020-11-23
* switch to R 4.0.

## gDRviz 0.0.44 - 2020-11-18
* update `gDRsearch` version.

## gDRviz 0.0.43 - 2020-11-17
* bug fix.

## gDRviz 0.0.42 - 2020-11-17
* update `gDRsearch`.

## gDRviz 0.0.41 - 2020-10-12
* add color control in Grid View. Users can now select/build color palettes or pick colors for curves one by one.
* critical bug fix.

## gDRviz 0.0.40 - 2020-10-11
* add data set number 10.
* slight changes to UI in selection of example data sets.

## gDRviz 0.0.39 - 2020-11-06
* add highlighting curves to grid view.

## gDRviz 0.0.38 - 2020-10-30
* add data set number 9.

## gDRviz 0.0.37 - 2020-10-27
* add example data set number 8.
* extend access restriction to `prd` environment.

## gDRviz 0.0.36 - 2020-10-27
* update `gDRsearch`. 

## gDRviz 0.0.35 - 2020-10-26
* restore standalone Response Curve view.
* fix spacing of boxes/violins in Metric Distribution view.
* esthetic improvements:
    * human friendly variable names
    * correlation statistics in Metric Contrast view
    * axis and plot titles in Response Curve and Response Grid views
    * more grouping and coloring options in Metric Distribution view
    * GR50, ic50 and ec50 are plotted raw (not log transformed) but on log scale
    * linear fit in Metric Contrast view disabled for metrics plotted in log scale
* bug fixes.

## gDRviz 0.0.34 - 2020-10-14
* add DATA_SOURCE for loading data.

## gDRviz 0.0.33 - 2020-10-13
* restrict access to `dev` environment to Genentech users and dev team.

## gDRviz 0.0.32 - 2020-10-08
* enter beta mode. App deployed on `dev` and `prd` environments thanks to `RPlatform`.

## gDRviz 0.0.31 - 2020-10-08
* add new data. 

## gDRviz 0.0.30 - 2020-10-08
* add panel titles in Response Grid view.
* add titles to Metric Clustering and Metric Ranking views.
* expand metric selection in Distribution, Contrast, Clustering and Ranking views.
* change format of exported images to `.svg`.
* bug fixes.

## gDRviz 0.0.29 - 2020-10-02
* add Metric Ranking view.
* bug fixes.

## gDRviz 0.0.28 - 2020-10-01
* disambiguation extended to drug names.

## gDRviz 0.0.27 - 2020-09-30
* minor fixes.

## gDRviz 0.0.26 - 2020-09-28
* add drug MOA information to Metric Clustering view.

## gDRviz 0.0.25 - 2020-09-24
* in Response Grid a double click on a panel opens a pop-up with the Response Curve view.
* Response Curve view splits second plot into panels, one for each co-drug.
* Response Curve view remove as a separate tab, now only appears within the Response Grid view.

## gDRviz 0.0.24 - 2020-09-18
* add basic tooltips to Response Grid view.

## gDRviz 0.0.23 - 2020-09-17
* loaders add to all visualizations.
* change appearance of data summary.
* minor fixes.

## gDRviz 0.0.22 - 2020-09-16
* major update to accommodate co-treatments in the existing visualizations.
* in Response Curve view a separate curve is plotted for every co-drug at every concentration.
* in Metric Contrast view points get different colors and shapes depending on co-treatments.
* in Metric Clustering view every co-treatment is treated as a separate drug and is add to the heatmap as a separate row.
* co-drug name and concentrations are add to data summary.
* co-drug names and concentrations add to tooltips.
* plots in Response Curve view now also shows GRvalue or RelativeViability.
* duplicate names of cell lines in data sets are disambiguated by appending clids.
* minor changes.

## gDRviz 0.0.21 - 2020-09-15
* recalculate example data and add three more data sets.
* upon selection of a data set, a summary is displayed.

## gDRviz 0.0.20 - 2020-09-09
* add loaders to Response Curve and Response Grid views

## gDRviz 0.0.19 - 2020-09-07
* add Response Grid view. Data is split by drugs or cell lines and a panel is drawn for every group. Mouse-over highlights corresponding curves in all panels.

## gDRviz 0.0.18 - 2020-09-07
* enhance Response Curve view. Two plots are drawn, the first only shows curves and selects data for the second one, which can show other data types.

## gDRviz 0.0.17 - 2020-08-31
* internal changes.

## gDRviz 0.0.16 - 2020-08-28
* internal changes
* bug fixes

## gDRviz 0.0.15 - 2020-08-28
* log transformation restricted to metrics "ending with '50'".

## gDRviz 0.0.14 - 2020-08-28
* capping GR50 and ic50 to get rid of infinite values.

## gDRviz 0.0.13 - 2020-08-27
* changes to plotting response curves in Response Curve view.
* mask data points are dropped automatically. Button removed.

## gDRviz 0.0.12 - 2020-08-24
* add loader in the Response Curve view.

## gDRviz 0.0.11 - 2020-08-20
* change plotting tool.
* add highlights to Response Curve view.

## gDRviz 0.0.10 - 2020.08.12
* internal changes.

## gDRviz 0.0.9 - 2020.08.07
* add Metric Clustering view.

## gDRviz 0.0.8 - 2020.08.06
* tweaks to interface.

## gDRviz 0.0.7 - 2020.08.05
* add option to remove masked data in Response Curve view.
* adjust limits of Y axis in Response Curve view.

## gDRviz 0.0.6 - 2020.08.03
* add buttons in Response Curve view.

## gDRviz 0.0.5 - 2020.08.03
* new version of Response Curves view, shows only one panel and switches between 1 cell line and any drugs or vice versa.

## gDRviz 0.0.4 - 2020.08.03
* internal changes.

## gDRviz 0.0.3 - 2020.07.30
* new visualization: Response Curves. Plot can be split into panels.

## gDRviz 0.0.2 - 2020.07.30
* new dashboard for visualizations.
* pick one of 4 example experiments.
* plot Metric Distribution and Metric Contrast.
