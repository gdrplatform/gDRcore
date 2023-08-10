# Change to v.0.99.33 - 2023-08-10
- update annotation column names for cell line annotation as per changes in the gDRutils

# Change to v.0.99.32 - 2023-07-25
- extended logic for supporting cols with dash, e.g. additional perturbations with "-"

# Change to v.0.99.31 - 2023-07-19
- update the logic for handling warnings in the pipeline

# Change to v.0.99.30 - 2023-07-13
- fix issue with wrong merging of data.tables without nested confounders

# Change to v.0.99.29 - 2023-07-07
- add information about source type for cases without metric data
- refactor the logic for splitting raw data from metadata (get rid of iterative approach)

# Change to v.0.99.28 - 2023-07-05
- update logic for parallel computing

# Change to v.0.99.27 - 2023-06-29
- optimize unit tests

# Change to v.0.99.26 - 2023-06-29
- remove backward compatibility for old data model

# Change to v.0.99.25 - 2023-06-27
- fix bug with missing rownames in normalized assay

# Change to v.0.99.24 - 2023-06-19
- update logic for merging data.table objects

# Change to v.0.99.23 - 2023-06-13
- replace `order` with `data.table::setorder`

# Change to v.0.99.22 - 2023-06-09
- switch from `merge` to `[[` for `data.table` objects

# Change to v.0.99.21 - 2023-06-06
- switch from `aggregate` to `data.table`

# Change to v.0.99.20 - 2023-06-07
- switch from `zoo::rollmean` to `data.table::frollmean`

# Change to v.0.99.19 - 2023-06-06
- replaced reshape2 functions by functions from data.table

# Change to v.0.99.18 - 2023-05-31
- fix managing of mixed types of raw data

# Change to v.0.99.17 - 2023-05-29
- fix bug with subsetting data for calculating isobologram

# Change to v.0.99.16 - 2023-05-22
- format the vignette with BiocStyle

# Change to v.0.99.15 - 2023-05-16
- fix related with data.table

# Change to v.0.99.14 - 2023-05-15
- rename `excess` to `x` to unify colnames in assay data

# Change to v.0.99.13 - 2023-05-10
- refactor normalization_types in combo-specific assays

# Change to v.0.99.12 - 2023-05-09
- Utilize `gDRutils::apply_bumpy_function` in fit_SE

# Change to v.0.99.11 - 2023-05-05
- fix bug with swapping untreated/vehicle values

# Change to v.0.99.10 - 2023-05-04
- Fix bug with data.table

# Change to v.0.99.9 - 2023-04-21
- Utilize `gDRutils::apply_bumpy_function` in average_SE

# Change to v.0.99.8 - 2023-04-20
- switch to OSI license

# Change to v.0.99.7 - 2023-04-19
- Fix bug with replacing vehicle to untreated values

# Change to v.0.99.6 - 2023-04-19
- Moved wrapper fuctions from gDRtestData

# Change to v.0.99.5 - 2023-04-18
- Update dependencies
- Add fix for bioc-devel (correct sorting in merge test)

# Change to v.0.99.4 - 2023-04-17
- Fix namespacing issue in examples
- Add R 4.2 as dependency
- Fix examples for normalize_SE

# Change to v.0.99.3 - 2023-04-12
- add logic for retrieving raw data from assay data

# Change to v.0.99.2 - 2023-04-07
- update maintainer

# Change to v.0.99.1 - 2023-04-04
- bugfix for the logic in 'cleanup_metadata'

# Change to v.0.99.0 - 2023-03-24
- Make the package Bioc-compatible

# Change to v.1.3.22 - 2023-03-21
- Improve performance of 'map_df' with refactored logic for exact matches

# Change to v.1.3.21 - 2023-03-15
- Refactor pipeline

# Change to v.1.3.20 - 2023-03-09
- Address co-treatment fit by using the matrix data type instead

# Change to v.1.3.19 - 2023-03-07
- Add support for splitting normalization types

# Change to v.1.3.18 - 2023-03-06
- improve logic in functions used to generate isobolograms' data

# Change to v.1.3.17 - 2023-03-06
- Remove obsolete code

# Change to v.1.3.16 - 2023-02-21
- Add support for partial pipeline runs

# Change to v.1.3.15 - 2023-02-10
- Update path to annotation data

# Change to v.1.3.14 - 2023-02-10
- Clean-up code

# Change to v.1.3.13 - 2023-01-13
- Clean-up code

# Change to v.1.3.12 - 2023-01-12
- Refactor mapping function to properly handle drug3

# Change to v.1.3.11 - 2022-12-16
- Replace parallelize function with gDRutils::loop

# Change to v.1.3.10 - 2022-12-15
- Add assert for vehicle values in input data in `runDrugResponseProcessingPipeline`

# Change to v.1.3.9 - 2022-12-14
- Fix error-handling if conditions in average_SE

# Change to v.1.3.8 - 2022-11-07
- Fix invalid encapsulation in tests

# Change to v.1.3.7 - 2022-10-18
- Add missing namespacing

# Change to v.1.3.6 - 2022-10-10
- Add support for setting number of cores for `BiocParallel` based on the env variable

# Change to v.1.3.5 - 2022-10-07
- Remove global parameters for `BiocParallel`

# Change to v.1.3.4 - 2022-09-29
- Change the logic for using cores in `BiocParallel`

# Change to v.1.3.3 - 2022-09-27
- Update the logic for parallel computing

# Change to v.1.3.2 - 2022-08-17
- Update the logic for default nested_confounders in `create_SE` function

# Change to v.1.3.1 - 2022-07-08
- Refactor create_SE to support reverse single-agent data

# Change to v.1.3.0 - 2022-06-02
- Release

# Change to v.1.1.39 - 2022-05-30
- add missing namespace for get_env_identifiers

# Change to v.1.1.38 - 2022-05-09
- replace `NA` by 0 in Concentration loaded in manifest file

# Change to v.1.1.37 - 2022-05-06
- switch from `data.table` to `data.frame` in add_annotation* functions

# Change to v.1.1.36 - 2022-04-12
- update function for adding unknown cell line annotations
- update logic for using nested confounders
- remove `grr` from dependencies

# Change to v.1.1.35 - 2022-04-07
- get rid of gDRwrapper

# Change to v.1.1.34 - 2022-03-30
- add support for additional barcode identifiers

# Change to v.1.1.33 - 2022-03-18
- fix documentation in calculate_matrix_metric

# Change to v.1.1.32 - 2022-03-17
- switch from `catchr` to `purrr`

# Change to v.1.1.31 - 2022-03-02
- fix mapping reference values of inverted treatments
- remove duplication of single-agent data

# Change to v.1.1.30 - 2022-02-17
- wrap SE into MAE on the level of runDrugResponseProcessingPipeline

# Change to v.1.1.29 - 2022-02-14
- change identifier `drugname` to `drug_name`

# Change to v.1.1.28 - 2022-02-09
- switch from `gDRinternal` to `gDRinternalData` for internal annotations

# Change to v.1.1.27 - 2022-02-08
- issue with subsetting by list in R 4.2.0

# Change to v.1.1.26 - 2022-02-02
- unlist output of intersect as per R 4.2.0

# Change to v.1.1.25 - 2022-01-31
- align version criteria between `dependencies.yaml` and `DESCRIPTION` package versions

# Change to v.1.1.24 - 2022-01-25
- standardize/improve CI

# Change to v.1.1.23 - 2022-01-24
- switch from `cores` variable to `detect_cores()` function

# Change to v.1.1.22 - 2022-01-24
- fix wrong type of `NUM_CORES` env variable

# Change to v.1.1.21 - 2022-01-03
- speed-up functions for mapping treated and untreated cases

# Change to v.1.1.20 - 2021-12-30
- fix linter issues

# Change to v.1.1.19 - 2021-12-14
- use parallel computing as an alternative for `for` loops

# Change to v.1.1.18 - 2021-12-07
- update annotation script as per new csv annotation files

# Change to v.1.1.17 - 2021-12-07
- detect co-trt data and treat them as single-agent

# Change to v.1.1.16 - 2021-11-08
- set excess = NA for single-agent

# Change to v.1.1.15 - 2021-10-29
- solve rounding issues
- add new bliss metric

# Change to v.1.1.14 - 2021-10-25
- refactor isobolograms

# Change to v.1.1.13 - 2021-10-25
- move p_trt_keys to the proper place

# Change to v.1.1.12 - 2021-10-20
- add support for masked data in fit_SE.combinations.R

# Change to v.1.1.11 - 2021-10-14
- address issues in creating the SE for combo matrix experiments

# Change to v.1.1.10 - 2021-10-13
- refactor the logic for calculation Loewe when there is no `Concentration == 0`

# Change to v.1.1.9 - 2021-10-12
- refactor the logic for combo data

# Change to v.1.1.8 - 2021-09-27
- updated normalization_types in 'calculate_combo_matrix' and 'fit_SE'

# Change to v.1.1.7 - 2021-09-21
- calculate_GR_value by removing cl_name param

# Change to v.1.1.6 - 2021-08-25
- fix but with missing `nested_identifiers` variables for creating DataFrame for masked values

# Change to v.1.1.5 - 2021-08-04
- refactor reading annotations and add default parameters

# Change to v.1.1.4 - 2021-07-30
- fix bug for nested_confounders not present in ref_df

# Change to v.1.1.3 - 2021-07-27
- add support for nested_identifiers and nested_confounders 

# Change to v.1.1.2 - 2021-07-23
- remove obsolete dependencies from DESCRIPTION:Imports

# Change to v.1.0.9 - 2021-07-01
- function for testing synthetic data

# Change to v.1.0.8 - 2021-06-25
- move functions for importing template files to gDRimport package 

# Change to v.1.0.7 - 2021-06-23
- add linter
- remove obsolete functions from test.utils.R

# Change to v.1.0.6 - 2021-06-22
- move importing functions to gDRimport package

# Change to v.1.0.5 - 2021-06-14
- add unit tests for synthetic data

# Change to v.1.0.4 - 2021-06-14
- remove deprecated functions and unit tests

# Change to v.1.0.3 - 2021-06-14
- update namespace for 'metadata()'

# Change to v.1.0.2 - 2021-06-10
- change package name (gDR => gDRcore)

# Change to v.1.0.1 - 2021-06-04
- export/update docs for 'standardize_record_values'

# Change to v.1.0.0 - 2021-06-02
- release 1.0.0

# Change to v.0.1.52 - 2021-04-29
- remove all the gDRinternal-related files

# Change to v.0.1.51 - 2021-04-27
- fix wrong argument name in `fit_curves`

# Change to v.0.1.50 - 2021-04-23
- Ensure that dts in `assay(se,"Averaged")` are NULL when there are not treatments

# Change to v.0.1.49 - 2021-04-20
- add processing info metadata to SE

# Change to v.0.1.48 - 2021-04-14
- improve handling of nested_keys and override_controls

# Change to v.0.1.47 - 2021-04-09
- get rid of sorting index columns in getMetaData

# Change to v.0.1.46 - 2021-03-31
- move SE metadata getters and setters to gDRutils

# Change to v.0.1.45 - 2021-03-31
- fix a bug in map_df function

# Change to v.0.1.44 - 2021-03-26
- refactor processing functions in create_SE2, normalize_SE2, average_SE2, and fit_SE2 and friends

# Change to v.0.1.43 - 2021-02-22
- add drug_moa_2 to merged df

# Change to v.0.1.42 - 2021-02-17
- fix a bug with adding annotation to codrugs

# Change to v.0.1.41 - 2021-02-04
- Add script for benchmarking [normalize/average/metrics]-SE functions

# Change to v.0.1.40 - 2021-02-02
- fix unit tests
- include discard keys info in getMetaData for rowData

# Change to v.0.1.39 - 2021-01-29

- utilize refactored and renamed fitting function RVGRfit to fit_curves

# Change to v.0.1.38 - 2021-01-27

- export some functions for use in gDRutils

# Change to v.0.1.37 - 2021-01-20

- move df_to_assay and df_to_bm_assay to gDRutils

# Change to v.0.1.36 - 2021-01-18

- updated rowData and colData in SE with additional columns

# Change to v.0.1.35 - 2021-01-14

- update positional header calls to explicit identifier calls
- minor changes to comply with gDRstyle
- minor changes to use identifiers over hard coded DrugName requirements

# Change to v.0.1.34 - 2021-01-12

- update scripts related to annotation of drugs and cell lines

# Change to v.0.1.33 - 2021-01-05
- anonymize test data

# Change to v.0.1.32 - 2020-12-24
- switch from GeneDataScreenR to gDRinternal from QCS to gDR functions

# Change to v.0.1.31 - 2020-12-22
- make createSE returning assays as 'BumpyMatrix' objects (previously 'matrix' objects only)
- add df_to_bm_assay function  (returning BumpyMatrix object from raw experiment data) 
- add initial tests for createSE function

# Change to v.0.1.30 - 2020-12-21
- minor clean-up

# Change to v.0.1.29 - 2020-12-15
- remove dplyr
- use the latest version of gDRutils

# Change to v.0.1.28 - 2020-12-07
- error messages

# Change to v.0.1.27 - 2020-12-02
- refactor analyze_data.R with taking account codilution combo

# Change to v.0.1.26 - 2020-11-24
- bug with reading raw data

# Change to v.0.1.25 - 2020-11-08
- clean-up repository from duplicated functions

# Change to v.0.1.25 - 2020-11-06
- switch to R 4.0

# Change to v.0.1.24 - 2020-10-01
- bugs in function for metrics calculation as per QCS dataset
- QCS data 
- scripts for processing and pushing QCS data on Rosalind

# Change to v.0.1.23 - 2020-09-29
- updated README.md

# Change to v.0.1.22 - 2020-09-22
- bug in update_experiment_metadata

# Change to v.0.1.21 - 2020-09-03
- logic as per new db model

# Change to v.0.1.20 - 2020-08-26
- RDS files for the other GDS data (GDS5, GDS6, GDS7)

# Change to v.0.1.19 - 2020-08-03
- support an import of single agent data from GDS
  * including the masked field to be able to remove the masked data from averages
  * add a function to annotate drug combinations (add_codrug_group)
  * update the tests that contained combination experiments accordingly.
  * add new tests for the GDS data
  * minor refactor to use the 'masked' column. If it is not part of the initial dataframe it is added automatically (and set to FALSE)

# Change to v.0.1.18 - 2020-07-29
- fixed issue in update_experiment_metadata.R
- unit test for update_experiment_metadata.R

# Change to v.0.1.17 - 2020-07-14
- function for updating metadata in SE

# Change to v.0.1.16 - 2020-07-14
- dependencies: kcloak, gDRwrapper

# Change to v.0.1.15 - 2020-07-08
- dependencies
