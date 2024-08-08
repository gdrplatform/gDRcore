## gDRcore 1.3.8 - 2024-08-08
* fix issue with mapping overrides untreated controls

## gDRcore 1.3.7 - 2024-08-05
* fix issue with non-avearaged concentration data

## gDRcore 1.3.6 - 2024-07-23
* fix issue with providing empty nested_confounder

## gDRcore 1.3.5 - 2024-07-17
* allow using custom functions for calculating HSA and Bliss scores for combination data

## gDRcore 1.3.4 - 2024-07-15
* refactor logic for calculating standard deviation for single values

## gDRcore 1.3.3 - 2024-07-10
* update unit tests

## gDRcore 1.3.2 - 2024-06-04
* switch to `get_supported_experiments`

## gDRcore 1.3.1 - 2024-05-27
* synchronize Bioconductor and GitHub versioning

## gDRcore 1.1.22 - 2024-05-16
* move `define_matrix_grid_positions` and `round_concentration` functions to `gDRutils` package

## gDRcore 1.1.21 - 2024-05-13
* add functions for retrieving annotation data from input dt

## gDRcore 1.1.20 - 2024-05-08
* fix typo

## gDRcore 1.1.19 - 2024-04-23
* add vignette with documentation for data annotation

## gDRcore 1.1.18 - 2024-04-22
* add support for external annotation specified in the env var

## gDRcore 1.1.17 - 2024-03-26
* fix bug with converting mae into raw_data

## gDRcore 1.1.16 - 2024-03-15
* remove unstable tests

## gDRcore 1.1.15 - 2024-03-14
* cleanup package

## gDRcore 1.1.14 - 2024-03-12
* update function description

## gDRcore 1.1.13 - 2024-02-26
* improve pkgdown site
  * improved references
  * valid NEWS.md

## gDRcore 1.1.12 - 2024-02-14
* fix issue with retrieving unique records from mix of control and treated samples

## gDRcore 1.1.11 - 2024-02-14
* make documentation compatible with pkdgdown
* rename 'matches' to 'grr_matches'

## gDRcore 1.1.10 - 2024-02-12
* fix unit tests for GRAN

## gDRcore 1.1.9 - 2024-02-07
* simplify logic of assays for combination data
* rename `matrix` into `combination`

## gDRcore 1.1.8 - 2024-02-06
* add support for internal source of annotation

## gDRcore 1.1.7 - 2024-02-06
* fix bug with converting standardize MAE into raw data

## gDRcore 1.1.6 - 2024-02-05
* add vignette for a data model

## gDRcore 1.1.5 - 2024-02-01
* update wrappers for co-dilution data

## gDRcore 1.1.4 - 2024-01-22
* add new description fields

## gDRcore 1.1.3 - 2024-01-04
* improve logic for normalization and identification of single-agent and matrix data with Drug3

## gDRcore 1.1.2 - 2023-12-15
* fix issue with wrong assignment of `untreated` records

## gDRcore 1.1.1 - 2023-11-22
* sync master with devel branch
* add support for unifying duplicates in combo matrix data
* add "Treatment" as template identifier

## gDRcore 1.1.0 - 2023-10-24
* release Bioc 3.18

## gDRcore 1.0.0 - 2023-10-24
* prerelease Bioc 3.18

## gDRcore 0.99.43 - 2023-10-17
* adjust NEWS to Bioc format

## gDRcore 0.99.42 - 2023-10-05
* bump version of gDRtestData
* fix bug with merging controls in triple combo with additional perturbations

## gDRcore 0.99.41 - 2023-09-25
* add support for adding custom annotations inside input files
* improve the performance

## gDRcore 0.99.40 - 2023-09-25
* fix bug with subsetting wrong combo matrix value

## gDRcore 0.99.39 - 2023-09-19
* extend the logic for matching missing controls

## gDRcore 0.99.38 - 2023-09-12
* set Drug3 as an official tertiary drug in the experiment

## gDRcore 0.99.37 - 2023-09-04
* fill NA by average values when there is no match with plate

## gDRcore 0.99.36 - 2023-09-01
* fill NA during aggregation of ref and trt data with mean

## gDRcore 0.99.35 - 2023-08-17
* fix issue with missing subsetting Day0 data

## gDRcore 0.99.34 - 2023-08-16
* update logic for supporting manifest and template files sharing the same column

## gDRcore 0.99.33 - 2023-08-10
* update annotation column names for cell line annotation as per changes in the gDRutils

## gDRcore 0.99.32 - 2023-07-25
* extend logic for supporting cols with dash, e.g. additional perturbations with "-"

## gDRcore 0.99.31 - 2023-07-19
* update the logic for handling warnings in the pipeline

## gDRcore 0.99.30 - 2023-07-13
* fix issue with wrong merging of data.tables without nested confounders

## gDRcore 0.99.29 - 2023-07-07
* add information about source type for cases without metric data
* refactor the logic for splitting raw data from metadata - get rid of iterative approach

## gDRcore 0.99.28 - 2023-07-05
* update logic for parallel computing

## gDRcore 0.99.27 - 2023-06-29
* optimize unit tests

## gDRcore 0.99.26 - 2023-06-29
* remove backward compatibility for old data model

## gDRcore 0.99.25 - 2023-06-27
* fix bug with missing rownames in normalized assay

## gDRcore 0.99.24 - 2023-06-19
* update logic for merging data.table objects

## gDRcore 0.99.23 - 2023-06-13
* replace `order` with `data.table::setorder`

## gDRcore 0.99.22 - 2023-06-09
* switch from `merge` to `[[` for `data.table` objects

## gDRcore 0.99.21 - 2023-06-06
* switch from `aggregate` to `data.table`

## gDRcore 0.99.20 - 2023-06-07
* switch from `zoo::rollmean` to `data.table::frollmean`

## gDRcore 0.99.19 - 2023-06-06
* replace reshape2 functions by functions from data.table

## gDRcore 0.99.18 - 2023-05-31
* fix managing of mixed types of raw data

## gDRcore 0.99.17 - 2023-05-29
* fix bug with subsetting data for calculating isobologram

## gDRcore 0.99.16 - 2023-05-22
* format the vignette with BiocStyle

## gDRcore 0.99.15 - 2023-05-16
* fix related with data.table

## gDRcore 0.99.14 - 2023-05-15
* rename `excess` to `x` to unify colnames in assay data

## gDRcore 0.99.13 - 2023-05-10
* refactor normalization_types in combo-specific assays

## gDRcore 0.99.12 - 2023-05-09
* utilize `gDRutils::apply_bumpy_function` in fit_SE

## gDRcore 0.99.11 - 2023-05-05
* fix bug with swapping untreated/vehicle values

## gDRcore 0.99.10 - 2023-05-04
* fix bug with data.table

## gDRcore 0.99.9 - 2023-04-21
* utilize `gDRutils::apply_bumpy_function` in average_SE

## gDRcore 0.99.8 - 2023-04-20
* switch to OSI license

## gDRcore 0.99.7 - 2023-04-19
* fix bug with replacing vehicle to untreated values

## gDRcore 0.99.6 - 2023-04-19
* moved wrapper fuctions from gDRtestData

## gDRcore 0.99.5 - 2023-04-18
* update dependencies
* add fix for bioc-devel - correct sorting in merge test

## gDRcore 0.99.4 - 2023-04-17
* fix namespacing issue in examples
* add R 4.2 as dependency
* fix examples for normalize_SE

## gDRcore 0.99.3 - 2023-04-12
* add logic for retrieving raw data from assay data

## gDRcore 0.99.2 - 2023-04-07
* update maintainer

## gDRcore 0.99.1 - 2023-04-04
* bugfix for the logic in 'cleanup_metadata'

## gDRcore 0.99.0 - 2023-03-24
* make the package Bioc-compatible

## gDRcore 0.1.3.22 - 2023-03-21
* improve performance of 'map_df' with refactored logic for exact matches

## gDRcore 0.1.3.21 - 2023-03-15
* refactor pipeline

## gDRcore 0.1.3.20 - 2023-03-09
* address co-treatment fit by using the matrix data type instead

## gDRcore 0.1.3.19 - 2023-03-07
* add support for splitting normalization types

## gDRcore 0.1.3.18 - 2023-03-06
* improve logic in functions used to generate isobolograms' data

## gDRcore 0.1.3.17 - 2023-03-06
* remove obsolete code

## gDRcore 0.1.3.16 - 2023-02-21
* add support for partial pipeline runs

## gDRcore 0.1.3.15 - 2023-02-10
* update path to annotation data

## gDRcore 0.1.3.14 - 2023-02-10
* clean-up code

## gDRcore 0.1.3.13 - 2023-01-13
* clean-up code

## gDRcore 0.1.3.12 - 2023-01-12
* refactor mapping function to properly handle drug3

## gDRcore 0.1.3.11 - 2022-12-16
* replace parallelize function with gDRutils::loop

## gDRcore 0.1.3.10 - 2022-12-15
* add assert for vehicle values in input data in `runDrugResponseProcessingPipeline`

## gDRcore 0.1.3.9 - 2022-12-14
* fix error-handling if conditions in average_SE

## gDRcore 0.1.3.8 - 2022-11-07
* fix invalid encapsulation in tests

## gDRcore 0.1.3.7 - 2022-10-18
* add missing namespacing

## gDRcore 0.1.3.6 - 2022-10-10
* add support for setting number of cores for `BiocParallel` based on the env variable

## gDRcore 0.1.3.5 - 2022-10-07
* remove global parameters for `BiocParallel`

## gDRcore 0.1.3.4 - 2022-09-29
* change the logic for using cores in `BiocParallel`

## gDRcore 0.1.3.3 - 2022-09-27
* update the logic for parallel computing

## gDRcore 0.1.3.2 - 2022-08-17
* update the logic for default nested_confounders in `create_SE` function

## gDRcore 0.1.3.1 - 2022-07-08
* refactor create_SE to support reverse single-agent data

## gDRcore 0.1.3.0 - 2022-06-02
* release 1.3.0

## gDRcore 0.1.1.39 - 2022-05-30
* add missing namespace for get_env_identifiers

## gDRcore 0.1.1.38 - 2022-05-09
* replace `NA` by 0 in Concentration loaded in manifest file

## gDRcore 0.1.1.37 - 2022-05-06
* switch from `data.table` to `data.frame` in add_annotation* functions

## gDRcore 0.1.1.36 - 2022-04-12
* update function for adding unknown cell line annotations
* update logic for using nested confounders
* remove `grr` from dependencies

## gDRcore 0.1.1.35 - 2022-04-07
* get rid of gDRwrapper

## gDRcore 0.1.1.34 - 2022-03-30
* add support for additional barcode identifiers

## gDRcore 0.1.1.33 - 2022-03-18
* fix documentation in calculate_matrix_metric

## gDRcore 0.1.1.32 - 2022-03-17
* switch from `catchr` to `purrr`

## gDRcore 0.1.1.31 - 2022-03-02
* fix mapping reference values of inverted treatments
* remove duplication of single-agent data

## gDRcore 0.1.1.30 - 2022-02-17
* wrap SE into MAE on the level of runDrugResponseProcessingPipeline

## gDRcore 0.1.1.29 - 2022-02-14
* change identifier `drugname` to `drug_name`

## gDRcore 0.1.1.28 - 2022-02-09
* switch from `gDRinternal` to `gDRinternalData` for internal annotations

## gDRcore 0.1.1.27 - 2022-02-08
* issue with subsetting by list in R 4.2.0

## gDRcore 0.1.1.26 - 2022-02-02
* unlist output of intersect as per R 4.2.0

## gDRcore 0.1.1.25 - 2022-01-31
* align version criteria between `dependencies.yaml` and `DESCRIPTION` package versions

## gDRcore 0.1.1.24 - 2022-01-25
* standardize/improve CI

## gDRcore 0.1.1.23 - 2022-01-24
* switch from `cores` variable to `detect_cores- ` function

## gDRcore 0.1.1.22 - 2022-01-24
* fix wrong type of `NUM_CORES` env variable

## gDRcore 0.1.1.21 - 2022-01-03
* speed-up functions for mapping treated and untreated cases

## gDRcore 0.1.1.20 - 2021-12-30
* fix linter issues

## gDRcore 0.1.1.19 - 2021-12-14
* use parallel computing as an alternative for `for` loops

## gDRcore 0.1.1.18 - 2021-12-07
* update annotation script as per new csv annotation files

## gDRcore 0.1.1.17 - 2021-12-07
* detect co-trt data and treat them as single-agent

## gDRcore 0.1.1.16 - 2021-11-08
* set excess = NA for single-agent

## gDRcore 0.1.1.15 - 2021-10-29
* solve rounding issues
* add new bliss metric

## gDRcore 0.1.1.14 - 2021-10-25
* refactor isobolograms

## gDRcore 0.1.1.13 - 2021-10-25
* move p_trt_keys to the proper place

## gDRcore 0.1.1.12 - 2021-10-20
* add support for masked data in fit_SE.combinations.R

## gDRcore 0.1.1.11 - 2021-10-14
* address issues in creating the SE for combo matrix experiments

## gDRcore 0.1.1.10 - 2021-10-13
* refactor the logic for calculation Loewe when there is no `Concentration == 0`

## gDRcore 0.1.1.9 - 2021-10-12
* refactor the logic for combo data

## gDRcore 0.1.1.8 - 2021-09-27
* updated normalization_types in 'calculate_combo_matrix' and 'fit_SE'

## gDRcore 0.1.1.7 - 2021-09-21
* calculate_GR_value by removing cl_name param

## gDRcore 0.1.1.6 - 2021-08-25
* fix but with missing `nested_identifiers` variables for creating DataFrame for masked values

## gDRcore 0.1.1.5 - 2021-08-04
* refactor reading annotations and add default parameters

## gDRcore 0.1.1.4 - 2021-07-30
* fix bug for nested_confounders not present in ref_df

## gDRcore 0.1.1.3 - 2021-07-27
* add support for nested_identifiers and nested_confounders 

## gDRcore 0.1.1.2 - 2021-07-23
* remove obsolete dependencies from DESCRIPTION:Imports

## gDRcore 0.1.0.9 - 2021-07-01
* function for testing synthetic data

## gDRcore 0.1.0.8 - 2021-06-25
* move functions for importing template files to gDRimport package 

## gDRcore 0.1.0.7 - 2021-06-23
* add linter
* remove obsolete functions from test.utils.R

## gDRcore 0.1.0.6 - 2021-06-22
* move importing functions to gDRimport package

## gDRcore 0.1.0.5 - 2021-06-14
* add unit tests for synthetic data

## gDRcore 0.1.0.4 - 2021-06-14
* remove deprecated functions and unit tests

## gDRcore 0.1.0.3 - 2021-06-14
* update namespace for 'metadata- '

## gDRcore 0.1.0.2 - 2021-06-10
* change package name - gDR => gDRcore

## gDRcore 0.1.0.1 - 2021-06-04
* export/update docs for 'standardize_record_values'

## gDRcore 0.1.0.0 - 2021-06-02
* release 1.0.0

## gDRcore 0.0.1.52 - 2021-04-29
* remove all the gDRinternal-related files

## gDRcore 0.0.1.51 - 2021-04-27
* fix wrong argument name in `fit_curves`

## gDRcore 0.0.1.50 - 2021-04-23
* ensure that dts in `assay- se,"Averaged"` are NULL when there are not treatments

## gDRcore 0.0.1.49 - 2021-04-20
* add processing info metadata to SE

## gDRcore 0.0.1.48 - 2021-04-14
* improve handling of nested_keys and override_controls

## gDRcore 0.0.1.47 - 2021-04-09
* get rid of sorting index columns in getMetaData

## gDRcore 0.0.1.46 - 2021-03-31
* move SE metadata getters and setters to gDRutils

## gDRcore 0.0.1.45 - 2021-03-31
* fix a bug in map_df function

## gDRcore 0.0.1.44 - 2021-03-26
* refactor processing functions in create_SE2, normalize_SE2, average_SE2, and fit_SE2 and friends

## gDRcore 0.0.1.43 - 2021-02-22
* add drug_moa_2 to merged df

## gDRcore 0.0.1.42 - 2021-02-17
* fix a bug with adding annotation to codrugs

## gDRcore 0.0.1.41 - 2021-02-04
* add script for benchmarking [normalize/average/metrics]-SE functions

## gDRcore 0.0.1.40 - 2021-02-02
* fix unit tests
* include discard keys info in getMetaData for rowData

## gDRcore 0.0.1.39 - 2021-01-29
* utilize refactored and renamed fitting function RVGRfit to fit_curves

## gDRcore 0.0.1.38 - 2021-01-27
* export some functions for use in gDRutils

## gDRcore 0.0.1.37 - 2021-01-20
* move df_to_assay and df_to_bm_assay to gDRutils

## gDRcore 0.0.1.36 - 2021-01-18
* updated rowData and colData in SE with additional columns

## gDRcore 0.0.1.35 - 2021-01-14
* update positional header calls to explicit identifier calls
* minor changes to comply with gDRstyle
* minor changes to use identifiers over hard coded DrugName requirements

## gDRcore 0.0.1.34 - 2021-01-12
* update scripts related to annotation of drugs and cell lines

## gDRcore 0.0.1.33 - 2021-01-05
* anonymize test data

## gDRcore 0.0.1.32 - 2020-12-24
* switch from GeneDataScreenR to gDRinternal from QCS to gDR functions

## gDRcore 0.0.1.31 - 2020-12-22
* make createSE returning assays as 'BumpyMatrix' objects - previously 'matrix' objects only
* add df_to_bm_assay function  - returning BumpyMatrix object from raw experiment data 
* add initial tests for createSE function

## gDRcore 0.0.1.30 - 2020-12-21
* minor clean-up

## gDRcore 0.0.1.29 - 2020-12-15
* remove dplyr
* use the latest version of gDRutils

## gDRcore 0.0.1.28 - 2020-12-07
* error messages

## gDRcore 0.0.1.27 - 2020-12-02
* refactor analyze_data.R with taking account codilution combo

## gDRcore 0.0.1.26 - 2020-11-24
* bug with reading raw data

## gDRcore 0.0.1.25 - 2020-11-08
* clean-up repository from duplicated functions

## gDRcore 0.0.1.25 - 2020-11-06
* switch to R 4.0

## gDRcore 0.0.1.24 - 2020-10-01
* bugs in function for metrics calculation as per QCS dataset
* QCS data 
* scripts for processing and pushing QCS data on Rosalind

## gDRcore 0.0.1.23 - 2020-09-29
* updated README.md

## gDRcore 0.0.1.22 - 2020-09-22
* bug in update_experiment_metadata

## gDRcore 0.0.1.21 - 2020-09-03
* logic as per new db model

## gDRcore 0.0.1.20 - 2020-08-26
* RDS files for the other GDS data - GDS5, GDS6, GDS7

## gDRcore 0.0.1.19 - 2020-08-03
* support an import of single agent data from GDS
  * including the masked field to be able to remove the masked data from averages
  * add a function to annotate drug combinations - add_codrug_group
  * update the tests that contained combination experiments accordingly.
  * add new tests for the GDS data
  * minor refactor to use the 'masked' column. If it is not part of the initial dataframe it is added automatically - and set to FALSE

## gDRcore 0.0.1.18 - 2020-07-29
* fixed issue in update_experiment_metadata.R
* unit test for update_experiment_metadata.R

## gDRcore 0.0.1.17 - 2020-07-14
* function for updating metadata in SE

## gDRcore 0.0.1.16 - 2020-07-14
* dependencies: kcloak, gDRwrapper

## gDRcore 0.0.1.15 - 2020-07-08
* dependencies
