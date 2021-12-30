<h3 align = "left"><strong>Changelog</strong></h3>

All notable changes to this project will be documented in this file.

#### [1.1.20] - 2021-12-30
#### Update
- fix linter issues

#### [1.1.19] - 2021-12-14
#### Update
- use parallel computing as an alternative for `for` loops

#### [1.1.18] - 2021-12-07
#### Update
- update annotation script as per new csv annotation files

#### [1.1.17] - 2021-12-07
#### Update
- detect co-trt data and treat them as single-agent

#### [1.1.16] - 2021-11-08
#### Update
- set excess = NA for single-agent

#### [1.1.15] - 2021-10-29
#### Update
- solve rounding issues
- add new bliss metric

#### [1.1.14] - 2021-10-25
#### Update
- refactor isobolograms

#### [1.1.13] - 2021-10-25
#### Update
- move p_trt_keys to the proper place

#### [1.1.12] - 2021-10-20
#### Update
- add support for masked data in fit_SE.combinations.R

#### [1.1.11] - 2021-10-14
#### Update
- address issues in creating the SE for combo matrix experiments

#### [1.1.10] - 2021-10-13
#### Update
- refactor the logic for calculation Loewe when there is no `Concentration == 0`

#### [1.1.9] - 2021-10-12
#### Update
- refactor the logic for combo data

#### [1.1.8] - 2021-09-27
#### Update
- updated normalization_types in 'calculate_combo_matrix' and 'fit_SE'

#### [1.1.7] - 2021-09-21
#### Update
- calculate_GR_value by removing cl_name param

#### [1.1.6] - 2021-08-25
#### Update
- fix but with missing `nested_identifiers` variables for creating DataFrame for masked values

#### [1.1.5] - 2021-08-04
#### Update
- refactor reading annotations and add default parameters

#### [1.1.4] - 2021-07-30
#### Update
- fix bug for nested_confounders not present in ref_df

#### [1.1.3] - 2021-07-27
#### Update
- add support for nested_identifiers and nested_confounders 

#### [1.1.2] - 2021-07-23
#### Update
- remove obsolete dependencies from DESCRIPTION:Imports

#### [1.0.9] - 2021-07-01
#### Update
- function for testing synthetic data

#### [1.0.8] - 2021-06-25
#### Update
- move functions for importing template files to gDRimport package 

#### [1.0.7] - 2021-06-23
#### Update
- add linter
- remove obsolete functions from test.utils.R

#### [1.0.6] - 2021-06-22
#### Update
- move importing functions to gDRimport package

#### [1.0.5] - 2021-06-14
#### Update
- add unit tests for synthetic data

#### [1.0.4] - 2021-06-14
#### Update
- remove deprecated functions and unit tests

#### [1.0.3] - 2021-06-14
#### Update
- update namespace for 'metadata()'

#### [1.0.2] - 2021-06-10
#### Update
- change package name (gDR => gDRcore)

#### [1.0.1] - 2021-06-04
#### Update
- export/update docs for 'standardize_record_values'

#### [1.0.0] - 2021-06-02
#### Update
- release 1.0.0

#### [0.1.52] - 2021-04-29
#### Update
- remove all the gDRinternal-related files

#### [0.1.51] - 2021-04-27
#### Update
- fix wrong argument name in `fit_curves`

#### [0.1.50] - 2021-04-23
#### Update
- Ensure that dts in `assay(se,"Averaged")` are NULL when there are not treatments

#### [0.1.49] - 2021-04-20
#### Update
- add processing info metadata to SE

#### [0.1.48] - 2021-04-14
#### Update
- improve handling of nested_keys and override_controls

#### [0.1.47] - 2021-04-09
#### Update
- get rid of sorting index columns in getMetaData

#### [0.1.46] - 2021-03-31
#### Update
- move SE metadata getters and setters to gDRutils

#### [0.1.45] - 2021-03-31
#### Update
- fix a bug in map_df function


#### [0.1.44] - 2021-03-26
#### Update
- refactor processing functions in create_SE2, normalize_SE2, average_SE2, and fit_SE2 and friends

#### [0.1.43] - 2021-02-22
#### Update
- add drug_moa_2 to merged df

#### [0.1.42] - 2021-02-17
#### Update
- fix a bug with adding annotation to codrugs

#### [0.1.41] - 2021-02-04
#### Update
- Add script for benchmarking [normalize/average/metrics]-SE functions

#### [0.1.40] - 2021-02-02
#### Update
- fix unit tests
- include discard keys info in getMetaData for rowData

#### [0.1.39] - 2021-01-29
#### Update
- utilize refactored and renamed fitting function RVGRfit to fit_curves

#### [0.1.38] - 2021-01-27
#### Update
- export some functions for use in gDRutils

#### [0.1.37] - 2021-01-20
#### Update
- move df_to_assay and df_to_bm_assay to gDRutils

#### [0.1.36] - 2021-01-18
#### Update
- updated rowData and colData in SE with additional columns

#### [0.1.35] - 2021-01-14
#### Update
- update positional header calls to explicit identifier calls
- minor changes to comply with gDRstyle
- minor changes to use identifiers over hard coded DrugName requirements

#### [0.1.34] - 2021-01-12
#### Update
- update scripts related to annotation of drugs and cell lines

#### [0.1.33] - 2021-01-05
#### Update
- anonymize test data

#### [0.1.32] - 2020-12-24
#### Update
- switch from GeneDataScreenR to gDRinternal from QCS to gDR functions

#### [0.1.31] - 2020-12-22
#### Added
- make createSE returning assays as 'BumpyMatrix' objects (previously 'matrix' objects only)
- add df_to_bm_assay function  (returning BumpyMatrix object from raw experiment data) 
- add initial tests for createSE function

#### [0.1.30] - 2020-12-21
#### Improved
- minor clean-up

#### [0.1.29] - 2020-12-15
#### Improved
- remove dplyr
- use the latest version of gDRutils

#### [0.1.28] - 2020-12-07
#### Improved
- error messages

#### [0.1.27] - 2020-12-02
#### Improved
- refactor analyze_data.R with taking account codilution combo

#### [0.1.26] - 2020-11-24

##### Fixed
- bug with reading raw data


#### [0.1.25] - 2020-11-08

##### Improved
- clean-up repository from duplicated functions

#### [0.1.25] - 2020-11-06

##### Update
- switch to R 4.0

#### [0.1.24] - 2020-10-01

##### Update
- bugs in function for metrics calculation as per QCS dataset

#### Added
- QCS data 
- scripts for processing and pushing QCS data on Rosalind

#### [0.1.23] - 2020-09-29

##### Improved
- updated README.md


#### [0.1.22] - 2020-09-22

##### Fixed
- bug in update_experiment_metadata

#### [0.1.21] - 2020-09-03

##### Update
- logic as per new db model

#### [0.1.20] - 2020-08-26

##### Added
- RDS files for the other GDS data (GDS5, GDS6, GDS7)

#### [0.1.19] - 2020-08-03

##### Updated
- support an import of single agent data from GDS
  * including the masked field to be able to remove the masked data from averages
  * add a function to annotate drug combinations (add_codrug_group)
  * update the tests that contained combination experiments accordingly.
  * add new tests for the GDS data
  * minor refactor to use the 'masked' column. If it is not part of the initial dataframe it is added automatically (and set to FALSE)

#### [0.1.18] - 2020-07-29

##### Updated
- fixed issue in update_experiment_metadata.R
##### Added
- unit test for update_experiment_metadata.R

#### [0.1.17] - 2020-07-14

##### Added
- function for updating metadata in SE

#### [0.1.16] - 2020-07-14

##### Added
- dependencies: kcloak, gDRwrapper

#### [0.1.15] - 2020-07-08

##### Updated
- dependencies
