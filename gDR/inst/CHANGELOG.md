<h3 align = "left"><strong>Changelog</strong></h3>

All notable changes to this project will be documented in this file.
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
