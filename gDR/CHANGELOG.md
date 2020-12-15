<h3 align = "left"><strong>Changelog</strong></h3>

All notable changes to this project will be documented in this file.
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
