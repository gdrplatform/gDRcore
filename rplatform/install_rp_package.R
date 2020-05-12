
## set cores
if (!require(parallel)) install.packages("parallel")

options(Ncpus = parallel::detectCores())

## install pkgs required for 'rp' package 
### new version of git2r causes errors with installing by ssh
MRAN_SNAPSHOT_DATE <- Sys.getenv("MRAN_SNAPSHOT_DATE")
install.packages(c("desc", "devtools", "stringr", "withr", "DelayedMatrixStats"))
devtools::install_version("git2r", version = "0.25.2", repos = "https://cloud.r-project.org")

## install rp directly from bitbucket
ssh_keys <- git2r::cred_ssh_key(file.path("/home/rstudio/.ssh/id_rsa.pub"), file.path("/home/rstudio/.ssh/id_rsa"))

devtools::install_git(
  url = "ssh://git@stash.intranet.roche.com:7999/rp/rp-package.git",
  ref = "master",
  credentials = ssh_keys
)