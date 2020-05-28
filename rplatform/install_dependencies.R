# RP package template version >= 0.0.77

# The base image makes use of the stable GRAN repository corresponding to image's R version as a source of packages. 
# If you wish to install custom package versions from other source you'll need to
# modify the "repos" option or provide the repository url explicitly during package installation.

# To install and verify packages use rp::installAndVerify function.
# Version requirements can be defined with "requirement" argument:
# -- install = install.packages
# '*' - any version
# '==0.1' - package version equal to 0.1
# '>=0.1' - package version greater than or equal to 0.1
# -- install = install_github
# 'r:tag_name' - reference github tag
# 's:sha1' - github SHA1

# Install packages using multiple cores
if (!require(parallel)) rp::installAndVerify(package = "parallel", requirement = "*")
options(Ncpus = parallel::detectCores())

.wd <- "/mnt/vol"

# don't install these packages - they will be installed separately
git_pkgs <- yaml::read_yaml(file.path(.wd, "rplatform", "git_dependencies.yml"))
dont.install <- c(
  names(git_pkgs$pkgs)
) 

# Extract dependencies from DESCRIPTION file
deps <- yaml::read_yaml(file.path(.wd, "rplatform", "DESCRIPTION_dependencies.yaml"))
deps <- deps[!names(deps) %in% dont.install]

# packages needed in gDRshiny
deps <- c(deps, 
          "DT" = "*", 
          "shiny" = "*", 
          "shinyjs" = "*", 
          "shinyBS" = "*", 
          "shinyalert" = "*", 
          "shinydashboard" = "*", 
          "plotly" = "*")

### new version of mclust (eg. 5.4.6) causes errors:
####  *** caught illegal operation ***
#### address 0x7f3bedee8c18, cause 'illegal operand'
#devtools::install_version("mclust", version = "4.0", repos = "https://cloud.r-project.org")                               

rp::installAndVerify(install = install.packages,
                     package = names(deps),
                     requirement = deps)


