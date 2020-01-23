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

# pkgs_to_install <- c(
#  # Add your dependencies here
# )
# rp::installAndVerify(package = pkgs_to_install)

.wd <- "/mnt/vol"

# don't install these packages - they will be installed separately
git_pkgs <- yaml::read_yaml(file.path(.wd, "rplatform", "git_dependencies.yml"))
dont.install <- c(
  names(git_pkgs$pkgs)
) 

# Extract dependencies from DESCRIPTION file
deps <- desc::desc_get_deps(file.path(.wd, "gDR/DESCRIPTION"))
deps <- deps[!(deps$package %in% dont.install), ]
# packages needed in gDRshiny
deps <- rbind(deps, data.frame(type = rep("Suggests", 6),
                               package = c("DT", "shiny", "shinyjs", "shinyBS", "shinyalert","shinydashboard"),
                               version = rep("*", 6)))

rp::installAndVerify(install = install.packages,
                     package = deps$package,
                     requirement = deps$version)


