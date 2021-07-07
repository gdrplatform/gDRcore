# Settings
repos <- c(CRAN = "https://cran.r-project.org")
options(repos = repos)
essential_pkgs <- list(
  list(name = "git2r", version = "0.28.0"),
  list(name = "yaml", version = "2.2.1"),
  list(name = "BiocManager", version = "1.30.16")
)
deps_yaml <- "/mnt/vol/dependencies.yaml"
use_ssh <- FALSE
# ssh_key_pub <- "/home/rstudio/.ssh/id_rsa.pub"
# ssh_key_priv <- "/home/rstudio/.ssh/id_rsa"

# Auxiliary functions
verify_version <- function(name, required_version) {
  pkg_version <- packageVersion(name)
  ## '>=1.2.3' => '>= 1.2.3'
  required_version <-
    gsub("^([><=]+)([0-9.]+)$", "\\1 \\2", required_version, perl = TRUE)
  if (!remotes:::version_satisfies_criteria(pkg_version, required_version)) {
    stop(sprintf(
      "Invalid version of %s. Installed: %s, required %s.",
      name,
      pkg_version,
      required_version
    ))
  }
}

# Install {remotes}
if (!"remotes" %in% installed.packages()) {
  install.packages(pkgs = "remotes")
}

# Install essential tools
for (pkg in essential_pkgs) {
  if (!pkg$name %in% installed.packages()) {
    remotes::install_version(
      package = pkg$name,
      version = pkg$version
    )
  }
}

# Use SSH keys
keys <- if (isTRUE(use_ssh)) {
  git2r::cred_ssh_key(
    publickey = ssh_key_pub,
    privatekey = ssh_key_priv
  )
}

# Install all dependencies
deps <- yaml::read_yaml(deps_yaml)$pkgs
for (name in names(deps)) {
  pkg <- deps[[name]]
  if (is.null(pkg$source)) { pkg$source <- "Git" }
  switch(toupper(pkg$source),

    ## CRAN installation
    "CRAN" = {
      if (is.null(pkg$repos)) { pkg$repos <- repos }
      remotes::install_version(
        package = name,
        version = pkg$ver,
        repos = pkg$repos
      )
    },

    ## Bioconductor installation
    "BIOC" = {
      if (is.null(pkg$ver)) { pkg$ver <- BiocManager::version() }
      BiocManager::install(
        pkgs = name,
        update = FALSE,
        version = pkg$ver  ## Bioc version or 'devel'
      )
    },
    
    ## GitHub installation
    "GITHUB" = {
      if (is.null(pkg$ref)) { pkg$ref <- "HEAD" }
      remotes::install_github(
        repo = pkg$url,
        ref = pkg$ref,
        subdir = pkg$subdir
      )
      verify_version(name, pkg$ver)
    },
    
    ## Git installation
    "GIT" = {
      remotes::install_git(
        url = pkg$url,
        subdir = pkg$subdir,
        ref = pkg$ref,
        credentials = keys
      )
      verify_version(name, pkg$ver)
    }
  )
}
