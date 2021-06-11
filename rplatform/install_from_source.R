# RP package template version >= 0.0.77

## Uncomment following code to install package from source directory
# rp::installAndVerify(
#   install = devtools::install,
#   requirement = sprintf("== %s", desc::desc_get_version("/mnt/vol/package_source_dir")),
#   package = "/mnt/vol/package_source_dir"
# )

install.packages("git2r", repos = paste0("https://mran.microsoft.com/snapshot/", Sys.getenv("MRAN_SNAPSHOT_DATE")))

## Uncomment following code to install package(s) directly from Bitbucket repository
## SSH keys should be copied to container before installation
ssh_keys <- git2r::cred_ssh_key(file.path("/home/rstudio/.ssh/id_rsa.pub"), file.path("/home/rstudio/.ssh/id_rsa"))
.wd <- "/mnt/vol"
.deps <- rp:::collectDependencies(desc.files = file.path(.wd, "gDRcore/DESCRIPTION"))
pkgs <- yaml::read_yaml(file.path(.wd, "rplatform", "git_dependencies.yml"))$pkgs


git2r::libgit2_features()

for (nm in names(pkgs)){
  
  #TODO: fix/remove this temporary fix
  if (!nm %in% c("gDRutils", "gDRwrapper")){
    rp::installAndVerify(
      install = devtools::install_git,
      url = pkgs[[nm]]$url,
      ref = pkgs[[nm]]$ref,
      credentials = ssh_keys,
      package = nm,
      # version requirement is taken from DESCRIPTION if not specified manually in yaml
      requirement = if (!is.null(pkgs[[nm]][["ver"]])) pkgs[[nm]][["ver"]] else .deps[[nm]],
      subdir = pkgs[[nm]]$subdir,
      upgrade = "never"
    )
  } else {
    #TODO: fix/remove it
    repo_path <- file.path(tempdir(), nm)
    git2r::clone(pkgs[[nm]]$url, local_path = repo_path)
    git2r::checkout(object = repo_path, branch = pkgs[[nm]]$ref)
    devtools::install(file.path(repo_path, nm), upgrade = "never")
  }  
}
