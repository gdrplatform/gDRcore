# RP package template version >= 0.0.77

## Uncomment following code to install package from source directory
# rp::installAndVerify(
#   install = devtools::install,
#   requirement = sprintf("== %s", desc::desc_get_version("/mnt/vol/package_source_dir")),
#   package = "/mnt/vol/package_source_dir"
# )

## Uncomment following code to install package(s) directly from Bitbucket repository
## SSH keys should be copied to container before installation
ssh_keys <- git2r::cred_ssh_key(file.path("/home/rstudio/.ssh/id_rsa.pub"), file.path("/home/rstudio/.ssh/id_rsa"))
.wd <- "/mnt/vol"
.deps <- rp:::collectDependencies(desc.files = file.path(.wd, "gDR/DESCRIPTION"))
pkgs <- yaml::read_yaml(file.path(.wd, "rplatform", "git_dependencies.yml"))$pkgs

for (nm in names(pkgs))
  rp::installAndVerify(
    install = devtools::install_git,
    url = pkgs[[nm]]$url,
    ref = pkgs[[nm]]$ref,
    credentials = .ssh_keys,
    package = nm,
    # version requirement is taken from DESCRIPTION if not specified manually in yaml
    requirement = if (!is.null(pkgs[[nm]][["ver"]])) pkgs[[nm]][["ver"]] else .deps[[nm]], 
    subdir = pkgs[[nm]]$subdir,
    upgrade = "never"
  )

