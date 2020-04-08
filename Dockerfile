# NOTE:
# This Dockerfile is stored under "rplatform" directory but it will be
# moved to top level project directory during image build by RP so
# all local paths should be relative to your project's top level
# directory.
#
# NOTE2:
# Base images with tag 3.5.1_rp0.0.75 and later are built as user 'rstudio'
# which belongs to sudoers, so every command requiring root permissions should be
# preceded with 'sudo' - otherwise add layer 'USER root' to run all commands as root.
# It is recommended to not to change 'rstudio' user due to permissions issues
# within Docker container, because container's RStudio Server is run as 'rstudio'.

FROM registry.rplatform.org:5000/rocker-rstudio-stable:3.5.1_rp0.0.78

# ------ Be aware that any changes in following may cause issue with RPlatform and CBS ---------------------------

LABEL MAINTAINER="scigocki.dariusz@gene.com"
LABEL NAME=gdr_core
LABEL GENERATE_SINGULARITY_IMAGE=false
LABEL production=false
LABEL VERSION=0.0.0.9000
LABEL CACHE_IMAGE="registry.rplatform.org:5000/githubroche/gdrplatform/gdr_core"

# ----------------------------------------------------------------------------------------------------------------

#================= Define your system dependencies in this Dockerfile

# install mysql client
RUN sudo apt-get update && sudo apt-get install -y \
    libmariadb-client-lgpl-dev

## Define your system dependencies in this Dockerfile

COPY gDR/DESCRIPTION /mnt/vol/gDR/DESCRIPTION
COPY rplatform/install_dependencies.R /mnt/vol/rplatform/install_dependencies.R
COPY rplatform/git_dependencies.yml /mnt/vol/rplatform/git_dependencies.yml
RUN R -f /mnt/vol/rplatform/install_dependencies.R

## Uncomment following if installing package(s) from Bitbucket. It requires adding rsa keys to repo if it's not public
COPY rplatform/ssh_keys/id_rsa /home/rstudio/.ssh/id_rsa
COPY rplatform/ssh_keys/id_rsa.pub /home/rstudio/.ssh/id_rsa.pub
RUN sudo chown rstudio:rstudio -R /home/rstudio/.ssh && \
    sudo chmod 400 /home/rstudio/.ssh/id_rsa && \
    eval `ssh-agent -s` && ssh-add /home/rstudio/.ssh/id_rsa

## Uncomment following to install package(s) from source and replace
## `/path/to/package_source_dir` with path to R package directory in your git repository
## `/mnt/vol/package_source_dir` with path to R package directory in docker container (must correspond with
## the one from `install_from_source.R` if using it)
#COPY /path/to/package_source_dir /mnt/vol/package_source_dir

## Uncomment following to install package(s) from source dir or Bitbucket
COPY rplatform/install_from_source.R /mnt/vol/rplatform/install_from_source.R
RUN R -f /mnt/vol/rplatform/install_from_source.R

RUN sudo rm -rf /home/rstudio/.ssh
RUN sudo rm -rf /mnt/vol/* 