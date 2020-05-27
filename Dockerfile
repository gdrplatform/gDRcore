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

FROM registry.rplatform.org:5000/rocker-rstudio-uat:3.6.1_rp0.0.79

# ------ Be aware that any changes in following may cause issue with RPlatform and CBS ---------------------------

LABEL MAINTAINER="scigocki.dariusz@gene.com"
LABEL NAME=gdr_core
LABEL GENERATE_SINGULARITY_IMAGE=false
LABEL production=false
LABEL VERSION=0.0.0.9100
LABEL CACHE_IMAGE="registry.rplatform.org:5000/githubroche/gdrplatform/gdr_core"
ARG MRAN_SNAPSHOT_DATE="2019-12-12"

# ----------------------------------------------------------------------------------------------------------------

# install system dependencies
RUN sudo apt-get update && sudo apt-get install -y \
    libmariadb-client-lgpl-dev \
    libmariadbclient-dev 

#================= copy Rprofile.site - set repos and other options
COPY rplatform/Rprofile.site /tmp/Rprofile.site
RUN sudo echo 'Sys.setenv(MRAN_SNAPSHOT_DATE = "'$MRAN_SNAPSHOT_DATE'")' "$(cat /tmp/Rprofile.site)" > /tmp/Rprofile.site
RUN mv /tmp/Rprofile.site $(R RHOME)/etc/Rprofile.site

#================= copy ssh keys
COPY rplatform/ssh_keys/id_rsa /home/rstudio/.ssh/id_rsa
COPY rplatform/ssh_keys/id_rsa.pub /home/rstudio/.ssh/id_rsa.pub

#================= install rp R package
COPY rplatform/install_rp_package.R /mnt/vol/rplatform/install_rp_package.R
RUN R -f /mnt/vol/rplatform/install_rp_package.R

COPY rplatform/DESCRIPTION_dependencies.yaml /mnt/vol/rplatform/DESCRIPTION_dependencies.yaml
COPY gDR/DESCRIPTION /mnt/vol/gDR/DESCRIPTION
COPY rplatform/install_dependencies.R /mnt/vol/rplatform/install_dependencies.R
COPY rplatform/git_dependencies.yml /mnt/vol/rplatform/git_dependencies.yml
RUN R -f /mnt/vol/rplatform/install_dependencies.R

## Uncomment following to install package(s) from source dir or Bitbucket
COPY rplatform/install_from_source.R /mnt/vol/rplatform/install_from_source.R
RUN R -f /mnt/vol/rplatform/install_from_source.R

#============ Disable login requirement for Rstudio
ENV DISABLE_AUTH=true

RUN sudo rm -rf /home/rstudio/.ssh
RUN sudo rm -rf /mnt/vol/* 
