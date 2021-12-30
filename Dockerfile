FROM dockerreg.kubemeainfra.science.roche.com/cbs/githubroche/gdrplatform/gdr_shiny:master-c274a040-build-8

# ------ Be aware that any changes in following may cause issue with RPlatform and CBS

LABEL MAINTAINER="scigocki.dariusz@gene.com"
LABEL NAME=gdr_core
LABEL GENERATE_SINGULARITY_IMAGE=false
LABEL production=false
LABEL VERSION=0.0.0.9100
LABEL CACHE_IMAGE="registry.rplatform.org:5000/githubroche/gdrplatform/gdr_core"
ARG MRAN_SNAPSHOT_DATE="2020-06-01"

# ----------------------------------------------------------------------------------------------------------------

#================= Install dependencies
COPY rplatform/dependencies.yaml rplatform/.github_access_token.txt* /mnt/vol/
COPY rplatform/install_all_deps.R /mnt/vol/install_all_deps.R
RUN R -f /mnt/vol/install_all_deps.R

#================= Build package
COPY gDRcore/ /tmp/gDRcore/
COPY rplatform/install_repo.R /mnt/vol/
RUN R -f /mnt/vol/install_repo.R 

#================= Clean up
RUN sudo rm -rf /mnt/vol/* /tmp/gDRcore/
