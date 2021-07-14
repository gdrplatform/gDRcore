FROM rocker/rstudio:4.1.0

# ------ Be aware that any changes in following may cause issue with RPlatform and CBS

LABEL MAINTAINER="scigocki.dariusz@gene.com"
LABEL NAME=gdr_core
LABEL GENERATE_SINGULARITY_IMAGE=false
LABEL production=false
LABEL VERSION=0.0.0.9100
LABEL CACHE_IMAGE="registry.rplatform.org:5000/githubroche/gdrplatform/gdr_core"
ARG MRAN_SNAPSHOT_DATE="2020-06-01"

# ----------------------------------------------------------------------------------------------------------------

#================= Install system dependencies
RUN sudo apt-get update && sudo apt-get install -y \
    libssl-dev \
    libsasl2-dev \
    libxml2-dev \
    libicu-dev \
    bzip2 \
    liblzma-dev \
    libbz2-dev \
    subversion \
    curl \
    libmariadbclient-dev \
    libv8-dev \
    procps \
    systemd \
    libmagick++-dev \
    libssh2-1-dev \
    ssh \
    openssl \
    supervisor \
    passwd \
    vim

#================= Copy Rprofile.site - set repos and other options
COPY rplatform/Rprofile.site /usr/local/lib/R/etc/Rprofile.site

#=================  passwordless sudo
RUN echo 'ALL ALL = (ALL) NOPASSWD: ALL' >> /etc/sudoers

#================ Add rstudio user into sudo,staff,root groups
RUN usermod -a -G sudo,staff,root rstudio

#================= Add Roche certs
RUN sudo wget -O /usr/local/share/ca-certificates/Roche_G3_Root_CA.crt  http://certinfo.roche.com/rootcerts/Roche%20G3%20Root%20CA.crt 
RUN sudo update-ca-certificates 

#================= Remove openssl settings (two short keys)
#TODO: contact auth team regarding this issue
RUN sudo grep -v "^CipherString = DEFAULT@SECLEVEL=2" /etc/ssl/openssl.cnf > /tmp/openssl.fixed.cnf
RUN sudo mv /tmp/openssl.fixed.cnf /etc/ssl/openssl.cnf 

#================= Install dependencies
COPY rplatform/dependencies.yaml rplatform/.github_access_token.txt* /mnt/vol/
COPY rplatform/install_all_deps.R /mnt/vol/install_all_deps.R
RUN R -f /mnt/vol/install_all_deps.R

#================= Check & build package
COPY gDRcore/ /tmp/gDRcore/
RUN R -e 'remotes::install_deps(pkgdir = "/tmp/gDRcore/", dependencies = TRUE, repos = "https://cran.r-project.org")' && \
    R CMD INSTALL /tmp/gDRcore/ 

#================= Clean up
RUN sudo rm -rf /mnt/vol/* /tmp/gDRcore/
