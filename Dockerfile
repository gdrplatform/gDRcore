ARG BASE_IMAGE=marcinkam/gdrshiny:1.0
FROM ${BASE_IMAGE}

# temporary fix
# GitHub token for downloading private dependencies
# Need to be defined after FROM as it flushes ARGs
ARG GITHUB_TOKEN

#================= Install dependencies
RUN mkdir -p /mnt/vol
COPY rplatform/dependencies.yaml rplatform/.github_access_token.txt* /mnt/vol
RUN echo "$GITHUB_TOKEN" >> /mnt/vol/.github_access_token.txt
RUN Rscript -e 'BiocManager::install("BiocStyle")'
#' @importFrom data.table := .SD .I
RUN Rscript -e "gDRstyle::installAllDeps()"

#================= Check & build package
COPY ./ /tmp/gDRcore/
RUN Rscript -e "gDRstyle::installLocalPackage('/tmp/gDRcore')"

#================= Clean up
RUN sudo rm -rf /mnt/vol/* /tmp/gDRcore/
