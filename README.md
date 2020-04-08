Functions to be used by the gDRshiny app ( https://github.roche.com/hafnerm6/gDRshiny )

SOP (Google Doc): https://docs.google.com/document/d/1om9rC3ynsCygfbeAFN1XpX4IPkEyPY8W7iDKmpYRw5w/edit#

Design of the mySQL database: https://drive.google.com/file/d/1h1U-lY_s_l1gnQHItIulCCvOcXnhKNMQ/view?usp=sharing

## How to launch gDR app with Docker container
### Prerequisites
* Docker
We recommend to use Docker Desktop available [here](https://www.docker.com/get-started)
* git
### Download gDR repositories
All the necessary files are available within gdrplatform organization on Roche GitHub
In order to download all the necessary files just use the following commands:
```
cd
mkdir -p git
cd git
git clone https://github.roche.com/gdrplatform/gDRcore.git
cd gDRcore
git clone https://github.roche.com/gdrplatform/gDRwrapper.git
git clone https://github.roche.com/gdrplatform/gDRshiny.git
git clone https://github.roche.com/gdrplatform/gdr_rp_test.git
git clone https://github.roche.com/gdrplatform/gDRsearch.git
git clone https://github.roche.com/gdrplatform/DataFrameMatrix.git
git clone https://github.roche.com/gdrplatform/gDRviz.git
```
**Note**: all the repositories will be cloned into the git directory on your home directory.

### Download, bind and run docker container
```
docker run -d -p 8788:8787 -v ~/git/gDRcore:/mnt/vol/gDR registry.rplatform.org:5000/githubroche/gdrplatform/gdr-shiny:0.0.0.9000
```

### Run dockerized RStudio
In your browser type: `http://localhost:8788/`. Then RStudio will be launched.

### Run Shipy app
In R console of dockerized RStudio type the following commands:
```
setwd("/mnt/vol/gDR/gDRshiny/gDRshiny/")
file.copy("inst/extdata/.Renviron", "~")
Sys.setenv("GDR_SHINY_ROOT_DIR" = "/home/rstudio/.gDR")
devtools::load_all("/mnt/vol/gDR/gDR")
devtools::load_all("/mnt/vol/gDR/gDRviz/gDRviz/")
devtools::load_all("/mnt/vol/gDR/gDRwrapper/gDRwrapper/")

devtools::load_all()
shiny::runApp("inst/shiny")
```
After that, Shiny app of gDR should start.
