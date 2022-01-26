#/bin/sh

echo "Executing $0"
echo "Environment: ${rp_env}"
echo "Working directory: `pwd`"
echo "Working directory contains: `ls | tr '\n' ' '`"

# exit when any command fails
set -e

echo ">>>>>>>> Running linter"
Rscript -e "gDRstyle::lintPkgDirs('/mnt/vol/gDRcore')"

echo ">>>>> RUNNING UNIT TESTS"
Rscript -e "testthat::test_local(path = '/mnt/vol/gDRcore', stop_on_failure = TRUE)"

# TODO: fix the issue with R CMD CHECK
#echo ">>>>> RUNNING DEVTOOLS::CHECK()"
#sudo R CMD check --no-build-vignettes --no-manual --no-tests /mnt/vol/gDRcore

echo ">>>>>>>> RUNNING CHECK DEPENDENCIES"
Rscript -e "gDRstyle::checkDependencies(desc_path='/mnt/vol/gDRcore/DESCRIPTION', dep_path='/mnt/vol/rplatform/dependencies.yaml')"
