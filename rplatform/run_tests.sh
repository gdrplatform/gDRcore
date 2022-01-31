#/bin/sh
repo_path="$1"

echo "Executing $0"
echo "Environment: ${rp_env}"
echo "Working directory: `pwd`"
echo "Working directory contains: `ls | tr '\n' ' '`"

# exit when any command fails
set -e

echo ">>>>>>>> Running linter"
Rscript -e "gDRstyle::lintPkgInDir('$repo_path')"

echo ">>>>> RUNNING UNIT TESTS"
Rscript -e "testthat::test_local(path = '$repo_path', stop_on_failure = TRUE)"

# TODO: fix the issue with R CMD CHECK
#echo ">>>>> RUNNING DEVTOOLS::CHECK()"
#sudo R CMD check --no-build-vignettes --no-manual --no-tests /mnt/vol/gDRwrapper
