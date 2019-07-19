#/bin/sh

echo "Executing $0"
echo "Environment: ${rp_env}"
echo "Working directory: `pwd`"
echo "Working directory contains: `ls | tr '\n' ' '`"

# exit when any command fails
set -e

echo ">>>>> RUNNING UNIT TESTS"
Rscript -e "devtools::test(pkg = '/mnt/vol/gDR', stop_on_failure = TRUE)"

echo ">>>>> RUNNING DEVTOOLS::CHECK()"
Rscript -e "devtools::check(pkg = '/mnt/vol/gDR', error_on = 'error', build_args = '--no-build-vignettes --no-manual')"
