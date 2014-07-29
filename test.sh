#!/usr/bin/env bash
set -e

BASEDIR=$(readlink -f $(dirname $0))

setup_environment() {
    export PYTHONPATH=${BASEDIR}:${BASEDIR}/testdata:${PYTHONPATH}
}

create_filelist(){
    local readonly filelist_name=/tmp/filelist
    find ${BASEDIR}/testdata -name '*.phot' | sed 's/\.phot$//' > ${filelist_name} 2>/dev/null
    echo "${filelist_name}"
}

run_test() {
    local readonly filelist_name=$(create_filelist)
    python ./bin/ZLP_create_outfile.py \
        --outdir ${BASEDIR}/testdata \
        ${filelist_name} \
        --apsize 2 \
        --nproc 1
}

assert_output() {
    python - <<EOF
import fitsio
import numpy as np
with fitsio.FITS("testdata/output.fits") as infile:
    flux = infile['flux'].read()

assert np.min(flux) == np.nan, "Error with min flux, should be nan, got {}".format(np.min(flux))
assert np.nanmin(flux) > 0, "Error with nanmin flux, should be > 0, got {}".format(np.nanmin(flux))
EOF
}

main() {
    (cd ${BASEDIR}
    setup_environment
    run_test
    assert_output
    )
}

main
