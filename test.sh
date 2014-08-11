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
    for hdu in ['flux', 'fluxerr', 'hjd']:
        data = infile[hdu].read()

        assert np.isnan(np.min(data)), "Error with min data ({}), should be nan, got {}".format(hdu, np.min(data))
        assert np.nanmin(data) > 0, "Error with nanmin data ({}), should be > 0, got {}".format(hdu, np.nanmin(data))
EOF
}

main() {
    (cd ${BASEDIR}
    setup_environment

    echo -n "Running test... "
    set +e
    run_test 2>/dev/null >/dev/null
    set -e
    assert_output
    if [[ "$?" == "0" ]]; then
        echo "Pass"
    else
        echo "Fail"
        exit 1
    fi
    )
}

main
