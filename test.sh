#!/usr/bin/env bash
set -e

abspath() {
    python -c "import os; print os.path.realpath('${1}')"
}

BASEDIR=$(abspath $(dirname $0))
PIPELINEDIR="${BASEDIR}/../zlp-script"
STDERRFILE=/tmp/test.stderr
STDOUTFILE=/tmp/test.stdout


setup_environment() {
    export PYTHONPATH=${BASEDIR}:${BASEDIR}/testdata:${PIPELINEDIR}/scripts:${PYTHONPATH}
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
    assert_npts_correct
}

assert_npts_correct() {
python - <<EOF
import fitsio
with fitsio.FITS("testdata/output.fits") as infile:
    catalogue = infile['catalogue']
    keys = catalogue.get_colnames()
    value_ind = keys.index('NPTS')
    nrows = catalogue.get_nrows()
    flux = infile['flux'].read()
    assert flux.shape[0] == nrows

    for (lc, cat_row) in zip(flux, catalogue):
        value = cat_row[value_ind]
        target = lc[lc == lc].size
        assert value == target, (value, target)
EOF
}

# XXX Deprecated, do not use
assert_nans_present() {
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

check_for_failure() {
    if [[ "$?" == "0" ]]; then
        echo "Pass"
    else
        echo "Fail"
        cat ${STDERRFILE}
        exit 1
    fi
}

main() {
    (cd ${BASEDIR}
    setup_environment

    echo -n "Running test... "
    set +e
    run_test 2>${STDERRFILE} >${STDOUTFILE}
    check_for_failure
    set -e
    echo -n "Running test... "
    assert_output
    check_for_failure
    )
}

main
