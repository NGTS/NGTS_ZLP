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

main() {
    (cd ${BASEDIR}
    setup_environment
    run_test
    )
}

main
