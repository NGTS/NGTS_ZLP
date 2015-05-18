#!/usr/bin/env bash

set -e

if [[ ! "$(hostname -s)" =~ ngts.+ ]]; then
    echo "This script must be run on an ngts node, sorry!" >&2
    exit 1
fi

. testinit.sh

create_filelist(){
    local readonly filelist_name=/tmp/filelist
    find -L ${BASEDIR}/longtestdata -name '*.phot' | sed 's/\.phot$//' > ${filelist_name}
    echo "${filelist_name}"
}

main
