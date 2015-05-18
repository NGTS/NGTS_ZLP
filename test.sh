#!/usr/bin/env bash
set -e

. testinit.sh


create_filelist(){
    local readonly filelist_name=/tmp/filelist
    # Randomise the file names
    find ${BASEDIR}/testdata -name '*.phot' | sed 's/\.phot$//' | sort -R > ${filelist_name} 2>/dev/null
    echo "${filelist_name}"
}

main
