#!/usr/bin/env bash

if [[ "$#" -ne 4 ]]; then
    echo "Usage: $0 <casuin> <catname> <chainname> <catsrc>" >&2
    echo "See documentation: http://ngts.warwick.ac.uk/twiki/bin/view/Main/ZLPDocAstrometrySolver" >&2
    exit 1
fi

CASUIN=$1
MYCATNAME=$2
CHAIN_NAME=$3
CATSRC=$4

imcore ${CASUIN} noconf ${MYCATNAME} 2 2

python $(dirname $0)/emcee_catmatch.py ${CASUIN} ${MYCATNAME} ${CHAIN_NAME} ${CATSRC}
