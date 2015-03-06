#!/usr/bin/env bash
set -e

# Arguments:
#   casuin
#   mycatname
#   chain_name
#   catsrc

main() {
    emcee_catmatch.py \
        emcee-testdata/casuin.fits \
        emcee-testdata/mycatname.cat \
        emcee-testdata/chain_name.fits \
        emcee-testdata/catcache \
        --nwalkers 30 \
        --nruns 2 \
        --nthreads 1
}

main
