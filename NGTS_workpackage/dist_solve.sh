#!/usr/bin/env bash

CASUIN=$1
MYCATNAME=$2
CHAIN_NAME=$3
CATSRC=$4

imcore ${CASUIN} noconf ${MYCATNAME} 2 2

python emcee_catmatch.py ${CASUIN} ${MYCATNAME} ${CHAIN_NAME} ${CATSRC}