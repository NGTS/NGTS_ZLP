#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from NGTS_workpackage.wcs_fitting import (initialise_wcs_cache,
                                          casu_solve,
                                          )
try:
    import cPickle as pickle
except ImportError:
    import pickle


def main(args):
    with open(args.solution) as infile:
        dist_map = pickle.load(infile)

    catpath = 'catcache'
    initialise_wcs_cache(args.filename, catpath, args.wcsref,
                         thresh=20.0, verbose=False)
    casu_solve(args.filename, args.wcsref, dist_map=dist_map, catpath=catpath)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('-s', '--solution', required=True)
    parser.add_argument('-w', '--wcsref', required=True)
    main(parser.parse_args())
