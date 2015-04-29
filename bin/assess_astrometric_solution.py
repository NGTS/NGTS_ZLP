#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from NGTS_workpackage.wcs_fitting import (initialise_wcs_cache, casu_solve, )
import fitsio
import shutil
import os
try:
    import cPickle as pickle
except ImportError:
    import pickle
import json


def image_size(fname):
    '''
    Return the image size from the Primary HDU
    '''
    return fitsio.read(fname).shape


def remove_overscan_strips(fname):
    '''
    Given an image of dimensions 2048x2088 remove the two overscan strips so
    the solution will work on solved frames
    '''
    if not image_size(fname) == (2048, 2088):
        raise RuntimeError("Frame has dimensions {}, this is likely not a raw "
                           "NGTS frame".format(image_size(fname)))

    with fitsio.FITS(fname) as infile:
        image_data = infile[0].read()
        header = infile[0].read_header()

    without_overscan = image_data[:, 20:-20]
    out_fname = os.path.splitext(fname)[0] + '_without_overscan_strips.fits'
    with fitsio.FITS(out_fname, 'rw', clobber=True) as outfile:
        outfile.write(without_overscan, header=header)

    return out_fname


def extract_dist_map(filename):
    with open(filename) as infile:
        try:
            dist_map = json.load(infile)
        except ValueError as err:
            if 'No JSON object could be decoded' in str(err):
                dist_map = pickle.load(infile)
            else:
                raise
        finally:
            return dist_map


def main(args):
    dist_map = extract_dist_map(args.solution)
    if 'meta' in dist_map:
        print(json.dumps(dist_map['meta'], indent=2))
        dist_map = dist_map['wcs']

    if 'CD1_1' not in dist_map:
        raise KeyError("Cannot find valid wcs solution in map {}".format(
            json.dumps(dist_map)))

    catpath = os.path.join(os.getcwd(), 'catcache')
    if os.path.isdir(catpath):
        shutil.rmtree(catpath)

    if not args.reduced:
        fname = remove_overscan_strips(args.filename)
    else:
        fname = args.filename

    print('Frame {} will be altered'.format(fname))

    assert image_size(fname) == (2048, 2048), (
        "Image incorrect shape, "
        "should be 2048x2048, is {}".format(image_size(fname)))

    initialise_wcs_cache(fname, wcsref=args.wcsref,
                         thresh=20.0,
                         verbose=False,
                         force=True)
    casu_solve(fname, wcsref=args.wcsref, dist_map=dist_map)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('-s', '--solution', required=True)
    parser.add_argument('-w', '--wcsref', required=True)
    parser.add_argument('--reduced',
                        required=False,
                        action='store_true',
                        default=False,
                        help='Is the frame reduced?')
    main(parser.parse_args())
