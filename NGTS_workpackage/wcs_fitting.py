# -*- coding: utf-8 -*-

import tempfile
from catmatch import shift_wcs_axis
from catmatch import apply_correct, apply_correct_old
from catmatch import lmq_fit
from catmatch import calc_seps
from catmatch import load_wcs_from_keywords
from catmatch import correct_catfile
from multiprocessing import Pool as ThreadPool
from functools import partial
import casutools
from util import validate_headers
from astropy.io import fits
import os
from vector_plot import wcsf_QCheck
import numpy as np
from wcs_status import set_wcs_status
import fitsio
from astropy.io import fits as pf
from collections import namedtuple
try:
    import cPickle as pickle
except ImportError:
    import pickle
import json


class NullPool(object):

    def __init__(self, *args, **kwargs):
        pass

    def map(self, fn, args):
        return map(fn, args)


Catalogue = namedtuple('Catalogue', ['cat_name', 'ra_lims', 'dec_lims'])


def initialise_wcs_cache(fname, wcsref, thresh, verbose, force=False):
    print("Constructing initial wcs cache")
    catalogue_name = 'initial-catalogue.fits'
    casutools.imcore(fname, catalogue_name, threshold=thresh, verbose=verbose)
    casutools.wcsfit(fname, catalogue_name, catpath=wcsref, verbose=verbose)


def m_solve_images(filelist, outfile, dist_map, wcsref,
                   nproc=None,
                   thresh=20.0,
                   verbose=False,
                   catsrc='viz2mass'):

    infiles = []

    with open(filelist) as infile:
        for line in infile:
            image = line.strip()
            status_checks = ['ok', 'ok']

            if all(status == 'ok' for status in status_checks):
                infiles.append(image)

    initialise_wcs_cache(infiles[0], wcsref, thresh, verbose)

    fn = partial(handle_errors_in_casu_solve,
                 wcsref=wcsref,
                 dist_map=dist_map,
                 thresh=thresh,
                 verbose=verbose,
                 catsrc=catsrc)

    pool = ThreadPool(nproc)

    return pool.map(fn, infiles)


def handle_errors_in_casu_solve(casuin, *args, **kwargs):
    '''
  Catch any exceptions that may be thrown by the image solving routine to
  prevent the pipeline from crashing
  '''
    try:
        return_value = casu_solve(casuin, *args, **kwargs)
    except Exception as err:
        print "Exception handled in `casu_solve`: {}".format(str(err))
        set_wcs_status(casuin, succeeded=False)
    else:
        set_wcs_status(casuin, succeeded=True)
        return return_value


def casu_solve_old(casuin, wcsref,
                   dist_map={},
                   thresh=20,
                   verbose=False,
                   catsrc='viz2mass',
                   catpath=None):
    hdulist = fitsio.read_header(casuin)

    cen = [[dist_map['CRPIX1'], dist_map['CRPIX2']]]

    TEL_RA = hdulist['TEL_RA']
    TEL_DEC = hdulist['TEL_DEC']

    for key in dist_map:
        print key, dist_map[key], hdulist.get(key)

    apply_correct_old(dist_map, casuin, TEL_RA, TEL_DEC)

    catfile_name = casuin.replace('.fits', '.cat')
    casutools.imcore(casuin, catfile_name, threshold=thresh, verbose=verbose)

    cat_names = []
    RA_lims = []
    DEC_lims = []

    catpath = (catpath if catpath is not None
            else os.path.join(os.getcwd(), 'catcache'))
    for line in open(catpath + '/index'):
        vals = line.strip('\n').split(' ')
        cat_names += [vals[0]]
        RA_lims += [[float(vals[2]), float(vals[3])]]
        DEC_lims += [[float(vals[4]), float(vals[5])]]

    n = 0

    cat_name = cat_names[n]

    with pf.open(catpath + '/' + cat_name) as catd:
        catt = catd[1].data.copy()
    cat = {'ra': catt['ra'], 'dec': catt['dec'], 'Jmag': catt['Jmag']}

    apply_correct_old(dist_map, casuin, TEL_RA, TEL_DEC)

    with fitsio.FITS(catfile_name) as mycatt:
        mycat = {'Aper_flux_3': mycatt[1]['Aper_flux_3'][:]}
        my_X = mycatt[1]['x_coordinate'][:]
        my_Y = mycatt[1]['y_coordinate'][:]

    try:
        dist_map = shift_wcs_axis(dist_map, mycat, cat, RA_lims, DEC_lims, my_X, my_Y,
                                  TEL_RA, TEL_DEC,
                                  iters=10)
        dist_map = lmq_fit(dist_map, mycat, cat, RA_lims, DEC_lims, my_X, my_Y, TEL_RA,
                           TEL_DEC,
                           fitlist=['RA_s', 'DEC_s', 'CD1_1', 'CD2_2', 'CD1_2', 'CD2_1'])
    except IOError:
        print "Performing initial fit"
        casutools.wcsfit(casuin, catfile_name, catpath=wcsref, verbose=verbose)
        dist_map = shift_wcs_axis(casuin, catfile_name, thresh=thresh, iters=30)
        dist_map = lmq_fit(dist_map, mycat, cat, RA_lims, DEC_lims, my_X, my_Y, TEL_RA,
                           TEL_DEC,
                           fitlist=['RA_s', 'DEC_s', 'CD1_1', 'CD2_2', 'CD1_2', 'CD2_1'])

    apply_correct_old(dist_map, casuin, TEL_RA, TEL_DEC)

    # wcs keywords may have changed since imcore was done, so we have to update the RA and DEC values.
    correct_catfile(catfile_name, casuin, nstars=2000)

    # Now we're ready to solve wcs
    casutools.wcsfit(casuin, catfile_name, catpath=wcsref, verbose=verbose)

    # Do QC checks. should really break this out.

    plot = True
    wcsf_QCheck(mycat, casuin, os.path.basename(casuin).strip('.fits') + '.png', cat,
                RA_lims, DEC_lims, my_X, my_Y,
                plot=plot)

    return 'ok'


def casu_solve(casuin, wcsref, dist_map, thresh=3, verbose=False, catsrc='viz2mass'):

    hdulist = fits.getheader(casuin)

    apply_correct(dist_map, casuin)
    catpath = os.path.join(os.getcwd(), 'catcache')

    catfile_name = casuin.replace('.fits', '.cat')
    casutools.imcore(casuin, catfile_name, threshold=thresh, verbose=verbose,
            ipix=2, rcore=3)

    # Now we're ready to solve wcs
    casutools.wcsfit(casuin, catfile_name, catpath=wcsref, verbose=verbose)


    catalogue = compute_frame_limits(catpath)
    cat = reference_catalogue_objects(catalogue, catpath)

    with fits.open(catfile_name) as mycatt:
        mycatt_data = mycatt[1].data
        mycat = {'Aper_flux_3': mycatt_data['Aper_flux_3']}
        my_X = mycatt_data['x_coordinate']
        my_Y = mycatt_data['y_coordinate']

    # Do QC checks. should really break this out.
    wcsf_QCheck(mycat, casuin,
                os.path.basename(casuin).replace('.fits', '') + '.png', cat,
                catalogue.ra_lims, catalogue.dec_lims, my_X, my_Y,
                plot=True)

    return 'ok'


def extract_dist_map(filename):
    with open(filename) as infile:
        try:
            dist_map = json.load(infile)
        except ValueError as err:
            if 'No JSON object could be decoded' in str(err):
                infile.seek(0)
                dist_map = pickle.load(infile)
            else:
                raise

    if 'meta' in dist_map:
        print(json.dumps(dist_map['meta'], indent=2))
        dist_map = dist_map['wcs']

    if 'CD1_1' not in dist_map:
        raise KeyError("Cannot find valid wcs solution in map {}".format(
            json.dumps(dist_map)))

    return dist_map


def reference_catalogue_objects(catalogue, catpath):
    cat_name = os.path.join(catpath, catalogue.cat_name)

    with fits.open(cat_name) as catd:
        catt = catd[1].data.copy()

    return {'ra': catt['ra'], 'dec': catt['dec'], 'Jmag': catt['Jmag']}


def compute_frame_limits(catpath):
    index_filename = os.path.join(catpath, 'index')
    cat_names = []
    RA_lims = []
    DEC_lims = []

    for line in open(index_filename):
        vals = line.strip().split()
        cat_names += [vals[0]]
        RA_lims += [[float(vals[2]), float(vals[3])]]
        DEC_lims += [[float(vals[4]), float(vals[5])]]

    return Catalogue(cat_names[0], RA_lims, DEC_lims)
