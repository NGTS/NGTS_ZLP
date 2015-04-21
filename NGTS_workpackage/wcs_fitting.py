# -*- coding: utf-8 -*-

import tempfile
from catmatch import shift_wcs_axis
from catmatch import apply_correct
from catmatch import lmq_fit
from catmatch import calc_seps
from catmatch import load_wcs_from_keywords
from catmatch import correct_catfile
from multiprocessing import Pool as ThreadPool
from functools import partial
import casutools
from util import validate_headers
import fitsio
import astropy.io.fits as pf
import os
from vector_plot import wcsf_QCheck
import numpy as np
from wcs_status import set_wcs_status


class NullPool(object):

    def __init__(self, *args, **kwargs):
        pass

    def map(self, fn, args):
        return map(fn, args)


def initialise_wcs_cache(fname, catpath, wcsref, thresh, verbose, force=False):
    if force or not os.path.isdir(catpath):
        print("Constructing initial wcs cache")
        catalogue_name = 'initial-catalogue.fits'
        casutools.imcore(fname, catalogue_name, threshold=thresh, verbose=verbose)
        casutools.wcsfit(fname, catalogue_name, catpath=wcsref, verbose=verbose)


def m_solve_images(filelist, outfile, dist_map, wcsref,
                   nproc=None,
                   thresh=20.0,
                   verbose=False,
                   catsrc='viz2mass',
                   catpath=False):
    infiles = []

    with open(filelist) as infile:
        for line in infile:
            image = line.strip('\n')
            status_checks = ['ok', 'ok']

            if all(status == 'ok' for status in status_checks):
                infiles.append(image)

    initialise_wcs_cache(infiles[0], catpath, wcsref, thresh, verbose)

    fn = partial(handle_errors_in_casu_solve,
                 wcsref=wcsref,
                 dist_map=dist_map,
                 thresh=thresh,
                 verbose=verbose,
                 catsrc=catsrc,
                 catpath=catpath)

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


def casu_solve(casuin, wcsref,
               dist_map,
               thresh=20,
               verbose=False,
               catsrc='viz2mass',
               catpath=False):


    hdulist = fitsio.read_header(casuin)

    apply_correct(dist_map, casuin)

    with tempfile.NamedTemporaryFile(dir='.',
                                     suffix='.fits',
                                     prefix='catalogue.') as catfile:
        catfile_name = catfile.name

        casutools.imcore(casuin, catfile_name, threshold=thresh, verbose=verbose)
        catfile.seek(0)

        # Now we're ready to solve wcs
        casutools.wcsfit(casuin, catfile_name, catpath=wcsref, verbose=verbose)

        # Do QC checks. should really break this out.
        wcsf_QCheck(mycat, casuin, os.path.basename(casuin).strip('.fits') + '.png', cat,
                    RA_lims, DEC_lims, my_X, my_Y,
                    plot=True)

        return 'ok'
