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
from astropy.io import fits
import os
from vector_plot import wcsf_QCheck
import numpy as np
from wcs_status import set_wcs_status
from collections import namedtuple


class NullPool(object):

    def __init__(self, *args, **kwargs):
        pass

    def map(self, fn, args):
        return map(fn, args)


Catalogue = namedtuple('Catalogue', ['cat_name', 'ra_lims', 'dec_lims'])


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
            image = line.strip()
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


def casu_solve(casuin, wcsref, dist_map, thresh=20, verbose=False, catsrc='viz2mass'):

    hdulist = fits.getheader(casuin)

    apply_correct(dist_map, casuin)
    catpath = os.path.join(os.getcwd(), 'catcache')

    with tempfile.NamedTemporaryFile(dir='.',
                                     suffix='.fits',
                                     prefix='catalogue.') as catfile:
        catfile_name = catfile.name

        casutools.imcore(casuin, catfile_name, threshold=thresh, verbose=verbose)
        catfile.seek(0)

        # Now we're ready to solve wcs
        casutools.wcsfit(casuin, catfile_name, catpath=wcsref, verbose=verbose)

        catfile.seek(0)

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
