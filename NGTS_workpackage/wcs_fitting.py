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

    cen = [[dist_map['CRPIX1'], dist_map['CRPIX2']]]

    TEL_RA = hdulist['TEL_RA']
    TEL_DEC = hdulist['TEL_DEC']


    for key in dist_map:
        print key, dist_map[key], hdulist.get(key)

    apply_correct(dist_map, casuin, TEL_RA, TEL_DEC)


    with tempfile.NamedTemporaryFile(dir='.',
                                     suffix='.fits',
                                     prefix='catalogue.') as catfile:
        catfile_name = catfile.name

        casutools.imcore(casuin, catfile_name, threshold=thresh, verbose=verbose)
        catfile.seek(0)

        cat_names = []
        RA_lims = []
        DEC_lims = []

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

        apply_correct(dist_map, casuin, TEL_RA, TEL_DEC)

        with fitsio.FITS(catfile_name) as mycatt:
            mycat = {'Aper_flux_3': mycatt[1]['Aper_flux_3'][:]}
            my_X = mycatt[1]['x_coordinate'][:]
            my_Y = mycatt[1]['y_coordinate'][:]

        # We need the extra correction here before we do the wcsfit, because the TEL RA and DEC measurements are not
        # always precise enough for the fit to work. This simply shifts CRVALS to align with CRPIX

        try:
            dist_map = shift_wcs_axis(dist_map, mycat, cat, RA_lims, DEC_lims, my_X,
                                      my_Y, TEL_RA, TEL_DEC,
                                      iters=10)
            dist_map = lmq_fit(
                dist_map, mycat, cat, RA_lims, DEC_lims, my_X, my_Y, TEL_RA, TEL_DEC,
                fitlist=['RA_s', 'DEC_s', 'CD1_1', 'CD2_2', 'CD1_2', 'CD2_1'])
        except IOError:
            print "Performing initial fit"
            casutools.wcsfit(casuin, catfile_name, catpath=wcsref, verbose=verbose)
            dist_map = shift_wcs_axis(casuin, catfile_name, thresh=thresh, iters=30)
            dist_map = lmq_fit(
                dist_map, mycat, cat, RA_lims, DEC_lims, my_X, my_Y, TEL_RA, TEL_DEC,
                fitlist=['RA_s', 'DEC_s', 'CD1_1', 'CD2_2', 'CD1_2', 'CD2_1'])

        apply_correct(dist_map, casuin, TEL_RA, TEL_DEC)

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
