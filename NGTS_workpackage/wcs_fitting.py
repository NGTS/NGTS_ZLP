# -*- coding: utf-8 -*-

import tempfile
from catmatch import shift_wcs_axis
from catmatch import apply_correct
from catmatch import lmq_fit
from catmatch import calc_seps
from catmatch import load_wcs_from_keywords
from catmatch import correct_catfile
from multiprocessing.dummy import Pool as ThreadPool
from functools import partial
import casutools
from util import validate_headers
import fitsio
import astropy.io.fits as pf
import os
from vector_plot import plot_differences

def m_solve_images(filelist, outfile, nproc=None, thresh=20.0, verbose=False, reset=False):

  infiles = []
  with open(filelist) as infile:
    for line in infile:
      image = line.strip('\n')
      status_checks = ['ok','ok']

      if all(status == 'ok' for status in status_checks):
        infiles.append(image)

  fn = partial(casu_solve, thresh=thresh, verbose=verbose)

  pool = ThreadPool(nproc)
  return pool.map(fn, infiles)

def casu_solve(casuin, thresh=20, verbose=False):

  validate_headers(casuin)

  best_fit = {'CD2_1': 1.5361223669861037e-06, 'CD2_2': 0.0014111069172090588, 'RA_s': -0.45791803176241935, 'CD1_2': -1.8903458858157504e-06, 'CD1_1': 0.0014114252281619941, 'CRVAL2': 49.678409065669463, 'CRPIX1': 1032.6704835673609, 'CRPIX2': 965.55322003866229, 'CRVAL1': 285.43804393221615, 'PV2_1': 1.0, 'PV2_3': 8.6870348412920961, 'PV2_5': 271.2258736568354, 'PV2_7': 13711.331560137938, 'DEC_s': 0.48582466566705385, 'CTYPE2': 'DEC--ZPN', 'CTYPE1': 'RA---ZPN'}


  hdulist = fitsio.read_header(casuin)
  TEL_RA = hdulist['TEL_RA']
  TEL_DEC = hdulist['TEL_DEC']

  with tempfile.NamedTemporaryFile(dir='.', suffix='.fits', prefix='catalogue.') as catfile:
    catfile_name = catfile.name

    casutools.imcore(casuin, catfile_name, threshold=thresh, verbose=verbose)
    catfile.seek(0)

    cat_names = []
    RA_lims = []
    DEC_lims = []
    for line in open('catcache/index'):
      vals = line.strip('\n').split(' ')
      cat_names += [vals[0]]
      RA_lims += [[float(vals[2]),float(vals[3])]]
      DEC_lims += [[float(vals[4]),float(vals[5])]]

    n = 0

    cat_name = cat_names[n]

    with pf.open('catcache/'+cat_name) as catd:
      catt = catd[1].data.copy()
    cat = {'ra':catt['ra'],'dec':catt['dec'],'Jmag':catt['Jmag']}

    apply_correct(best_fit,casuin,TEL_RA,TEL_DEC)

    with fitsio.FITS(catfile_name) as mycatt:
      mycat = {'Aper_flux_3':mycatt[1]['Aper_flux_3'][:]}  
      my_X = mycatt[1]['x_coordinate'][:]
      my_Y = mycatt[1]['y_coordinate'][:]

  # We need the extra correction here before we do the wcsfit, because the TEL RA and DEC measurements are not
  # always precise enough for the fit to work. This simply shifts CRVALS to align with CRPIX
    try:
      best_fit = shift_wcs_axis(best_fit,mycat,cat,RA_lims,DEC_lims,my_X,my_Y,TEL_RA,TEL_DEC,iters=1)
#      best_fit = lmq_fit(best_fit,mycat,cat,RA_lims,DEC_lims,my_X,my_Y,TEL_RA,TEL_DEC,fitlist=['RA_s','DEC_s'])
#      best_fit = lmq_fit(best_fit,mycat,cat,RA_lims,DEC_lims,my_X,my_Y,TEL_RA,TEL_DEC)
    except IOError:
      print "Performing initial fit"
      casutools.wcsfit(casuin, catfile_name, verbose=verbose)      
      best_fit = shift_wcs_axis(best_fit,mycat,cat,RA_lims,DEC_lims,my_X,my_Y,TEL_RA,TEL_DEC,iters=1)

    apply_correct(best_fit,casuin,TEL_RA,TEL_DEC)

# wcs keywords may have changed since imcore was done, so we have to update the RA and DEC values.
    correct_catfile(catfile_name,casuin,nstars=2000)

# Now we're ready to solve wcs
    casutools.wcsfit(casuin, catfile_name, verbose=verbose)

# make the QC vector plot

#    plot_differences(mycat,casuin,casuin.strip('.fits')+'.png',cat,RA_lims,DEC_lims,my_X,my_Y)

    return 'ok'