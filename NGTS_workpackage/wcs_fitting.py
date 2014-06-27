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

def m_solve_images(filelist, outfile, nproc=None, thresh=20.0, verbose=False, catsrc='viz2mass', catpath=False):
  infiles = []
  with open(filelist) as infile:
    for line in infile:
      image = line.strip('\n')
      status_checks = ['ok','ok']

      if all(status == 'ok' for status in status_checks):
        infiles.append(image)

  fn = partial(casu_solve, thresh=thresh, verbose=verbose, catsrc=catsrc, catpath=catpath)

  pool = ThreadPool(nproc)
  return pool.map(fn, infiles)

def casu_solve(casuin, thresh=20, verbose=False,catsrc='viz2mass',catpath=False):

#  validate_headers(casuin)

  best_fit = {'CD2_1': 1.5361223669861037e-06, 'CD2_2': 0.0014111069172090588, 'RA_s': -0.045791803176241935, 'CD1_2': -1.8903458858157504e-06, 'CD1_1': 0.0014114252281619941, 'CRVAL2': 49.678409065669463, 'CRPIX1': 1032.6704835673609, 'CRPIX2': 965.55322003866229, 'CRVAL1': 285.43804393221615, 'PV2_1': 1.0, 'PV2_3': 8.6870348412920961, 'PV2_5': 271.2258736568354, 'PV2_7': 13711.331560137938, 'DEC_s': 0.048582466566705385, 'CTYPE2': 'DEC--ZPN', 'CTYPE1': 'RA---ZPN'}

  best_fit = {'CD2_1': -1.41054781057e-05, 'CD2_2': 0.00138875543313, 'RA_s': -0.045791803176241935, 'CD1_2': 1.41656423353e-05, 'CD1_1': 0.00138877871593, 'CRVAL2': 49.678409065669463, 'CRPIX1': 1024.0 - 30.358966, 'CRPIX2': 1024 - 1.76800, 'CRVAL1': 285.43804393221615, 'PV2_1': 0.99999, 'PV2_3': 8.11, 'PV2_5': 901.97, 'DEC_s': 0.048582466566705385, 'CTYPE2': 'DEC--ZPN', 'CTYPE1': 'RA---ZPN'}

  #best fit for june 11th
  best_fit = {'CD2_1': 1.5380488570791813e-06, 'CD2_2': 0.0014111274390500408, 'RA_s': -0.52908717651493298, 'CD1_2': -1.6695976562054481e-06, 'CD1_1': 0.0014111916879534069, 'CRVAL2': 49.596298794653762, 'CRPIX1': 1005.9291412514006, 'CRPIX2': 963.82708264715689, 'CRVAL1': 285.36671845522181, 'PV2_1': 1.0, 'PV2_3': 8.9692393356428699, 'PV2_5': 155.22228777185455, 'PV2_7': -33859.107455138263, 'DEC_s': 0.40376250716080447, 'CTYPE2': 'DEC--ZPN', 'CTYPE1': 'RA---ZPN'}	


  # best fit (so far) for june 19th

  best_fit = {'CD2_1': -1.763272435228194e-06, 'CD2_2': 0.0013887753673947855, 'RA_s': 0.0051376769131176378, 'CD1_2': 1.6937843402091391e-06, 'CD1_1': 0.0013885358074354759, 'CRVAL2': 49.096226173873745, 'CRPIX1': 1024.4073925114958, 'CRPIX2': 973.5786968881423, 'CRVAL1': 285.90102292994243, 'PV2_1': 1.0, 'PV2_3': 9.0780825035408235, 'PV2_5': 348.76108086340531, 'PV2_7': 21640.644597109309, 'DEC_s': -0.096267427744546094, 'CTYPE2': 'DEC--ZPN', 'CTYPE1': 'RA---ZPN'}

  # guess the offset here to get rid of a lot of night to night uncertainty, if we assume initial solution to be 'nearly' correct. 

  hdulist = fitsio.read_header(casuin)

#  cen = [[best_fit['CRPIX1'],best_fit['CRPIX2']]]

#  old_world = load_wcs_from_keywords(hdulist,cen)

  TEL_RA = hdulist['TEL_RA']
  TEL_DEC = hdulist['TEL_DEC']

#  best_fit['RA_s'] = (old_world[0][0] - TEL_RA)
#  best_fit['DEC_s'] = (old_world[0][1] - TEL_DEC) 

  apply_correct(best_fit,casuin,TEL_RA,TEL_DEC) 

  catpath = '/ngts/pipedev/AperturePhot/june_19th_test/'

  with tempfile.NamedTemporaryFile(dir='.', suffix='.fits', prefix='catalogue.') as catfile:
    catfile_name = catfile.name

    casutools.imcore(casuin, catfile_name, threshold=thresh, verbose=verbose)
    catfile.seek(0)

    cat_names = []
    RA_lims = []
    DEC_lims = []

    for line in open(catpath+'catcache/index'):
      vals = line.strip('\n').split(' ')
      cat_names += [vals[0]]
      RA_lims += [[float(vals[2]),float(vals[3])]]
      DEC_lims += [[float(vals[4]),float(vals[5])]]

    n = 0

    cat_name = cat_names[n]

    with pf.open(catpath+'catcache/'+cat_name) as catd:
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
#      best_fit = shift_wcs_axis(best_fit,mycat,cat,RA_lims,DEC_lims,my_X,my_Y,TEL_RA,TEL_DEC,iters=1)
      best_fit = lmq_fit(best_fit,mycat,cat,RA_lims,DEC_lims,my_X,my_Y,TEL_RA,TEL_DEC,fitlist=['RA_s','DEC_s'])
#      best_fit = lmq_fit(best_fit,mycat,cat,RA_lims,DEC_lims,my_X,my_Y,TEL_RA,TEL_DEC)
    except IOError:
      print "Performing initial fit"
      casutools.wcsfit(casuin, catfile_name, verbose=verbose)
      best_fit = shift_wcs_axis(casuin, catfile_name, thresh=thresh, iters=30)
      # make mag limited version should go in here

    apply_correct(best_fit,casuin,TEL_RA,TEL_DEC)

# wcs keywords may have changed since imcore was done, so we have to update the RA and DEC values.
    correct_catfile(catfile_name,casuin,nstars=2000)

# Now we're ready to solve wcs
    casutools.wcsfit(casuin, catfile_name, verbose=verbose)

# Do QC checks. plotting disabled for now.
    plot = True
    wcsf_QCheck(mycat,casuin,casuin.strip('.fits')+'.png',cat,RA_lims,DEC_lims,my_X,my_Y,plot=plot)

    return 'ok'

