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

def m_solve_images(filelist,outfile,dist_map,wcsref,nproc=None, thresh=20.0, verbose=False, catsrc='viz2mass', catpath=False):
  infiles = []

  with open(filelist) as infile:
    for line in infile:
      image = line.strip('\n')
      status_checks = ['ok','ok']

      if all(status == 'ok' for status in status_checks):
        infiles.append(image)

  if not os.path.isdir(catpath):
    print("Constructing initial wcs cache")
    catalogue_name = 'initial-catalogue.fits'
    casutools.imcore(infiles[0], catalogue_name, threshold=thresh, verbose=verbose)
    casutools.wcsfit(infiles[0], catalogue_name, catsrc='localfits', catpath=wcsref,
                     verbose=verbose)

  fn = partial(casu_solve,wcsref=wcsref,dist_map=dist_map,thresh=thresh, verbose=verbose, catsrc=catsrc, catpath=catpath)

  pool = ThreadPool(nproc)
  return pool.map(fn, infiles)

def casu_solve(casuin,wcsref,dist_map={},thresh=20, verbose=False,catsrc='viz2mass',catpath=False):

#  validate_headers(casuin)

#  dist_map = {'CD2_1': 1.5361223669861037e-06, 'CD2_2': 0.0014111069172090588, 'RA_s': -0.045791803176241935, 'CD1_2': -1.8903458858157504e-06, 'CD1_1': 0.0014114252281619941, 'CRVAL2': 49.678409065669463, 'CRPIX1': 1032.6704835673609, 'CRPIX2': 965.55322003866229, 'CRVAL1': 285.43804393221615, 'PV2_1': 1.0, 'PV2_3': 8.6870348412920961, 'PV2_5': 271.2258736568354, 'PV2_7': 13711.331560137938, 'DEC_s': 0.048582466566705385, 'CTYPE2': 'DEC--ZPN', 'CTYPE1': 'RA---ZPN'}

#  dist_map = {'CD2_1': -1.41054781057e-05, 'CD2_2': 0.00138875543313, 'RA_s': -0.045791803176241935, 'CD1_2': 1.41656423353e-05, 'CD1_1': 0.00138877871593, 'CRVAL2': 49.678409065669463, 'CRPIX1': 1024.0 - 30.358966, 'CRPIX2': 1024 - 1.76800, 'CRVAL1': 285.43804393221615, 'PV2_1': 0.99999, 'PV2_3': 8.11, 'PV2_5': 901.97, 'DEC_s': 0.048582466566705385, 'CTYPE2': 'DEC--ZPN', 'CTYPE1': 'RA---ZPN'}

  #best fit for june 5th

#  dist_map = {'CD2_1': 1.5361223669861037e-06, 'CD2_2': 0.0014111069172090588, 'RA_s': -0.45791803176241935, 'CD1_2': -1.8903458858157504e-06, 'CD1_1': 0.0014114252281619941, 'CRVAL2': 49.678568623942, 'CRPIX1': 1012.6704835673609, 'CRPIX2': 965.55322003866229, 'CRVAL1':285.438652515631, 'PV2_1': 1.0, 'PV2_3': 8.6870348412920961, 'PV2_5': 271.2258736568354, 'PV2_7': 13711.331560137938, 'DEC_s': 0.48582466566705385, 'CTYPE2': 'DEC--ZPN', 'CTYPE1': 'RA---ZPN'}


  #best fit for june 11th
#  dist_map = {'CD2_1': 1.5415410350820255e-06, 'CD2_2': 0.0014112047383328613, 'RA_s': -0.52982456599841854, 'CD1_2': -1.6950507735128331e-06, 'CD1_1': 0.0014112627531226202, 'CRVAL2': 49.59559537217558, 'CRPIX1': 1005.5871367339122, 'CRPIX2': 963.33376370419967, 'CRVAL1': 285.36598106573831, 'PV2_1': 1.0, 'PV2_2': 0.0, 'PV2_3': 9.3291677327369538, 'PV2_5': -542.20833721702797, 'PV2_7': 335143.48864371137, 'DEC_s': 0.40305908468262103, 'CTYPE2': 'DEC--ZPN', 'CTYPE1': 'RA---ZPN'}	


  # best fit for june 19th
#  dist_map = {'CD2_1': -9.3263106728144037e-07, 'CD2_2': 0.0013885458770754415, 'RA_s': -0.042574725586388078, 'CD1_2': 7.6821287599453613e-07, 'CD1_1': 0.0013884790887784609, 'CRVAL2': 49.099782962962607, 'CRPIX1': 1001.8384609589782, 'CRPIX2': 976.12370523564334, 'CRVAL1': 285.85331052744289, 'PV2_1': 1.0, 'PV2_3': 8.3820912409837796, 'PV2_5': 690.45326416643377, 'PV2_7': -12173.91504945435, 'DEC_s': -0.092710638655683769, 'CTYPE2': 'DEC--ZPN', 'CTYPE1': 'RA---ZPN'} 

  # guess the offset here to get rid of a lot of night to night uncertainty, if we assume initial solution to be 'nearly' correct. 

  # best fit for june 26th

#  dist_map = {'CD2_1': -1.4191656988099457e-08, 'CD2_2': 0.0013885222778780631, 'RA_s': 0.055870512902160017, 'CD1_2': -1.8513409749920829e-07, 'CD1_1': 0.0013883723486931344, 'CRVAL2': 49.414386141920922, 'CRPIX1': 1022.0594111585548, 'CRPIX2': 982.11907123259516, 'CRVAL1': 285.95176591229813, 'PV2_1': 1.0, 'PV2_2': 0.0, 'PV2_3': 7.898207770193908, 'PV2_5': 1392.9968242630193, 'PV2_7': -265152.6479961705, 'DEC_s': 0.22189846295346452, 'CTYPE2': 'DEC--ZPN', 'CTYPE1': 'RA---ZPN'}  

  hdulist = fitsio.read_header(casuin)

  cen = [[dist_map['CRPIX1'],dist_map['CRPIX2']]]

#  old_world = load_wcs_from_keywords(hdulist,cen)

  TEL_RA = hdulist['TEL_RA']
  TEL_DEC = hdulist['TEL_DEC']

#  dist_map['RA_s'] = (old_world[0][0] - TEL_RA)
#  dist_map['DEC_s'] = (old_world[0][1] - TEL_DEC) 

  for key in dist_map:
	print dist_map[key], hdulist.get(key)

  apply_correct(dist_map,casuin,TEL_RA,TEL_DEC) 

#  catpath = '/ngts/pipedev/AperturePhot/june_19th_test/catcache'

#  catpath = '/ngts/pipedev/AperturePhot/output/26th_june_output/catcache'

#  catpath = '/ngts/pipedev/AperturePhot/output/11th_june_output/catcache'

  with tempfile.NamedTemporaryFile(dir='.', suffix='.fits', prefix='catalogue.') as catfile:
    catfile_name = catfile.name

    casutools.imcore(casuin, catfile_name, threshold=thresh, verbose=verbose)
    catfile.seek(0)

    cat_names = []
    RA_lims = []
    DEC_lims = []

    for line in open(catpath+'/index'):
      vals = line.strip('\n').split(' ')
      cat_names += [vals[0]]
      RA_lims += [[float(vals[2]),float(vals[3])]]
      DEC_lims += [[float(vals[4]),float(vals[5])]]

    n = 0

    cat_name = cat_names[n]

    with pf.open(catpath+'/'+cat_name) as catd:
      catt = catd[1].data.copy()
    cat = {'ra':catt['ra'],'dec':catt['dec'],'Jmag':catt['Jmag']}

    apply_correct(dist_map,casuin,TEL_RA,TEL_DEC)

    with fitsio.FITS(catfile_name) as mycatt:
      mycat = {'Aper_flux_3':mycatt[1]['Aper_flux_3'][:]}  
      my_X = mycatt[1]['x_coordinate'][:]
      my_Y = mycatt[1]['y_coordinate'][:]

  # We need the extra correction here before we do the wcsfit, because the TEL RA and DEC measurements are not
  # always precise enough for the fit to work. This simply shifts CRVALS to align with CRPIX
    try:
      dist_map = shift_wcs_axis(dist_map,mycat,cat,RA_lims,DEC_lims,my_X,my_Y,TEL_RA,TEL_DEC,iters=3)
      dist_map = lmq_fit(dist_map,mycat,cat,RA_lims,DEC_lims,my_X,my_Y,TEL_RA,TEL_DEC,fitlist=['RA_s','DEC_s'])
#      dist_map = lmq_fit(dist_map,mycat,cat,RA_lims,DEC_lims,my_X,my_Y,TEL_RA,TEL_DEC)
    except IOError:
      print "Performing initial fit"
      casutools.wcsfit(casuin, catfile_name, catpath=wcsref, catsrc=catsrc, verbose=verbose)
      dist_map = shift_wcs_axis(casuin, catfile_name, thresh=thresh, iters=30)
      # make mag limited version should go in here

    apply_correct(dist_map,casuin,TEL_RA,TEL_DEC)

# wcs keywords may have changed since imcore was done, so we have to update the RA and DEC values.
    correct_catfile(catfile_name,casuin,nstars=2000)

# Now we're ready to solve wcs
    casutools.wcsfit(casuin, catfile_name, catpath=wcsref, catsrc=catsrc, verbose=verbose)

# Do QC checks. should really break this out.
    plot = True
    wcsf_QCheck(mycat,casuin,casuin.strip('.fits')+'.png',cat,RA_lims,DEC_lims,my_X,my_Y,plot=plot)

    return 'ok'

# vim: sw=2
