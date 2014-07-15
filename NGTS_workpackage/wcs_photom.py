# -*- coding: utf-8 -*-
from astropy.io import fits as pf
import os
import linecache
import threading
from os.path import isfile, join
import multiprocessing
from quality_checks import *
import casutools
from multiprocessing import Pool
from functools import partial
from super_sample import call_find_fwhm 

def m_wcs_photom(filelist,outlist,appsize,conf_file,cat_file,nproc=1,verbose=False):

  os.system('cp '+filelist+' '+outlist)

  infiles = []
  with open(filelist) as infile:
    for line in infile:
      image = line.strip('\n')

      status_checks = ['ok','ok']

      if all(status == 'ok' for status in status_checks):
        infiles.append(image)

  pool = Pool(nproc)

  fn = partial(wcs_photom, cat_file=cat_file, conf_file=conf_file, appsize=appsize, verbose=verbose)
  pool.map(fn, infiles)

  first_image = infiles[0] + '.phot'
  pf.setval(first_image,'SKY_MOVE',1,value=0)

  RA_shift, DEC_shift, tot_shift, RA, DEC = frame_shift(infiles[0],infiles[0])
  pf.setval(first_image,'RA_MOVE',1,value=RA_shift,comment='RA shift from previous image [arcseconds]')
  pf.setval(first_image,'DEC_MOVE',1,value=DEC_shift,comment='Dec shift from previous image [arcseconds]')
  pf.setval(first_image,'SKY_MOVE',1,value=tot_shift,comment='Total movement on sky [arcseconds]')

  pf.setval(first_image,'WCSF_RA',1,value=RA_shift,comment='RA center pix')
  pf.setval(first_image,'WCSF_DEC',1,value=DEC_shift,comment='Dec center pix')

  indexes = arange(1,len(infiles))
  fn = partial(m_frame_shift, infiles)
  pool.map(fn,indexes)

def wcs_photom(image,cat_file='nocat',conf_file='noconf',appsize=2.0,verbose=False):
  
  outname = image + '.phot'

  casutools.imcore_list(image, cat_file, outname, confidence_map=conf_file,rcore=appsize, noell=True,
            verbose=verbose)

   #    do some quality checks

  factor = 5
  size = 11
  stars = 100

  fwhm = call_find_fwhm(image,factor,size,stars,tag=image)
     
  cloud_status = cloud_check(image)

  pixel_fwhm = pf.getval(outname,'SEEING',1)
  plate_scale = pf.getval(outname,'CD1_1',1)
  seeing = round(plate_scale*pixel_fwhm*3600,2)

  pf.setval(outname,'CLOUD_S',1,value=round(cloud_status,2),comment='A measure of bulk structure in the image (S/N)')
  pf.setval(outname,'FWHM',1,value=pixel_fwhm,comment='[pixels] Average FWHM')
  pf.setval(outname,'SEEING',1,value=seeing,comment='[arcseconds] Average FWHM')

  return 'ok'
