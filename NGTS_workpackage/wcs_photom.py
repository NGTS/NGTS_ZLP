# -*- coding: utf-8 -*-
from astropy.io import fits as pf
import os
import linecache
import threading
from os.path import isfile, join
import multiprocessing
from quality_checks import *
import casutools
from multiprocessing.dummy import Pool as ThreadPool
from functools import partial


def m_wcs_photom(filelist,outlist,appsize,conf_file,cat_file,nproc=1,verbose=False):

  os.system('cp '+filelist+' '+outlist)

  infiles = []
  with open(filelist) as infile:
    for line in infile:
      image = line.strip('\n')

      status_checks = ['ok','ok']

      if all(status == 'ok' for status in status_checks):
        infiles.append(image)

  pool = ThreadPool(nproc)

  fn = partial(wcs_photom, cat_file=cat_file, conf_file=conf_file, appsize=appsize, verbose=verbose)
  pool.map(fn, infiles)

  first_image = infiles[0] + '.phot'
  pf.setval(first_image,'SHIFT',1,value=0)

  indexes = arange(1,len(infiles))
  fn = partial(m_frame_shift, infiles)
  pool.map(fn,indexes)

def wcs_photom(image,cat_file='nocat',conf_file='noconf',appsize=2.0,verbose=False):
  
  outname = image + '.phot'

  casutools.imcore_list(image, cat_file, outname, confidence_map=conf_file,rcore=appsize, noell=True,
            verbose=verbose)

   #    do some quality checks
     
  cloud_status = cloud_check(image)

  pixel_fwhm = pf.getval(outname,'SEEING',1)
  plate_scale = pf.getval(outname,'CD1_1',1)
  seeing = round(plate_scale*pixel_fwhm*3600,2)

  pf.setval(outname,'CLOUD_S',1,value=round(cloud_status,2),comment='A measure of bulk structure in the image (S/N)')
  pf.setval(outname,'FWHM',1,value=pixel_fwhm,comment='[pixels] Average FWHM')
  pf.setval(outname,'SEEING',1,value=seeing,comment='[arcseconds] Average FWHM')

  return 'ok'
