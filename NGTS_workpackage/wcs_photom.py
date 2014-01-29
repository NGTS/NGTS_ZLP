# -*- coding: utf-8 -*-
from astropy.io import fits as pf
import os
import linecache
import threading
from os.path import isfile, join
from util import thread_alloc, status_update
import multiprocessing

def m_wcs_photom(filelist,outlist,appsize,conf_file,cat_file,nproc=1,verbose=False):
  from quality_checks import *

  nfiles = 0
  for line in open(filelist):
    nfiles += 1

  os.system('cp '+filelist+' '+outlist)

  starts, ends = thread_alloc(nfiles,nproc)

  process = []
  for i in range(0,nproc):
    p = multiprocessing.Process(target=wcs_photom, args = (filelist,outlist,starts[i],ends[i],i+1,conf_file,cat_file,appsize,verbose))
    process.append(p)
  [x.start() for x in process]
  [x.join() for x in process]


  files_done = 0
  for name in open(outlist):
    files_done += 1

  first_image = linecache.getline(outlist,1).split(' ')[0] + '.phot'
  pf.setval(first_image,'SHIFT',1,value='0')

  starts, ends = thread_alloc(files_done-1,nproc)
  starts = starts + 1
  ends = ends + 1

  process = []
  for i in range(0,nproc):
    p = multiprocessing.Process(target=m_frame_shift, args=(outlist,starts[i],ends[i]))
    process.append(p)
  [x.start() for x in process]
  [x.join() for x in process]

def wcs_photom(filelist,outlist,minlen,maxlen,thread,conf_file,cat_file,appsize,verbose=False):
  from quality_checks import *
  
  first_frame = True

  for i in range(minlen,maxlen):
    percent = 100*((i+1)-minlen)/(maxlen-minlen)
    line = linecache.getline(filelist,i).rstrip('\n')
    image = line.split(' ')[0]
    outname = image + '.phot'
    status_check = line.split(' ')[1:]
    if(all([status == 'ok' for status in status_check])):
      phot_status = casu_photom(image,conf_file,cat_file,appsize,verbose)
      status_update(outlist,line,line+' '+phot_status)

     #    do some quality checks
      fwhm = get_fwhm(outname,int(appsize))
      cloud_status = cloud_check(image)

      pf.setval(outname,'CLOUD_S',1,value=cloud_status)
      pf.setval(outname,'FWHM',1,value=fwhm)
    else:
      status_update(outlist,line,line+' not_done')
    if verbose == True:
      print 'process ',thread,' ',percent,'% complete'

  linecache.clearcache()


def casu_photom(image,conf_file,cat_file,appsize,verbose=False):
    outname = image + '.phot'
    screen = 'imcore_list '+image+' '+conf_file+' '+cat_file+' '+outname+' --rcore='+str(appsize)+' --noell'
    if verbose == True:
      print screen
    os.system(screen)

    status = 'ok'

    return status

