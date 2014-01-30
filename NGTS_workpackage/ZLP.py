# -*- coding: utf-8 -*-
import os
import linecache
import threading
from os.path import isfile, join
from util import thread_alloc
import numpy
import sys
import astropy.io.fits as pf
from NGTS_workpackage import *

def m_call_ZLP(filelist,outlist,conf_file,cat_file,appsize,nproc=1,thresh=20.0,verbose=False):

  nfiles = 0
  for line in open(filelist):
    nfiles += 1

  starts, ends = thread_alloc(nfiles,nproc)

  threads = []
  for i in range(0,nproc):
    t = threading.Thread(target=call_ZLP, args = (filelist,outlist,conf_file,cat_file,appsize,starts[i],ends[i],i+1,thresh,verbose))
    threads.append(t)
  [x.start() for x in threads]
  [x.join() for x in threads]

def call_ZLP(filelist,outlist,conf_file,cat_file,appsize,minlen,maxlen,thread,thresh=20.0,verbose=False):

  for i in range(minlen,maxlen):
    line = linecache.getline(filelist,i).rstrip('\n')
    image = line.split(' ')[0]
    cal_status = line.split(' ')[1]

    wcs_status = 'not_done'
    phot_status = 'not_done'

    if cal_status == 'ok':
      wcs_status = casu_solve(image,thresh,thread=thread,verbose=verbose)
      if wcs_status == 'ok':
        phot_status = casu_photom(image,conf_file,cat_file,appsize,verbose)

    with open(outlist,'a') as outfile:
      outfile.write(image + '.phot ' + cal_status + ' ' + wcs_status + ' ' + phot_status + '\n')

    if verbose==True:
      percent = 100.0*((i+1)-minlen)/(maxlen-minlen)
      print 'Process ',thread,' is ',percent,'% complete'

  linecache.clearcache()
