# -*- coding: utf-8 -*-
import pyfits as pf
import os
import linecache
import threading
from os.path import isfile, join
from util import thread_alloc, status_update

def m_wcs_photom(filelist,outlist,appsize,conf_file,cat_file,nproc=1,verbose=False):


  nfiles = 0
  for line in open(filelist):
    nfiles += 1

  os.system('cp '+filelist+' '+outlist)

  starts, ends = thread_alloc(nfiles,nproc)

  threads = []
  for i in range(0,nproc):
    t = threading.Thread(target=wcs_photom, args = (filelist,outlist,starts[i],ends[i],i+1,conf_file,cat_file,appsize,verbose))
    threads.append(t)
  [x.start() for x in threads]
  [x.join() for x in threads]

def wcs_photom(filelist,outlist,minlen,maxlen,thread,conf_file,cat_file,appsize,verbose=False):

  for i in range(minlen,maxlen):
    percent = 100*((i+1)-minlen)/(maxlen-minlen)
    line = linecache.getline(filelist,i).rstrip('\n')
    image = line.split(' ')[0]
    status_check = line.split(' ')[1:]
    if(all([status == 'ok' for status in status_check])):
      phot_status = casu_photom(image,conf_file,cat_file,appsize,verbose)
      status_update(outlist,line,line+' '+phot_status)
    else:
      status_update(outlist,line,line+' not_done')
    if verbose == True:
      print 'process ',thread,' ',percent,'% complete'

  linecache.clearcache()

def casu_photom(image,conf_file,cat_file,appsize,verbose=False):
    outname = image+'.phot'
    screen = 'imcore_list '+image+' '+conf_file+' '+cat_file+' '+outname+' --rcore='+str(appsize)+' --noell'
    if verbose == True:
      print screen
    os.system(screen)

    status = 'ok'

    return status

