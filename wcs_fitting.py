# -*- coding: utf-8 -*-
import os
import linecache
import threading
from os.path import isfile, join
from util import thread_alloc, status_update
import numpy
import sys
import astropy.io.fits as pf
from catmatch import *

def m_solve_images(filelist,outfile,nproc=1,thresh=20.0,verbose=False):

  nfiles = 0
  for line in open(filelist):
    nfiles += 1

  os.system('cp '+filelist+' '+outfile)

  starts, ends = thread_alloc(nfiles,nproc)

  threads = []
  for i in range(0,nproc):
    t = threading.Thread(target=solve_images, args = (filelist,outfile,starts[i],ends[i],i+1,thresh,verbose))
    threads.append(t)
  [x.start() for x in threads]
  [x.join() for x in threads]

def solve_images(filelist,outfile,minlen,maxlen,thread,thresh=20.0,verbose=False):

  for i in range(minlen,maxlen):
    line = linecache.getline(filelist,i).rstrip('\n')
    image = line.split(' ')[0]
    status_check = line.split(' ')[1:]
    if(all([status=='ok' for status in status_check])):
      fit_status = casu_solve(image,thresh,thread=thread,verbose=verbose)
      status_update(outfile,line,line+' '+fit_status)
    else:
      status_update(outfile,line,line+' not_done')
    if verbose==True:
      percent = 100.0*((i+1)-minlen)/(maxlen-minlen)
      print 'Process ',thread,' is ',percent,'% complete'

  linecache.clearcache()

def casu_solve(casuin,thresh=20,thread='',verbose=False):

  command = 'imcore '+casuin+' noconf outputcat'+str(thread)+'.fits 2 '+str(thresh)+' --filtfwhm=1 --noell'
  if verbose == True:
    print command
  os.system(command)

  # quick correction factor because the central wcs axis is not always pointed in the right place at the central distortion axis
  shift_wcs_axis(casuin,'outputcat'+str(thread)+'.fits',thresh=thresh)

  command = 'wcsfit '+casuin+' outputcat'+str(thread)+'.fits --site cds'
  if verbose == True:
    print command
  os.system(command)
  os.system('rm outputcat'+str(thread)+'.fits')

  status = 'ok'

  return status
