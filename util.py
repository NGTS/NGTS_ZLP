# -*- coding: utf-8 -*-
from __future__ import division, print_function
import os
from numpy import *
def thread_alloc(nfiles, nproc):

  chunks = int(nfiles/nproc)
  leftovers = nfiles - chunks*nproc

  starts = []
  ends = []

  for i in range(1,nproc+1):
    starts += [(i-1)*chunks]
    ends += [i*chunks]

  ends = array(ends)+1
  starts = array(starts)+1

#  adds the leftovers to the processors evenly
  procn = 0
  while ends[-1] < (nfiles+1) and (procn<nproc):
    starts[procn+1:] += 1
    ends[procn:] += 1
    procn += 1

  return starts, ends

def genfilelist(directory_list,name):

# generates a sanitised file list suitible for use with the rest of the toolset

  fout = open(name,'w')
  i = 0

  for directory in directory_list:
    os.system('ls '+ directory +'> tmp')  
    for line in open('tmp'):
      if '.fits\n' in line:
	fout.write(directory + line)
	i += 1

  fout.close

  os.system('rm tmp')

  return i

def load_wcs_from_file(filename,pixcrd):
# Load the WCS information from a fits header, and use it
# to convert pixel coordinates to world coordinates.
    import numpy
    from astropy import wcs
    from astropy.io import fits
    import sys
   # Load the FITS hdulist using astropy.io.fits

    hdulist = fits.open(filename)

    # Parse the WCS keywords in the primary HDU
    w = wcs.WCS(hdulist[0].header)

    # Convert pixel coordinates to world coordinates
    # The second argument is "origin" -- in this case we're declaring we
    # have 1-based (Fortran-like) coordinates.
    world = w.wcs_pix2world(pixcrd, 1)

    return world

def status_update(file_path, pattern, subst):
  from tempfile import mkstemp
  from shutil import move
  from os import remove, close

  #Create temp file
  fh, abs_path = mkstemp()
  new_file = open(abs_path,'w')
  old_file = open(file_path)
  for line in old_file:
      new_file.write(line.replace(pattern, subst))
  #close temp file
  new_file.close()
  close(fh)
  old_file.close()
  #Remove original file
  remove(file_path)
  #Move new file
  move(abs_path, file_path)


def callsysrem(dir,repeats,modified = False):
  # theres an option to use simons modified sysrem now

  command = 'Sysrem'
  if modified == True:
    command = '/home/astro/phrfbf/work/PostDoc/NGTS/sysrem-altering'

  os.system(command)
  for i in range(repeats):
    os.system('mv output.fits tamout_'+str(i)+'.fits')
    os.system('mv tamout.fits output.fits')
    os.system(command)
  os.system('mv output.fits tamout_'+str(repeats)+'.fits')
  os.system('mv tamout.fits tamout_'+str(repeats + 1)+'.fits')    
