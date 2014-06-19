# -*- coding: utf-8 -*-
from astropy import wcs
from astropy.io import fits
import numpy as np

def thread_alloc(nfiles, nproc):

  chunks = int(nfiles/nproc)
  leftovers = nfiles - chunks*nproc

  starts = []
  ends = []

  for i in range(1,nproc+1):
    starts += [(i-1)*chunks]
    ends += [i*chunks]

  ends = np.array(ends)+1
  starts = np.array(starts)+1

#  adds the leftovers to the processors evenly
  procn = 0
  while ends[-1] < (nfiles+1) and (procn<nproc):
    starts[procn+1:] += 1
    ends[procn:] += 1
    procn += 1

  return starts, ends

def load_wcs_from_file(filename,pixcrd):
# Load the WCS information from a fits header, and use it
# to convert pixel coordinates to world coordinates.
   # Load the FITS hdulist using astropy.io.fits

    #hdulist = fits.open(filename)
    #fheader = hdulist[0].header
    fheader = fitsio.read_header(filename, 0)

    print(fheader)

    # Parse the WCS keywords in the primary HDU
    w = wcs.WCS(fheader)

    # Convert pixel coordinates to world coordinates
    # The second argument is "origin" -- in this case we're declaring we
    # have 1-based (Fortran-like) coordinates.
    world = w.wcs_pix2world(pixcrd, 1)

    return world

def load_wcs_from_keywords(fheader,pixcrd):
# Load the WCS information from a fits header, and use it
# to convert pixel coordinates to world coordinates.
   # Load the FITS hdulist using astropy.io.fits

    #hdulist = fits.open(filename)
    #fheader = hdulist[0].header

    # Parse the WCS keywords in the primary HDU
    w = wcs.WCS(fheader)

    # Convert pixel coordinates to world coordinates
    # The second argument is "origin" -- in this case we're declaring we
    # have 1-based (Fortran-like) coordinates.
    world = w.wcs_pix2world(pixcrd, 1)

    return world

def status_update(file_path, pattern, subst):
  i = 0
  while i < 100:
    try:
     s_update(file_path, pattern, subst)
     return
    except:
      print('attempt ',str(i),'!')
      time.sleep(0.01)
      i += 1
      
def validate_headers(image):
  fits.setval(image,'CTYPE1',value='RA---ZPN')
  fits.setval(image,'CTYPE2',value='DEC--ZPN')

def s_update(file_path, pattern, subst):

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
  os.remove(file_path)
  #Move new file
  shutil.move(abs_path, file_path)
