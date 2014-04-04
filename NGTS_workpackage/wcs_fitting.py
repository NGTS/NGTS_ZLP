# -*- coding: utf-8 -*-

import tempfile
from catmatch import shift_wcs_axis
from multiprocessing.dummy import Pool as ThreadPool
from functools import partial
import casutools
from util import validate_headers

def m_solve_images(filelist, outfile, nproc=None, thresh=20.0, verbose=False):
  infiles = []
  with open(filelist) as infile:
    for line in infile:
      image = line.strip('\n')
      status_checks = ['ok','ok']

      if all(status == 'ok' for status in status_checks):
        infiles.append(image)

  fn = partial(casu_solve, thresh=thresh, verbose=verbose)

  pool = ThreadPool(nproc)
  return pool.map(fn, infiles)

def casu_solve(casuin, thresh=20, verbose=False):

  validate_headers(casuin)

  with tempfile.NamedTemporaryFile(dir='.', suffix='.fits', prefix='catalogue.') as catfile:
    catfile_name = catfile.name

    casutools.imcore(casuin, catfile_name, threshold=thresh, verbose=verbose)
    catfile.seek(0)

    # quick correction factor because the central wcs axis is not always pointed in the right place at the central distortion axis
    try:
      shift_wcs_axis(casuin, catfile_name, thresh=thresh)
    except IOError:
      print "Performing initial fit"
      casutools.wcsfit(casuin, catfile_name, verbose=verbose)
      shift_wcs_axis(casuin, catfile_name, thresh=thresh)

    casutools.wcsfit(casuin, catfile_name, verbose=verbose)
    return 'ok'

