#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
 
Zero Level Pipeline apperture photometry 

Usage: 
  ZLP_app_photom [-h] [--verbose] [--nproc=NPROC] [--appsize=APPSIZE] [--s_thresh=S_THRESH] (--confmap=CONFMAP) (--catfile=CATFILE) (--filelist=FILELIST | INPUT ...) [--outlist=OUTLIST]

Options:
  -h --help              Show help text
  --verbose              Print more text
  --catfile=CATFILE      The catalog file for this field
  --confmap=CONFMAP      The confidence map (callibration product)
  --filelist=FILELIST    Specify a filelist to use instead of command line
  --outlist=OUTLIST      Specify the name of the list of completed files
  --nproc=NPROC          Enable multithreading if you're analysing a lot of files at once [default: 1]
  --appsize=APPSIZE      The radius of the apperture you wish to use in the photometry stage [default: 2]
  --s_thresh=S_THRESH    The detection threshold to use when WCS solving images - typically higher than when doing actual photometry [default: 20]

This is the apperture photometry portion of the pipeline. It can be driven either in a list mode
or on a single file
 
"""
from docopt import docopt
import sys
import linecache
from numpy import *
import threading
from os.path import isfile, join
from NGTS_workpackage import *

argv = docopt(__doc__)

#if you don't provide an outlist name i'll assume you just want to add _phot to the end
if not argv['--outlist']:
  argv['--outlist'] = argv['--filelist'] + '_phot'

outfile = open(argv['--outlist'],'w')
outfile.close()


if argv['--filelist']:
  m_solve_images(argv['--filelist'],argv['--outlist'],nproc=int(argv['--nproc']),thresh=int(argv['--s_thresh']),verbose=argv['--verbose'])
  m_wcs_photom(argv['--filelist'],argv['--outlist'],int(argv['--appsize']),argv['--confmap'],argv['--catfile'],nproc=int(argv['--nproc']),verbose=argv['--verbose'])
  m_condense_data(argv['--filelist'],int(argv['--nproc']),int(argv['--appsize']),verbose=argv['--verbose'])

if argv['INPUT']:
  for filename in argv['INPUT']:
    casu_solve(filename,argv['--s_thresh'],verbose=argv['--verbose'])
    casu_photom(filename,argv['--confmap'],argv['--catfile'],argv['--appsize'],verbose=argv['--verbose'])
