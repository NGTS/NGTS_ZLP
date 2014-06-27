#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
 
Zero Level Pipeline apperture photometry 

Usage: 
  ZLP_app_photom [options] (-c <CONFMAP> | --confmap <CONFMAP>) (-C <CATFILE> | --catfile <CATFILE>) (-f <FILELIST> | --filelist <FILELIST> | INPUT ...)

Options:
  -h --help              Show help text
  --verbose              Print more text
  --outlist=OUTLIST      Specify the name of the list of completed files
  --nproc=NPROC          Enable multithreading if you're analysing a lot of files at once [default: 1]
  --apsize=APSIZE        The radius of the apperture you wish to use in the photometry stage [default: 2]
  --s_thresh=S_THRESH    The detection threshold to use when WCS solving images - typically higher than when doing actual photometry [default: 7]
  --catsrc=CATSRC        What catalogue to use during catalog matching [default: viz2mass]
  --catpath=CATPATH      If you're using a local catalog for cat matching, where is it? [default: False]
  --outdir=OUTDIR        Where you would like the result files to go [default: ./]

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

filelist = argv['<FILELIST>']
m_condense_data(filelist,int(argv['--nproc']),int(argv['--apsize']),verbose=argv['--verbose'],outdir=argv['--outdir'])