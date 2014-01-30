#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
 
Zero Level Pipeline catalog generation

Usage: 
  ZLP_create_cat [-h] [--verbose] [--nproc=NPROC] [--nfiles=NFILES] [--s_thresh=S_THRESH] [--c_thresh=C_THRESH] [--stacklist=STACKLIST] (--confmap=CONFMAP) (--filelist=FILELIST) [--outname=OUTNAME]

Options:
  -h --help              Show help text
  --verbose              Print more text
  --confmap=CONFMAP      The confidence map (callibration product)
  --filelist=FILELIST    Specify a filelist to use instead of command line
  --outname=OUTNAME      Specify the name of the output catalog [default: catfile.fits]
  --stacklist=STACKLIST  The name of the file that stores the names of the images used in the stack [default: stackfilelist]
  --nproc=NPROC          Enable multithreading if you're analysing a lot of files at once [default: 16]
  --c_thresh=C_THRESH    The detection threshold to use when defining the input [default: 2]
  --s_thresh=S_THRESH    The detection threshold to use when WCS solving images - typically higher than when doing actual photometry [default: 20]
  --nfiles=NFILES        Maximum number of files to use in the stack [default: 16]

This is the catalog generation tool, requires a filelist input. need to work on being selective on the files used in input.
 
"""

from docopt import docopt
import sys
import linecache
from numpy import *
import threading
from os.path import isfile, join
from NGTS_workpackage import *

argv = docopt(__doc__)

if argv['--verbose'] == True:
  print 'Creating source catalogue from first ',argv['--nfiles'],' images...'

# for now we just pick the first n images, in future we might want to be more selective.

tmp = open('tmp','w')
i = 0
for line in open(argv['--filelist'],'r'):
    cal_stat = line.strip('\n').split(' ')[1]
    if cal_stat == 'ok':
      image = line.split(' ')[0]
      if i < int(argv['--nfiles']):
        tmp.write(line)
      i+=1
tmp.close()

m_solve_images('tmp','tmp',thresh=argv['--s_thresh'],nproc=int(argv['--nproc']),verbose=argv['--verbose'])

stacklist = open(argv['--stacklist'],'w')
for line in open('tmp'):
  image = line.strip('\n').split(' ')[0]
  status_check = line.strip('\n').split(' ')[1:]
  if all([status == 'ok' for status in status_check]):
    stacklist.write(image+'\n')
stacklist.close()

os.system('rm tmp')

command = 'casu_imstack @'+argv['--stacklist']+' '+argv['--confmap']+' "" outstack.fits outstackconf.fits'
if argv['--verbose'] == True:
  print(command)
os.system(command)

command = 'imcore outstack.fits outstackconf.fits '+argv['--outname']+' 2 '+argv['--c_thresh'] + ' --filtfwhm=1'
if argv['--verbose'] == True:
  print(command)
os.system(command)  

if argv['--verbose'] == True:
  print 'Catalogue complete'