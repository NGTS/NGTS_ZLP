#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Zero Level Pipeline catalog generation

Usage:
  ZLP_create_cat [options] (--incat <INPUT>)

Options:
  -h --help                                 Show help text
  --verbose                                 Print more text
  --outname=OUTNAME                         Specify the name of the output catalog [default: maglim.fits]
  --max=MAX_MAG                             Maximum mag [default: 13]
  --min=MIN_MAG                             Minimum mag [default: 8]

This is the catalog generation tool, requires a filelist input. need to work on being selective on the files used in input.

"""

from docopt import docopt
from astropy.io import fits
from numpy import *

def main(argv):
    if argv['--verbose'] == True:
        print 'Creating mag limited finder catalogue from first {} images...'

    outfile = open(argv['--outname'],'w')

    max_mag = float(argv['--max'])
    min_mag = float(argv['--min'])

    ra = array([])
    dec = array([])

    for line in open(argv['<INPUT>'],'r'):
	try:
  	  mag = float(line.split('\t')[5])         
	  if ((mag > min_mag) and (mag < max_mag)):
	    ra = append(ra,float(line.split('\t')[0]))
            dec = append(dec,float(line.split('\t')[1]))
	except:
	  print 'failure'

    if argv['--verbose'] == True:
        print 'Catalogue complete'

    print ra
    print dec

if __name__ == '__main__':
    main(docopt(__doc__))
