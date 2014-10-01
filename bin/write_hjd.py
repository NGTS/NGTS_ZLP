#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import fitsio
from numpy import cos, sin
from astropy import units as u
from astropy.units import cds
from NGTS_workpackage.hjd_correction import append_hjd_column

def main(args):
    append_hjd_column(args.filename, column_name=args.column_name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('-c', '--column_name', help='Column name to insert',
                        default='hjd', required=False)
    main(parser.parse_args())
