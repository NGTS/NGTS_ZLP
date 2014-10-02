from astropy import units as u
from astropy.units import cds
import numpy as np
import fitsio
import itertools


def compute_hjd(jd, ra, dec, sun_ra, sun_dec):
    r = u.AU.to(u.m)
    c = (1 * cds.c).si.value

    ra_rad, dec_rad, sun_ra_rad, sun_dec_rad = [
        np.radians(data).astype(np.float64) for data in [
            ra, dec, sun_ra, sun_dec
        ]
    ]

    first_term = np.sin(dec_rad) * np.sin(sun_dec_rad)
    second_term = np.cos(
        dec_rad) * np.cos(sun_dec_rad) * np.cos(ra_rad - sun_ra_rad)
    correction_seconds = (r / c) * (first_term + second_term)

    return -(correction_seconds / 86400.).astype(np.float64)


def compute_hjd_column(fname):
    with fitsio.FITS(fname) as infile:
        catalogue = infile[1]
        header = catalogue.read_header()

        ra, dec = [
            catalogue[key].read().astype(np.float64) for key in ['ra', 'dec']
        ]

    mjd, sun_ra, sun_dec = [header[key]
                            for key in ['mjd', 'sun_ra', 'sun_dec']]

    return compute_hjd(mjd, ra, dec, sun_ra, sun_dec)


def append_hjd_column(fname, column_name='hjd'):
    hjd_data = compute_hjd_column(fname)

    with fitsio.FITS(fname, 'rw') as outfile:
        outfile[1].insert_column(column_name, hjd_data)
