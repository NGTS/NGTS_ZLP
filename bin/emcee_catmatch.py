#!/usr/bin/env python
# -*- coding: utf-8 -*-

from NGTS_workpackage.catmatch import *
from NGTS_workpackage.emcee_tools import run_emcee
import numpy as np
from astropy.io import fits
import sys
import argparse
import multiprocessing as mp
import os
import emcee


def extract_coordinate_limits(filename):
    cat_names, RA_lims, DEC_lims = [], [], []
    with open(filename) as infile:
        for line in infile:
            vals = line.strip('\n').split(' ')
            cat_names += [vals[0]]
            RA_lims += [[float(vals[2]), float(vals[3])]]
            DEC_lims += [[float(vals[4]), float(vals[5])]]
    return cat_names, RA_lims, DEC_lims


def compute_separation(tel_ra, tel_dec, centre_ra, centre_dec):
    first = (tel_ra - centre_ra) * np.cos(np.radians(tel_dec))
    second = (tel_dec - centre_dec)
    sep = np.sqrt(first ** 2 + second ** 2)
    return sep


def find_nearest_catalogue(sourcedir, tel_ra, tel_dec):
    cat_names, RA_lims, DEC_lims = extract_coordinate_limits(
        os.path.join(sourcedir, 'index'))

    cen_RA = np.array([(f[0] + f[1]) / 2.0 for f in RA_lims])
    cen_DEC = np.array([(f[0] + f[1]) / 2.0 for f in DEC_lims])

    # Find the nearest catalogue
    sep = compute_separation(tel_ra, tel_dec, cen_RA, cen_DEC)

    cat_name = cat_names[np.argmin(sep)]

    return os.path.join(sourcedir, cat_name), RA_lims, DEC_lims


def extract_catalogue_data(cat_name, tel_ra, tel_dec, bright_mag=10, faint_mag=15):
    max_delta_coordinate = 2
    mag_key = 'Jmag'
    with fits.open(cat_name) as infile:
        data = infile[1].data

    ra = data['ra']
    dec = data['dec']
    mag = data[mag_key]

    ind = (mag < faint_mag) & (mag > bright_mag)
    ind &= ((ra > tel_ra - max_delta_coordinate) &
            (ra < tel_ra + max_delta_coordinate))
    ind &= ((dec > tel_dec - max_delta_coordinate) &
            (dec < tel_dec + max_delta_coordinate))

    return {
        'ra': ra[ind],
        'dec': dec[ind],
        mag_key: mag[ind],
    }


def main(args):
    nwalkers = args.nwalkers
    nruns = args.nruns
    nthreads = mp.cpu_count() if args.nthreads is None else args.nthreads
    burns = args.burns

    hdulist = fits.getheader(args.casuin)
    XVAL = hdulist['NAXIS1'] / 2
    YVAL = hdulist['NAXIS2'] / 2
    TEL_RA = hdulist['TEL_RA']
    TEL_DEC = hdulist['TEL_DEC']

    cat_name, RA_lims, DEC_lims = find_nearest_catalogue(
        args.catsrc, TEL_RA, TEL_DEC)
    cat = extract_catalogue_data(cat_name, TEL_RA, TEL_DEC)

    with fits.open(args.mycatname) as mycatt:
        catdata = mycatt[1].data

    # Only choose the bright objects
    flux_key = 'Aper_flux_3'
    ind = (catdata[flux_key] > 1E3) & (catdata[flux_key] < 1E5)
    catdata = catdata[ind]

    mycat = {flux_key: catdata[flux_key]}
    my_X = catdata['x_coordinate']
    my_Y = catdata['y_coordinate']
    pix_coords = [[my_X[i], my_Y[i]] for i in range(0, len(my_X))]

    # 7th order fit
    dicty = {
        'CRPIX1': 1.03259815e+03,
        'CRPIX2': 9.65505144e+02,
        'CD1_1': 1.41142333e-03,
        'CD2_2': 1.41109400e-03,
        'CD1_2': -1.89116218e-06,
        'CD2_1': 1.53342393e-06,
        'PV2_1': 1.0,
        'PV2_3': 8.68515702e+00,
        'PV2_5': 2.70336203e+02,
        'PV2_7': 1.37726138e+04,
        'RA_s': -0.45807896397,
        'DEC_s': 0.48575139999,
        'CTYPE1': 'RA---ZPN',
        'CTYPE2': 'DEC--ZPN',
    }
    name_list = ['CRPIX1', 'CRPIX2', 'CD1_1', 'CD2_2', 'CD1_2',
                 'CD2_1', 'PV2_3', 'PV2_5', 'PV2_7', 'RA_s', 'DEC_s']

    dicty['CRVAL1'] = TEL_RA + dicty['RA_s']
    dicty['CRVAL2'] = TEL_DEC + dicty['DEC_s']

    cen = [[dicty['CRPIX1'], dicty['CRPIX2']]]

    old_world = load_wcs_from_keywords(hdulist, cen)

    TEL_RA = hdulist['TEL_RA']
    TEL_DEC = hdulist['TEL_DEC']

    dicty['RA_s'] = (old_world[0][0] - TEL_RA)
    dicty['DEC_s'] = (old_world[0][1] - TEL_DEC)

    apply_correct(dicty, args.casuin, TEL_RA, TEL_DEC)

    prior = [hdulist[name] for name in name_list]

    prior[-1] = hdulist['CRVAL2'] - TEL_DEC
    prior[-2] = hdulist['CRVAL1'] - TEL_RA

    start_size = np.array([1e-3] * len(prior))

    prob_args = [args.casuin, mycat, cat, XVAL, YVAL, TEL_RA, TEL_DEC,
                 RA_lims, DEC_lims, my_X, my_Y, pix_coords, name_list, dicty]

    ndim = len(prior)
    p0 = []
    for i in range(0, nwalkers):
        shuffle = (10 ** (start_size * (np.random.rand(ndim) - 0.5)))
        p0 += [list(shuffle * prior)]

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=prob_args,
            threads=nthreads)
    sampler.run_mcmc(p0, nruns)


def lnprob(x, casuin, mycat, cat,
           XVAL, YVAL, TEL_RA, TEL_DEC, RA_lims, DEC_lims,
           my_X, my_Y, pix_coords, name_list, dicty):

    for i in range(0, len(x)):
        dicty[name_list[i]] = x[i]

    dicty['CRVAL1'] = TEL_RA + dicty['RA_s']
    dicty['CRVAL2'] = TEL_DEC + dicty['DEC_s']

    try:
        rms = fit_shift_wcs_axis(dicty, casuin, mycat, cat, XVAL, YVAL,
                                 TEL_RA, TEL_DEC, RA_lims, DEC_lims, my_X, my_Y, pix_coords)
    except:
        return -np.inf

    likelihood = -2000 * ((np.median(rms)) ** 2.0)

    lp = lnprior(dicty, rms)
    if not np.isfinite(lp):
        return -np.inf

    lnoutput = lp + likelihood

    print lnoutput, np.median(rms)

    return lnoutput


def lnprior(dicty, rms):

    if ((len(rms) > 1000) and
            (0.0012 < dicty['CD1_1'] < 0.0025) and
            (0.0012 < dicty['CD1_1'] < 0.0025) and
            ((abs(dicty['CD2_1']) < 1e-4)) and
            (abs(dicty['CD1_2']) < 1e-4) and
            (abs(dicty['RA_s']) < 1.0) and
            (abs(dicty['DEC_s'] < 1.0))):
        return 0.0

    return -np.inf

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('casuin')
    parser.add_argument('mycatname')
    parser.add_argument('chain_name')
    parser.add_argument('catsrc')
    parser.add_argument('--nwalkers', required=False, type=int, default=1000)
    parser.add_argument('--nruns', required=False, type=int, default=1E7)
    parser.add_argument('--nthreads', required=False, type=int,
                        default=None)
    parser.add_argument('--burns', required=False, type=int, default=0)
    main(parser.parse_args())
