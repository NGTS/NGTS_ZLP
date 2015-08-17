# -*- coding: utf-8 -*-
"""
 
Create a vector plot of the astrometric errors. The errors are greatly exagerated for clarity.
Returns the median astrometric error for diagnostic use.

Usage:
  vector_plot.py [options] (<CATALOG_NAME>) (<IMAGE_NAME>) (<PLOT_NAME>) (<PLOT_NAME>) (<CAT>) (<CAT_SRC>)

Options:
  -h --help  Show help text

"""
import argparse
import numpy as np
from astropy import wcs
from astropy.io import fits
import fitsio
import os
from NGTS_workpackage.catmatch import load_wcs_from_keywords
from NGTS_workpackage.catmatch import calc_seps
from astropy.io import fits as pf


def wcsf_QCheck(catalog_name, image_name, plot_name, cat,
                upscale_factor=500,
                plot=True):

    cat_names = []
    RA_lims = []
    DEC_lims = []

    catsrc = cat.rstrip(cat.split('/')[-1])

    for line in open(catsrc+'/index'):
        vals = line.strip('\n').split(' ')
        cat_names += [vals[0]]
        RA_lims += [[float(vals[2]),float(vals[3])]]
        DEC_lims += [[float(vals[4]),float(vals[5])]]

    n = 0

    cat_name = cat_names[n]

    with pf.open(catsrc+'/'+cat_name) as catd:
      catt = catd[1].data.copy()
    cat = {'ra':catt['ra'],'dec':catt['dec'],'Jmag':catt['Jmag']}


    with pf.open(catalog_name) as mycatd:
        mycatt = mycatd[1].data.copy()
        my_X = mycatd[1].data['x_coordinate']
        my_Y = mycatd[1].data['y_coordinate']

    mycat = {'Aper_flux_3':mycatt['Aper_flux_3']}

    print 'about to plot'

    plot_dir = os.path.dirname(image_name)

    im_header = fitsio.read_header(image_name)
    pix_coords = [[my_X[i], my_Y[i]] for i in range(0, len(my_X))]
    world = load_wcs_from_keywords(im_header, pix_coords)

    xs, ys, RA_sep, DEC_sep, x_sep, y_sep, sep_list = calc_seps(
        mycat, cat, RA_lims, DEC_lims, world, my_X, my_Y, im_header)

    rms = np.median(sep_list)

    #from scipy.optimize import leastsq
    #true_cen_xs = iter_poly(xs,RA_sep)
    #true_cen_ys = iter_poly(ys,DEC_sep)
    #print true_cen_xs, true_cen_ys

    h = fitsio.read_header(image_name)
    CRPIX1 = h['CRPIX1']
    CRPIX2 = h['CRPIX2']

    cen_X = h['NAXIS1'] / 2.0
    cen_Y = h['NAXIS2'] / 2.0

    cen = [[cen_X, cen_Y]]

    cen_world = load_wcs_from_keywords(im_header, cen)

    print 'cen is', cen_world

    left = [[1024, 0]]
    right = [[1024, 2048]]
    top = [[2048, 1024]]
    bottom = [[0, 1024]]

    wcs_left = load_wcs_from_keywords(left,1)
    wcs_right = load_wcs_from_keywords(right,1)
    wcs_top = load_wcs_from_keywords(top,1)
    wcs_bottom = load_wcs_from_keywords(bottom,1)

    x_sep = sky_sep(wcs_left[0],wcs_right[0])/3600
    y_sep = sky_sep(wcs_top[0],wcs_bottom[0])/3600

    #  axis.plot(true_cen_xs,true_cen_ys,'go',markersize=10)


    print 'FOV:',x_sep,y_sep
    print 'found ',len(sep_list),' sources'


    if plot == True:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axis = plt.subplots(figsize=(11, 8))
        for i in range(0, len(xs)):
            axis.plot([xs[i], (xs[i] + x_sep[i] * upscale_factor)],
                      [ys[i], (ys[i] + y_sep[i] * upscale_factor)], 'k-')
            axis.plot((xs[i]), (ys[i]), 'ko', markersize=3)

        axis.plot(CRPIX1, CRPIX2, 'ro', markersize=10)

        axis.set_ylim([0, 2048])
        axis.set_xlim([0, 2048])
        axis.set_aspect('equal')
        axis.set_xlabel(r'X')
        axis.set_ylabel(r'Y')
        axis.set_title(r'RMS: ' + str(rms) + '"')
        fig.tight_layout()
        fig.savefig(os.path.join(plot_dir, plot_name), bbox_inches='tight')
        plt.close(fig)

    print 'done plotting'

    with fitsio.FITS(image_name, 'rw') as fits:
        fits[0].write_key('wcsf_ns', len(sep_list))
        fits[0].write_key('wcsf_RMS', rms)
        fits[0].write_key('wcsf_RA', cen_world[0][0])
        fits[0].write_key('wcsf_DEC', cen_world[0][1])

    return rms


def iter_poly(xs, RA_sep):

    RA_sep = np.array(RA_sep)
    xs = np.array(xs)
    prior = [0.001, 0.00000001, -0.000001, 0.000001, 500]
    xs_copy = xs.copy()
    RA_sep_copy = RA_sep.copy()

    for i in range(0, 50):
        x, success = leastsq(fit_polynomial, prior,
                             args=(xs_copy, RA_sep_copy),
                             epsfcn=0.01)
        poly = polynomial(x, xs_copy)
        residuals = poly - RA_sep_copy
        std_r = np.std(residuals)
        xs_copy = xs_copy[abs(residuals) < 3.0 * std_r]
        RA_sep_copy = RA_sep_copy[abs(residuals) < 3.0 * std_r]

    #poly = polynomial(x,xs)

    #plt.plot(xs,poly,'b.')

    #plt.plot(xs,RA_sep,'r.')
    #plt.show()

    return x[-1]


def fit_polynomial(x, y, data):

    model = polynomial(x, y)

    resid = data - model

    return resid / 0.01


def polynomial(x, y):

    model = y * 0
    for i in range(0, len(x) - 1.0):
        model += x[i] * (y - x[-1]) ** ((i * 2) + 1)
    return model


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("catalog_name")
    parser.add_argument('image_name')
    parser.add_argument('cat')
    parser.add_argument('plot_name')

    args = parser.parse_args()
    med_sep = wcsf_QCheck(args.catalog_name, args.image_name,
                               args.plot_name, args.cat)
    print 'the median sep is', med_sep, ' arcseconds'
