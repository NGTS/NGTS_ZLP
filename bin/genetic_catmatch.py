#!/usr/bin/env python
# -*- coding: utf-8 -*-

from NGTS_workpackage.catmatch import *
from NGTS_workpackage.emcee_tools import run_emcee
from deap import base, creator, tools, algorithms
import numpy as np
import random
import fitsio
import astropy.io.fits as pf
import sys
import argparse
import multiprocessing as mp

def main(args):

    resume = False
    nwalkers = 1000
    nruns = 1e7
    nthreads = mp.cpu_count()
    burns = 0
    start_size = 1e-2

    hdulist = fitsio.read_header(args.casuin.rstrip('.fits')+'.new')
    XVAL = hdulist['NAXIS1'] / 2
    YVAL = hdulist['NAXIS2'] / 2
    TEL_RA = hdulist['TEL_RA']
    TEL_DEC = hdulist['TEL_DEC']

    cat_names = []
    RA_lims = []
    DEC_lims = []
    for line in open(args.catsrc + '/index'):
        vals = line.strip().split(' ')
        cat_names += [vals[0]]
        RA_lims += [[float(vals[2]), float(vals[3])]]
        DEC_lims += [[float(vals[4]), float(vals[5])]]

    cen_RA = np.array([(f[0] + f[1]) / 2.0 for f in RA_lims])
    cen_DEC = np.array([(f[0] + f[1]) / 2.0 for f in DEC_lims])

    sep = (((TEL_RA - cen_RA) * (np.cos(TEL_DEC * np.pi / 180.0))) ** 2.0 +
           (TEL_DEC - cen_DEC) ** 2.0) ** 0.5

    cat_name = cat_names[np.argmin(sep)]

    with pf.open(args.catsrc + '/' + cat_name) as catd:
        catt = catd[1].data.copy()
    cat = {'ra': catt['ra'], 'dec': catt['dec'], 'Jmag': catt['Jmag']}

    with fitsio.FITS(args.mycatname) as mycatt:
        mycat = {'Aper_flux_3': mycatt[1]['Aper_flux_3'][:]}
        my_X = mycatt[1]['x_coordinate'][:]
        my_Y = mycatt[1]['y_coordinate'][:]
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
        'CTYPE2': 'DEC--ZPN'
    }
    name_list = ['CRPIX1', 'CRPIX2', 'CD1_1', 'CD2_2', 'CD1_2', 'CD2_1',
                 'PV2_3', 'PV2_5', 'PV2_7', 'RA_s', 'DEC_s']

    dicty['CRVAL1'] = TEL_RA + dicty['RA_s']
    dicty['CRVAL2'] = TEL_DEC + dicty['DEC_s']

    cen = [[dicty['CRPIX1'], dicty['CRPIX2']]]

    old_world = load_wcs_from_keywords(hdulist, cen)

    TEL_RA = hdulist['TEL_RA']
    TEL_DEC = hdulist['TEL_DEC']

    dicty['RA_s'] = (old_world[0][0] - TEL_RA)
    dicty['DEC_s'] = (old_world[0][1] - TEL_DEC)

    apply_correct_old(dicty, args.casuin, TEL_RA, TEL_DEC)

    prior = []
    for i in name_list:
        try:
            prior += [hdulist[i]]
        except:
            prior += [dicty[i]]


    prior[-1] = hdulist['CRVAL2'] - TEL_DEC
    prior[-2] = hdulist['CRVAL1'] - TEL_RA

    prior = [1036.403014148268, 1094.799578289627, -0.0014024689793905526, -0.001400992952741607, 2.434072414406271e-05, -2.3811448821252857e-05, 10.077655889076954, 278.6726383959758, 13717.617866838911, -0.002732639996949601, -0.11274463084725035]

    prior = [1036.7659160556957, 1096.8252100341174, -0.001401848440110451, -0.0014013722072105585, 2.3887290238287583e-05, -2.42612367087932e-05, 6.94957709088958, 443.9094722656311, 15713.234246634003, -0.003301766749383896, -0.1156524974929483]

    start_size = np.array([1e-3] * len(prior))

    for i in [2, 3, 4, 5]:
        start_size[i] = 1e-3

    for i in [6, 7, 8, 9]:
        start_size[i] = 1e-3

    my_args = [args.casuin, mycat, cat, XVAL, YVAL, TEL_RA, TEL_DEC, RA_lims,DEC_lims, my_X, my_Y, pix_coords, name_list, dicty]

    creator.create("FitnessMax", base.Fitness, weights=(-1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)

    toolbox = base.Toolbox()

    # this method sets up the initial population with a guassian scatter around the best guess values with
    # standard deviation of 0.1
    toolbox.register("my_prior", scatter_prior,prior,0.01)

    toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.my_prior)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    # tell the code where the fitness function is
    toolbox.register("evaluate", lnprob, args.casuin, mycat, cat, XVAL, YVAL, TEL_RA, TEL_DEC, RA_lims,DEC_lims, my_X, my_Y, pix_coords, name_list, dicty)

#    details of evolution - requires tweaking!
    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("mutate", scatter_individual, sigma=0.01, indp=1.0,scale_jump_p=0.1,tweak_jump_p=0.2)
    toolbox.register("select", tools.selTournament, tournsize=3)

    pop = toolbox.population(n=1000)
    hof = tools.HallOfFame(1)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("min", np.min)
    stats.register("max", np.max)
    
#    old_best = [872.7225757887486, 1113.8052862722673, 0.001235502346442564, 0.00071341762270345, -1.8620475874018094e-06, 9.392021451322169e-07, 8.696428851660816, 281.6029470630021, 19648.54124845112, 0.865776751183994, -2.023251331235531]

    for i in range(0,10000):
        #required to be in loop to plug a memory leak... non optimal.
        pool = mp.Pool(nthreads)
        toolbox.register("map", pool.map)

        pop, logbook = algorithms.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.5, ngen=10, stats=stats, halloffame=hof, verbose=True)
        pool.close()

        x = hof[0]
        print x

        for i in range(0, len(x)):
            dicty[name_list[i]] = x[i]

        dicty['CRVAL1'] = TEL_RA + dicty['RA_s']
        dicty['CRVAL2'] = TEL_DEC + dicty['DEC_s']

        rms = fit_shift_wcs_axis(dicty, args.casuin, mycat, cat, XVAL, YVAL, TEL_RA,
                                 TEL_DEC, RA_lims, DEC_lims, my_X, my_Y,
                                 pix_coords,update=True)
        print np.median(rms)


def scatter_prior(prior,sigma):
    scattered_prior = []
    for i in range(0,len(prior)):
        new_val = np.random.normal(prior[i],np.abs(sigma*prior[i]))
        scattered_prior += [new_val]
    return scattered_prior

def scatter_individual(individual,sigma=0.01,indp=0.01,scale_jump_p=0.01,tweak_jump_p=0.01):
    scattered_prior = []
    for i in range(0,len(individual)):
        if np.random.uniform(0,1) < indp:
            if np.random.uniform(0,1) < scale_jump_p:
                individual[i] = 10**(np.random.normal(np.log10(individual[i]),np.abs(sigma*np.log10(individual[i]))))
            elif np.random.uniform(0,1) < tweak_jump_p:
                individual[i] = np.random.normal(individual[i],np.abs(0.001*sigma*individual[i]))
            else:
                individual[i] = np.random.normal(individual[i],np.abs(sigma*individual[i]))
    return individual,

def lnprob(casuin, mycat, cat, XVAL, YVAL, TEL_RA, TEL_DEC, RA_lims, DEC_lims, my_X, my_Y, pix_coords, name_list, dicty,x):

    for i in range(0, len(x)):
        dicty[name_list[i]] = x[i]

    dicty['CRVAL1'] = TEL_RA + dicty['RA_s']
    dicty['CRVAL2'] = TEL_DEC + dicty['DEC_s']

    try:
        rms = fit_shift_wcs_axis(dicty, casuin, mycat, cat, XVAL, YVAL, TEL_RA,
                                 TEL_DEC, RA_lims, DEC_lims, my_X, my_Y,
                                 pix_coords)
    except:
        return np.inf,

    likelyhood = -2000 * ((np.median(rms)) ** 2.0)

    im_header = fitsio.read_header(casuin)
    lp = lnprior(dicty, rms, im_header)
    lnoutput = lp + likelyhood

    if not np.isfinite(lp):
        return np.inf,

#    print lnoutput, np.median(rms)

#    return lnoutput,
    return np.median(rms),


def lnprior(dicty, rms,im_header):

#    if np.all((len(rms) > 1000) and (0.0012 < dicty['CD1_1'] < 0.0025) and
#              (0.0012 < dicty['CD1_1'] < 0.0025) and (
#                  (abs(dicty['CD2_1']) < 1e-4)) and
#              (abs(dicty['CD1_2']) < 1e-4) and (abs(dicty['RA_s']) < 1.0) and
#              (abs(dicty['DEC_s'] < 1.0))):
    left = [[1024, 0]]
    right = [[1024, 2048]]
    top = [[2048, 1024]]
    bottom = [[0, 1024]]

    w = wcs.WCS(dicty)
    wcs_left = w.wcs_pix2world(left,1)
    wcs_right = w.wcs_pix2world(right,1)
    wcs_top = w.wcs_pix2world(top,1)
    wcs_bottom = w.wcs_pix2world(bottom,1)

    x_sep = sky_sep(wcs_left[0],wcs_right[0])/3600
    y_sep = sky_sep(wcs_top[0],wcs_bottom[0])/3600

    if np.all((len(rms) > 100) & (x_sep > 2.75) & (y_sep > 2.75) & (x_sep < 2.95) & (y_sep < 2.95)):
        return 0.0
    return np.inf


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('casuin')
    parser.add_argument('mycatname')
    parser.add_argument('chain_name')
    parser.add_argument('catsrc')
    main(parser.parse_args())
