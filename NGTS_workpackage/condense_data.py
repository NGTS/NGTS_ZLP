# -*- coding: utf-8 -*-
from astropy.io import fits as pf
import os
import numpy as np
import linecache
import threading
import multiprocessing
#import scipy.optimize as opt
from os import listdir
from os.path import isfile, join
from util import thread_alloc
from numpy import *
import fitsio
import sys
from functools import partial
from pipeutils import NullPool
import os


def m_condense_data(filelist, nproc, appsize, verbose=False, outdir='./'):
    with open(filelist) as infile:
        files = sorted([line.strip() for line in infile])

    condense_data(files, os.path.join(outdir, 'output.fits'), appsize, verbose)


def condense_data(files, outfile_name, appsize, verbose):

    flux = []
    flux_err = []

    flux_grid = []
    flux_err_grid = []

    Skylev = []
    Skyrms = []
    sky = []
    xpos = []
    ypos = []
    ALT = []
    AZ = []
    TEL_RA = []
    TEL_DEC = []
    time = []
    centerra = []
    centerdec = []
    meanbias = []
    stdbias = []
    skylevel = []
    std_image = []
    T = []
    coolstat = []
    ADU_DEV = []
    fwhm = []

    fwhma_1 = []
    fwhmb_1 = []
    fwhmt_1 = []

    fwhma_2 = []
    fwhmb_2 = []
    fwhmt_2 = []

    fwhma_3 = []
    fwhmb_3 = []
    fwhmt_3 = []

    fwhma_4 = []
    fwhmb_4 = []
    fwhmt_4 = []

    fwhma_5 = []
    fwhmb_5 = []
    fwhmt_5 = []

    fwhma_6 = []
    fwhmb_6 = []
    fwhmt_6 = []

    fwhma_7 = []
    fwhmb_7 = []
    fwhmt_7 = []

    fwhma_8 = []
    fwhmb_8 = []
    fwhmt_8 = []

    fwhma_9 = []
    fwhmb_9 = []
    fwhmt_9 = []

    SKY_MED = []
    CLOUDS = []
    SHIFT = []
    exposure = []
    seeing = []
    imid = []
    hjd_hist = []
    airmass = []

    #the .copy() addition is necessary when dealing with long filelists - by default python list optimization keeps all files open otherwise,
    #leading to a crash from too many open files

    npix = np.pi * appsize ** 2.0

    first_frame = True

    for fname in files:
        try:
            with pf.open(fname+'.phot') as photdata:
                ambient = photdata[1].header.get('WXTEMP', 30.0)
                cloud_status = photdata[1].header['CLOUD_S']
                fwhm_frame = photdata[1].header['FWHM']
                psf_a_1 = photdata[1].header['PSF_a_1']
                psf_b_1 = photdata[1].header['PSF_b_1']
                psf_t_1 = photdata[1].header['PSF_t_1']

                psf_a_2 = photdata[1].header['PSF_a_2']
                psf_b_2 = photdata[1].header['PSF_b_2']
                psf_t_2 = photdata[1].header['PSF_t_2']

                psf_a_3 = photdata[1].header['PSF_a_3']
                psf_b_3 = photdata[1].header['PSF_b_3']
                psf_t_3 = photdata[1].header['PSF_t_3']

                psf_a_4 = photdata[1].header['PSF_a_4']
                psf_b_4 = photdata[1].header['PSF_b_4']
                psf_t_4 = photdata[1].header['PSF_t_4']

                psf_a_5 = photdata[1].header['PSF_a_5']
                psf_b_5 = photdata[1].header['PSF_b_5']
                psf_t_5 = photdata[1].header['PSF_t_5']

                psf_a_6 = photdata[1].header['PSF_a_6']
                psf_b_6 = photdata[1].header['PSF_b_6']
                psf_t_6 = photdata[1].header['PSF_t_6']

                psf_a_7 = photdata[1].header['PSF_a_7']
                psf_b_7 = photdata[1].header['PSF_b_7']
                psf_t_7 = photdata[1].header['PSF_t_7']

                psf_a_8 = photdata[1].header['PSF_a_8']
                psf_b_8 = photdata[1].header['PSF_b_8']
                psf_t_8 = photdata[1].header['PSF_t_8']

                psf_a_9 = photdata[1].header['PSF_a_9']
                psf_b_9 = photdata[1].header['PSF_b_9']
                psf_t_9 = photdata[1].header['PSF_t_9']

                frame_shift = photdata[1].header['SKY_MOVE']
                seeing_frame = photdata[1].header['SEEING']

                imid += [photdata[1].header['IMAGE_ID']]
                SHIFT += [frame_shift]

                CLOUDS += [cloud_status]
                SKY_MED += [photdata[1].header['SKYLEVEL']]
                ALT +=[photdata[1].header['TEL_ALT']]
                AZ +=[photdata[1].header['TEL_AZ']]
                TEL_RA +=[photdata[1].header['TEL_RA']]
                TEL_DEC +=[photdata[1].header['TEL_DEC']]
                airmass += [photdata[1].header.get('AIRMASS')]
                exposure += [photdata[1].header['EXPOSURE']]
                try:
                    ADU_DEV +=[photdata[1].header['ADU_DEV']]
                    skylevel +=[photdata[1].header['skylevel']]
                    meanbias += [photdata[1].header['BIASMEAN']]
                except:
                    imagedata = pf.open(line)
                    ADU_DEV +=[std(imagedata[0].data)]
                    skylevel +=[median(imagedata[0].data)]      
                    biasstrip = append(imagedata[0].data[:,:20],imagedata[0].data[:,-20:])
                    meanbias += [mean(biasstrip)]
                centerra +=[photdata[1].header['WCSF_RA']]
                centerdec +=[photdata[1].header['WCSF_DEC']]
                # correcting for airmass - the newer fits files have an airmass term, so just use that instead perhaps
                sky += [photdata[1].data['Sky_level'].copy()]
                frame_xpos = photdata[1].data['X_coordinate'].copy()
                frame_ypos = photdata[1].data['Y_coordinate'].copy()
                xpos += [frame_xpos]
                ypos += [frame_ypos]
                utc = photdata[1].header['OBSSTART'].split('T')
                yr, month, day = utc[0].split('-')
                hr, min, sec = utc[1].split(':')
                fwhm += [fwhm_frame]
                fwhma_1 += [psf_a_1]
                fwhmb_1 += [psf_b_1]
                fwhmt_1 += [psf_t_1]

                fwhma_2 += [psf_a_2]
                fwhmb_2 += [psf_b_2]
                fwhmt_2 += [psf_t_2]

                fwhma_3 += [psf_a_3]
                fwhmb_3 += [psf_b_3]
                fwhmt_3 += [psf_t_3]

                fwhma_4 += [psf_a_4]
                fwhmb_4 += [psf_b_4]
                fwhmt_4 += [psf_t_4]

                fwhma_5 += [psf_a_5]
                fwhmb_5 += [psf_b_5]
                fwhmt_5 += [psf_t_5]

                fwhma_6 += [psf_a_6]
                fwhmb_6 += [psf_b_6]
                fwhmt_6 += [psf_t_6]

                fwhma_7 += [psf_a_7]
                fwhmb_7 += [psf_b_7]
                fwhmt_7 += [psf_t_7]

                fwhma_8 += [psf_a_8]
                fwhmb_8 += [psf_b_8]
                fwhmt_8 += [psf_t_8]

                fwhma_9 += [psf_a_9]
                fwhmb_9 += [psf_b_9]
                fwhmt_9 += [psf_t_9]

                seeing += [seeing_frame]
                rawflux = photdata[1].data['Aper_flux_3'].copy()
                Skylev += [photdata[1].data['Sky_level'].copy()]
                Skyrms += [photdata[1].data['Sky_rms'].copy()]
                flux += [rawflux]
                flux_err += [photdata[1].data['Aper_flux_3_err'].copy()]

                all_appertures = []
                all_err_apps = []
                for aper in [2, 4, 5, 6, 7]:
                    flux_key = 'Aper_flux_{}'.format(aper)
                    fluxerr_key = '{}_err'.format(flux_key)
                    rawflux = photdata[1].data[flux_key].copy()
                    all_appertures += [rawflux]
                    all_err_apps += [photdata[1].data[fluxerr_key].copy()]

                flux_grid += [all_appertures]
                flux_err_grid += [all_err_apps]

                T +=[photdata[1].header['CCDTEMP']]
                coolstat +=[photdata[1].header['COOLSTAT']]
                mjd = photdata[1].header['MJD']
                hjd_hist += [mjd + photdata[1].data['hjd_correction'].copy()]
                time +=[[mjd]*len(flux[0])]
                if verbose == True:
                    print shape(time), line.split(' ')[0]+'.phot', thread_no

        except Exception as err:
            sys.stderr.write('Error analysing file {}, original error: {}\n'.format(
                image + '.phot', str(err)))

    # generate time of mid exposure array
    tarray = np.array(time)
    tmid = tarray[:, 0] + ((array(exposure) / 2.0) / (3600.0 * 24.0))

    flux_grid = np.array(flux_grid)
    error_grid = np.array(flux_err_grid)

    #get some data that should be common to all frames
    with pf.open(files[0] + '.phot') as photdata:
        RA = photdata[1].data['RA']
        DEC = photdata[1].data['DEC']

    zeros = tmid * 0

    objid = np.arange(len(flux[0])) + 1

    fluxarray = np.array(flux).T
    flux_err_array = np.array(flux_err).T
    hjdarray = np.array(hjd_hist).T

    meanflux = []
    for index in range(0, len(fluxarray[:, 0])):
        meanflux += [np.mean(fluxarray[index,:])]

    meanarray = np.array(meanflux)

    npts = np.ones_like(meanarray) * len(files)

    c1 = pf.Column(name='OBJ_ID', format='26A', array=objid)
    c2 = pf.Column(name='FLUX_MEAN',
                   format='1D',
                   unit='Counts',
                   array=meanarray)
    c3 = pf.Column(name='BLEND_FRACTION', format='1D', array=(meanarray) * 0)
    c4 = pf.Column(name='NPTS', format='1J', array=npts)
    c5 = pf.Column(name='RA', format='1D', array=RA)
    c6 = pf.Column(name='DEC', format='1D', array=DEC)

    a1 = pf.Column(name='HICOUNT', format='1J', array=zeros)
    a2 = pf.Column(name='LOCOUNT', format='1J', array=zeros)
    a3 = pf.Column(name='TMID', format='1D', unit='JD', array=tmid)
    a4 = pf.Column(name='ALT', format='1D', array=ALT)
    a5 = pf.Column(name='AZ', format='1D', array=AZ)
    a6 = pf.Column(name='RAPOS', format='1D', array=centerra)
    a7 = pf.Column(name='DECPOS', format='1D', array=centerdec)
    a8 = pf.Column(name='T', format='1D', array=T)
    a9 = pf.Column(name='MEANBIAS', format='1D', array=meanbias)
    a10 = pf.Column(name='ADU_DEV', format='1D', array=ADU_DEV)
    a11 = pf.Column(name='SKYLEVEL', format='1D', array=skylevel)
    a12 = pf.Column(name='FWHM', format='1D', array=fwhm)
    a13 = pf.Column(name='SEEING', format='1D', array=seeing)
    a14 = pf.Column(name='CLOUDS', format='1D', array=CLOUDS)
    a15 = pf.Column(name='SHIFT', format='1D', array=SHIFT)
    a16 = pf.Column(name='EXPOSURE', format='1D', array=exposure)
    a17 = pf.Column(name='IMAGE_ID',format='1K',array=imid)
    a18 = pf.Column(name='AIRMASS', format='1D', array=airmass)

    a19 = pf.Column(name='PSF_a_1', format='1D', array=fwhma_1)
    a20 = pf.Column(name='PSF_b_1',format='1D',array=fwhmb_1)
    a21 = pf.Column(name='PSF_t_1', format='1D', array=fwhmt_1)

    a22 = pf.Column(name='PSF_a_2', format='1D', array=fwhma_2)
    a23 = pf.Column(name='PSF_b_2',format='1D',array=fwhmb_2)
    a24 = pf.Column(name='PSF_t_2', format='1D', array=fwhmt_2)

    a25 = pf.Column(name='PSF_a_3', format='1D', array=fwhma_3)
    a26 = pf.Column(name='PSF_b_3',format='1D',array=fwhmb_3)
    a27 = pf.Column(name='PSF_t_3', format='1D', array=fwhmt_3)

    a28 = pf.Column(name='PSF_a_4', format='1D', array=fwhma_4)
    a29 = pf.Column(name='PSF_b_4',format='1D',array=fwhmb_4)
    a30 = pf.Column(name='PSF_t_4', format='1D', array=fwhmt_4)

    a31 = pf.Column(name='PSF_a_5', format='1D', array=fwhma_5)
    a32 = pf.Column(name='PSF_b_5',format='1D',array=fwhmb_5)
    a33 = pf.Column(name='PSF_t_5', format='1D', array=fwhmt_5)

    a34 = pf.Column(name='PSF_a_6', format='1D', array=fwhma_6)
    a35 = pf.Column(name='PSF_b_6',format='1D',array=fwhmb_6)
    a36 = pf.Column(name='PSF_t_6', format='1D', array=fwhmt_6)

    a37 = pf.Column(name='PSF_a_7', format='1D', array=fwhma_7)
    a38 = pf.Column(name='PSF_b_7',format='1D',array=fwhmb_7)
    a39 = pf.Column(name='PSF_t_7', format='1D', array=fwhmt_7)

    a40 = pf.Column(name='PSF_a_8', format='1D', array=fwhma_8)
    a41 = pf.Column(name='PSF_b_8',format='1D',array=fwhmb_8)
    a42 = pf.Column(name='PSF_t_8', format='1D', array=fwhmt_8)

    a43 = pf.Column(name='PSF_a_9', format='1D', array=fwhma_9)
    a44 = pf.Column(name='PSF_b_9',format='1D',array=fwhmb_9)
    a45 = pf.Column(name='PSF_t_9', format='1D', array=fwhmt_9)

    hducatalogue=pf.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6])

    hduimagelist=pf.BinTableHDU.from_columns([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45])

    hduprime = pf.PrimaryHDU(np.array(flux).T)

    hduflux = pf.ImageHDU(fluxarray)
    hdufluxerr = pf.ImageHDU(flux_err_array)
    hduxpos = pf.ImageHDU(np.array(xpos).T)
    hduypos = pf.ImageHDU(np.array(ypos).T)
    hdutime = pf.ImageHDU(hjdarray)
    hdunullq = pf.ImageHDU((np.array(time).T)*0 + 1)
    hduskylev = pf.ImageHDU(np.array(Skylev).T)
    hduskyrms = pf.ImageHDU(np.array(Skyrms).T)

    hduflux_1 = pf.ImageHDU(flux_grid[:,0].T)
    hduflux_2 = pf.ImageHDU(flux_grid[:,1].T)
    hduflux_3 = pf.ImageHDU(flux_grid[:,2].T)
    hduflux_4 = pf.ImageHDU(flux_grid[:,3].T)
    hduflux_5 = pf.ImageHDU(flux_grid[:,4].T)

    hduerror_1 = pf.ImageHDU(error_grid[:,0].T)
    hduerror_2 = pf.ImageHDU(error_grid[:,1].T)
    hduerror_3 = pf.ImageHDU(error_grid[:,2].T)
    hduerror_4 = pf.ImageHDU(error_grid[:,3].T)
    hduerror_5 = pf.ImageHDU(error_grid[:,4].T)

    hdulist = pf.HDUList([hduprime] + [hducatalogue] + [hduimagelist] + [hdutime] + [hduflux] +
                                             [hdufluxerr] + [hdunullq] + [hduxpos] + [hduypos] + [hduskylev] +
                                             [hduskyrms] + [hduflux_1] + [hduflux_2] + [hduflux_3] + [hduflux_4] + [hduflux_5] + [hduerror_1] + [hduerror_2] + [hduerror_3] + [hduerror_4] + [hduerror_5])

    hdulist[0].name = 'Primary'
    hdulist[1].name = 'CATALOGUE'
    hdulist[2].name = 'IMAGELIST'
    hdulist[3].name = 'HJD'
    hdulist[4].name = 'FLUX'
    hdulist[5].name = 'FLUXERR'
    hdulist[6].name = 'QUALITY'
    hdulist[7].name = 'CCDX'
    hdulist[8].name = 'CCDY'
    hdulist[9].name = 'SKYBKG'
    hdulist[10].name = 'SKYRMS'
    hdulist[11].name = 'FLUX_1'
    hdulist[12].name = 'FLUX_2'
    hdulist[13].name = 'FLUX_3'
    hdulist[14].name = 'FLUX_4'
    hdulist[15].name = 'FLUX_5'
    hdulist[16].name = 'ERROR_1'
    hdulist[17].name = 'ERROR_2'
    hdulist[18].name = 'ERROR_3'
    hdulist[19].name = 'ERROR_4'
    hdulist[20].name = 'ERROR_5'

    hdulist.writeto(outfile_name, clobber=True)


def stitch(filelist, appsize, outname):

    #combine the sub output files into a single master file, preserving the data types

    hdulist = []

    for filen in filelist:
        hdulist += [pf.open(filen)]

    headername = 'FLUX'
    catalogue = hdulist[0][headername].data
    combine = []
    for i in range(0, len(hdulist)):
        combine += list(hdulist[i][headername].data.T)
        print shape(combine)

    combine = array(combine).T

    fluxmean = mean(combine, axis=0)
    hduflux = pf.ImageHDU(combine)
    hduprime = pf.PrimaryHDU()

    a = []
    headername ='IMAGELIST'
    for column in hdulist[0][headername].columns:
        combine =[]
        colname = column.name
        for i in range(0,len(hdulist)):
            combine += list(hdulist[i][headername].data[colname])
        a += [pf.Column(name=colname, format=column.format, array=combine)]

    c = []
    headername = 'CATALOGUE'
    for column in hdulist[0][headername].columns:
        colname = column.name
        combine = list(hdulist[i][headername].data[colname])
        c += [pf.Column(name=colname, format=column.format, array=combine)]

    c[1].array = fluxmean

    headername_list = ['HJD', 'FLUXERR', 'QUALITY', 'CCDX', 'CCDY', 'SKYBKG',
                       'SKYRMS', 'FLUX_1', 'FLUX_2', 'FLUX_3', 'FLUX_4',
                       'FLUX_5', 'ERROR_1', 'ERROR_2', 'ERROR_3', 'ERROR_4',
                       'ERROR_5']
    dicty = {}

    for headername in headername_list:
        catalogue = hdulist[0][headername].data
        combine = []
        for i in range(0, len(hdulist)):
            combine += list(hdulist[i][headername].data.T)
        combine = array(combine).T
        dicty[headername] = pf.ImageHDU(combine)

    hduimagelist = pf.BinTableHDU.from_columns(a)
    hducatalogue = pf.BinTableHDU.from_columns(c)

    c1 = pf.Column(name='FLUX_MEAN',
                   format='1D',
                   unit='Counts',
                   array=fluxmean)

    new_hdulist = pf.HDUList(
        [hduprime] + [hducatalogue] + [hduimagelist] + [dicty['HJD']] +
        [hduflux] + [dicty['FLUXERR']] + [dicty['QUALITY']] + [dicty['CCDX']] +
        [dicty['CCDY']] + [dicty['SKYBKG']] + [dicty['SKYRMS']] + [
            dicty['FLUX_1']
        ] + [dicty['FLUX_2']] + [dicty['FLUX_3']] + [dicty['FLUX_4']] +
        [dicty['FLUX_5']] + [dicty['ERROR_1']] + [dicty['ERROR_2']] +
        [dicty['ERROR_3']] + [dicty['ERROR_4']] + [dicty['ERROR_5']])

    new_hdulist[0].name = 'Primary'
    new_hdulist[1].name = 'CATALOGUE'
    new_hdulist[2].name = 'IMAGELIST'
    new_hdulist[3].name = 'HJD'
    new_hdulist[4].name = 'FLUX'
    new_hdulist[5].name = 'FLUXERR'
    new_hdulist[6].name = 'QUALITY'
    new_hdulist[7].name = 'CCDX'
    new_hdulist[8].name = 'CCDY'
    new_hdulist[9].name = 'SKYBKG'
    new_hdulist[10].name = 'SKYRMS'
    new_hdulist[11].name = 'FLUX_1'
    new_hdulist[12].name = 'FLUX_2'
    new_hdulist[13].name = 'FLUX_3'
    new_hdulist[14].name = 'FLUX_4'
    new_hdulist[15].name = 'FLUX_5'
    new_hdulist[16].name = 'ERROR_1'
    new_hdulist[17].name = 'ERROR_2'
    new_hdulist[18].name = 'ERROR_3'
    new_hdulist[19].name = 'ERROR_4'
    new_hdulist[20].name = 'ERROR_5'

    mult = [1.0, 1.0, 0.5, np.sqrt(2.0), 2.0, 2.0 * np.sqrt(2.0), 4.0, 0.5,
            np.sqrt(2.0), 2.0, 2.0 * np.sqrt(2.0), 4.0]
    ind = [4, 5, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

    for i in range(0, len(mult)):
        new_hdulist[ind[i]].header['AP_MULT'] = (mult[i],
                                                 'Multiple of core radius')
        new_hdulist[ind[i]].header['AP_SIZE'] = (mult[i] * appsize,
                                                 'Aperture radius [Pixels]')

    new_hdulist.writeto(outname, clobber=True)


def compute_final_statistics(fname):
    '''
  Given a filename, for each aperture in the file compute the number of points in the lightcurve
  which are not NaNs

  This function is destructive - it changes the ouptut file
  '''
    with fitsio.FITS(fname, 'rw') as outfile:
        flux = outfile['flux'].read()
        fluxerr = outfile['fluxerr'].read()
        original_catalogue = outfile['catalogue']
        keys = original_catalogue.get_colnames()
        original_catalogue_data = original_catalogue.read()

        out_npts, out_flux_mean = [], []
        for (lc, lcerr, cat_row) in zip(flux, fluxerr, original_catalogue_data):
            ind = np.isfinite(lc)
            npts = lc[ind].size

            if npts > 0:
                flux_mean = np.average(lc[ind], weights=1. / lcerr[ind] ** 2)
            else:
                flux_mean = np.nan

            out_npts.append(npts)
            out_flux_mean.append(flux_mean)

        original_catalogue_data['NPTS'] = np.array(out_npts)
        original_catalogue_data['FLUX_MEAN'] = np.array(out_flux_mean)
        original_catalogue.write(original_catalogue_data)
