# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import

from astropy.io import fits as pf
import os
import sys
import linecache
import threading
from os.path import isfile, join
import multiprocessing
from multiprocessing import Pool
from functools import partial

from NGTS_workpackage.wcs_status import wcs_succeeded
from NGTS_workpackage.hjd_correction import append_hjd_correction_column
from NGTS_workpackage.quality_checks import *
from NGTS_workpackage.super_sample import call_find_fwhm
from NGTS_workpackage import casutools


def m_wcs_photom(filelist, outlist, appsize, conf_file, cat_file,
                 nproc=1,
                 verbose=False):

    os.system('cp ' + filelist + ' ' + outlist)

    infiles = []
    with open(filelist) as infile:
        for line in infile:
            image = line.strip()

            status_checks = ['ok', 'ok']

            if all(status == 'ok' for status in status_checks):
                infiles.append(image)

    pool = Pool(nproc)

    fn = partial(handle_errors_in_wcs_photom,
                 cat_file=cat_file,
                 conf_file=conf_file,
                 appsize=appsize,
                 verbose=verbose)
    pool.map(fn, infiles)

    first_image = infiles[0] + '.phot'
    pf.setval(first_image, 'SKY_MOVE', 1, value=0)

    RA_shift, DEC_shift, tot_shift, RA, DEC = frame_shift(infiles[0],
                                                          infiles[0])
    pf.setval(first_image, 'RA_MOVE', 1,
              value=RA_shift,
              comment='RA shift from previous image [arcseconds]')
    pf.setval(first_image, 'DEC_MOVE', 1,
              value=DEC_shift,
              comment='Dec shift from previous image [arcseconds]')
    pf.setval(first_image, 'SKY_MOVE', 1,
              value=tot_shift,
              comment='Total movement on sky [arcseconds]')

    pf.setval(first_image, 'WCSF_RA', 1, value=RA, comment='RA center pix')
    pf.setval(first_image, 'WCSF_DEC', 1, value=DEC, comment='Dec center pix')

    indexes = arange(1, len(infiles))
    fn = partial(m_frame_shift, infiles)
    pool.map(fn, indexes)


def handle_errors_in_wcs_photom(image, *args, **kwargs):
    try:
        return wcs_photom(image, *args, **kwargs)
    except Exception as err:
        print("Exception handled in wcs_photom: {}".format(str(err)),
                file=sys.stderr)

def wcs_photom(image,cat_file='nocat',conf_file='noconf',appsize=2.0,verbose=False):

    if not wcs_succeeded(image):
        return 'failed'

    outname = image + '.phot'

    casutools.imcore_list(image, cat_file, outname, confidence_map=conf_file,rcore=appsize, noell=True,
                        verbose=verbose)

     #      do some quality checks

    factor = 5
    size = 11
    stars = 1000
    side = 3 #how many subdevisions of the chip along each axis. must be odd.


    fwhm_a, fwhm_b, t = call_find_fwhm(image,cat_file,factor,size,stars,tag=image,side=side)
         
    cloud_status = cloud_check(image)

    pixel_fwhm = pf.getval(outname,'SEEING',1)
    plate_scale = 5.0 
    seeing = round(plate_scale*pixel_fwhm*3600,2)

    pf.setval(outname,'CLOUD_S',1,value=round(cloud_status,2),comment='A measure of bulk structure in the image (S/N)')
    pf.setval(outname,'FWHM',1,value=pixel_fwhm,comment='[pixels] Average FWHM')
    pf.setval(outname,'SEEING',1,value=seeing,comment='[arcseconds] Average FWHM')
    positions = ['Top left.','Top middle.','Top right.','Middle left.','Center.','Middle right.','Bottom left.','Bottom middle.','Bottom right.']

    for val in range(1,10):
        pf.setval(outname,'PSF_a_'+str(val),1,value=fwhm_a['f_'+str(val)][0],comment='[pixels] psf long axis FWHM. '+positions[val-1])
        pf.setval(outname,'PSF_b_'+str(val),1,value=fwhm_b['f_'+str(val)][0],comment='[pixels] psf short axis FWHM. '+positions[val-1])
        pf.setval(outname,'PSF_t_'+str(val),1,value=t['f_'+str(val)][0],comment='[degrees] psf rotation angle. '+positions[val-1])

    # Compute the HJD values
    append_hjd_correction_column(outname)

    return 'ok'
