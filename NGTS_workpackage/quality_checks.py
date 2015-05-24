# -*- coding: utf-8 -*-
from astropy.io import fits as pf
import os
import linecache
import threading
import multiprocessing
from os import listdir
from os.path import isfile, join
from util import thread_alloc
from numpy import *


def cum_guassian_fit(p, x, data):
    f = cum_guassian_func(p, x)
    return data - f


def cum_guassian_func(p, x):
    f = p[0] * norm.cdf(x / p[1])
    return f


def compute_frame_signal_to_noise(image_name):
    with pf.open(image_name) as imagedata:
        gain = 1.
        SNimage = (mean(imagedata[0].data) /
                   std(imagedata[0].data)) * sqrt(gain)
    return SNimage


def m_frame_shift(imagelist, index):
    try:
        image1 = imagelist[index - 1]
        image2 = imagelist[index]
        RA_shift, DEC_shift, tot_shift, RA, DEC = frame_shift(image1, image2)
        pf.setval(image2 + '.phot', 'RA_MOVE', 1,
                  value=RA_shift,
                  comment='RA shift from previous image [arcseconds]')
        pf.setval(image2 + '.phot', 'DEC_MOVE', 1,
                  value=DEC_shift,
                  comment='Dec shift from previous image [arcseconds]')
        pf.setval(image2 + '.phot', 'SKY_MOVE', 1,
                  value=tot_shift,
                  comment='Total movement on sky [arcseconds]')

        pf.setval(image2 + '.phot', 'WCSF_RA', 1,
                  value=RA_shift,
                  comment='RA center pix')
        pf.setval(image2 + '.phot', 'WCSF_DEC', 1,
                  value=DEC_shift,
                  comment='Dec center pix')
    except Exception as err:
        print "Exception handled in m_frame_shift: {}".format(str(err))


def frame_shift(image1, image2):

    print image1

    with pf.open(image1) as photdata:
        RA_prev = photdata[0].header['WCSF_RA']
        DEC_prev = photdata[0].header['WCSF_DEC']

    with pf.open(image2) as photdata:
        RA = photdata[0].header['WCSF_RA']
        DEC = photdata[0].header['WCSF_DEC']

    RA_shift = 3600 * (RA_prev - RA)
    DEC_shift = 3600 * (DEC_prev - DEC)

    tot_shift = ((RA_shift * cos(DEC * pi / 180)) ** 2 + DEC_shift ** 2) ** 0.5

    print tot_shift

    return RA_shift, DEC_shift, tot_shift, RA, DEC
