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

def get_fwhm(fname,appsize):


  # quick estimate of the median FWHM of the IQR range of the image by looking at the curve of growth.
  # defunct - replaced by reading 'seeing' value from CASU output

  photdata = pf.open(fname)


  mean_fluxes = photdata[1].data['Core3_flux']

  IQR = [(mean_fluxes < (median(mean_fluxes[mean_fluxes > median(mean_fluxes)]))) & (mean_fluxes > (median(mean_fluxes[mean_fluxes < median(mean_fluxes)])))]

  core = photdata[1].data['Core_flux'][IQR]
  core1 = photdata[1].data['Core1_flux'][IQR]
  core2 = photdata[1].data['Core2_flux'][IQR]
  core3 = photdata[1].data['Core3_flux'][IQR]
  core4 = photdata[1].data['Core4_flux'][IQR]
  core5 = photdata[1].data['Core5_flux'][IQR]
  total = core5
  annulus1 = (core1/total)
  annulus = (core/total)
  annulus2 = (core2/total) - (annulus)
  annulus3 = (core3/total) - (annulus + annulus2)
  annulus4 = (core4/total) - (annulus + annulus2 + annulus3)
  annulus5 = (core5/total) - (annulus + annulus2 + annulus3 + annulus4)

  profile = array([annulus1,annulus,annulus2,annulus3,annulus4])
  cum_profile = array([core1/total,core/total, core2/total, core3/total, core4/total, core5/total])

  median_profile = median(profile, axis = 1)
  median_cum_profile = median(cum_profile, axis = 1)
  rad = appsize*array([0.5,1.0,sqrt(2.0),2.0,2.0*sqrt(2.0),4])
  x1 = [0.25,2.0]

  rad = rad[:3]
  median_cum_profile = median_cum_profile[:3]

  x, success = opt.leastsq(cum_guassian_fit, x1, args=(rad,median_cum_profile))
  x[1] = abs(x[1])
  fwhm = (sqrt(2*log(2)))*x[1]

  return fwhm

def cum_guassian_fit(p,x,data):
  f = cum_guassian_func(p,x)
  return data - f

def cum_guassian_func(p,x):
  f = p[0]*norm.cdf(x/p[1])
  return f

def cloud_check(image_name):
  with pf.open(image_name) as imagedata:
    gain = imagedata[0].header['GAINFACT']
    SNimage = (mean(imagedata[0].data)/std(imagedata[0].data))*sqrt(gain)
  return SNimage

def m_frame_shift(imagelist,index):
  image1 = imagelist[index-1] + '.phot'
  image2 = imagelist[index] + '.phot'
  shift = frame_shift(image1,image2)
  plate_scale = pf.getval(image2,'CD1_1',1)
  pf.setval(image2,'SHIFT',1,value=shift,comment='[pixels] shift from previous image')
  pf.setval(image2,'SKYSHIFT',1,value=shift*plate_scale*3600,comment='[arcseconds] shift from previous image')

def frame_shift(image1,image2):

  with pf.open(image1) as photdata:
    xpos_prev = photdata[1].data['X_coordinate'].copy()
    ypos_prev = photdata[1].data['Y_coordinate'].copy()

  with pf.open(image2) as photdata:
    xpos = photdata[1].data['X_coordinate'].copy()
    ypos = photdata[1].data['Y_coordinate'].copy()

  xshift = sqrt(mean((xpos - xpos_prev)**2))
  yshift = sqrt(mean((ypos - ypos_prev)**2))
  shift = sqrt(xshift**2 + yshift**2)
  return shift
