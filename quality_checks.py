# -*- coding: utf-8 -*-
import pyfits as pf
import os
import linecache
import threading
import multiprocessing
import scipy.optimize as opt
from os import listdir
from os.path import isfile, join
from util import thread_alloc
from numpy import *

def get_fwhm(photdata,appsize):

  from pylab import *

  # quick estimate of the median FWHM of the IQR range of the image by looking at the curve of growth.

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

  print median(core), median(core5)
  fit = cum_guassian_func(x,rad)
  plot(rad,median_cum_profile)
  plot(rad,fit)
  show()
  print x
  print fwhm

  return fwhm

def cum_guassian_fit(p,x,data):
  f = cum_guassian_func(p,x)
  return data - f

def cum_guassian_func(p,x):
  from scipy.stats import norm
  f = p[0]*norm.cdf(x/p[1])
  return f

def cloud_check(image_name):
  with pf.open(image_name) as imagedata:
    gain = imagedata[0].header['GAINFACT']
    SNimage = (mean(imagedata[0].data)/std(imagedata[0].data))*sqrt(gain)
  return SNimage
