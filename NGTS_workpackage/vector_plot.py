# -*- coding: utf-8 -*-

"""
 
Create a vector plot of the astrometric errors. The errors are greatly exagerated for clarity.
Returns the median astrometric error for diagnostic use.

Usage:
  vector_plot.py [options] (<CATALOG_NAME>) (<IMAGE_NAME>) (<PLOT_NAME>)

Options:
  -h --help  Show help text

"""
from docopt import docopt
import numpy as np
from astropy import wcs
from astropy.io import fits
import fitsio
from NGTS_workpackage.catmatch import load_wcs_from_keywords
from NGTS_workpackage.catmatch import calc_seps

def wcsf_QCheck(catalog_name,image_name,plot_name,cat,RA_lims,DEC_lims,my_X,my_Y,upscale_factor=500,plot=True):

  print 'about to plot'

  plot_dir = ''

  im_header = fitsio.read_header(image_name)
  pix_coords = [[my_X[i],my_Y[i]] for i in range(0,len(my_X))]
  world = load_wcs_from_keywords(im_header,pix_coords)

  xs,ys,RA_sep,DEC_sep,x_sep,y_sep,sep_list = calc_seps(catalog_name,cat,RA_lims,DEC_lims,world,my_X,my_Y,im_header)

  rms = np.median(sep_list)

  #from scipy.optimize import leastsq
  #true_cen_xs = iter_poly(xs,RA_sep)
  #true_cen_ys = iter_poly(ys,DEC_sep)
  #print true_cen_xs, true_cen_ys

  h = fitsio.read_header(image_name)
  CRPIX1 =  h['CRPIX1']
  CRPIX2 =  h['CRPIX2']

  cen_X = h['NAXIS1']/2.0
  cen_Y = h['NAXIS2']/2.0

  cen = [[cen_X,cen_Y]]

  cen_world = load_wcs_from_keywords(im_header,cen)

  print 'cen is',cen_world

  #  axis.plot(true_cen_xs,true_cen_ys,'go',markersize=10)

  if plot == True:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axis = plt.subplots(figsize=(11, 8))
    for i in range(0,len(xs)):
      axis.plot([xs[i],(xs[i]+x_sep[i]*upscale_factor)],[ys[i],(ys[i]+y_sep[i]*upscale_factor)],'k-')
      axis.plot((xs[i]),(ys[i]),'ko',markersize=3)



    axis.plot(CRPIX1,CRPIX2,'ro',markersize=10)

    axis.set_ylim([0,2048])
    axis.set_xlim([0,2048])
    axis.set_aspect('equal')
    axis.set_xlabel(r'X')
    axis.set_ylabel(r'Y')
    axis.set_title(r'RMS: '+str(rms)+'"')
    fig.tight_layout()
    fig.savefig(plot_dir+plot_name, bbox_inches='tight')

    plt.close()

  print 'done plotting'

  with fitsio.FITS(image_name,'rw') as fits:
    fits[0].write_key('wcsf_ns',len(sep_list))
    fits[0].write_key('wcsf_RMS',rms)
    fits[0].write_key('wcsf_RA',cen_world[0][0])
    fits[0].write_key('wcsf_DEC',cen_world[0][1])

def iter_poly(xs,RA_sep):

  RA_sep = np.array(RA_sep)
  xs = np.array(xs)
  prior = [0.001,0.00000001,-0.000001,0.000001,500]
  xs_copy = xs.copy()
  RA_sep_copy = RA_sep.copy()

  for i in range(0,50):
    x, success = leastsq(fit_polynomial,prior,args=(xs_copy,RA_sep_copy),epsfcn=0.01)
    poly = polynomial(x,xs_copy)
    residuals = poly - RA_sep_copy
    std_r = np.std(residuals)
    xs_copy = xs_copy[abs(residuals) < 3.0*std_r]
    RA_sep_copy = RA_sep_copy[abs(residuals) < 3.0*std_r]

  #poly = polynomial(x,xs)

  #plt.plot(xs,poly,'b.')

  #plt.plot(xs,RA_sep,'r.')
  #plt.show()

  return x[-1]

def fit_polynomial(x,y,data):

  model = polynomial(x,y)

  resid = data - model

  return resid/0.01

def polynomial(x,y):

  model = y*0
  for i in range(0,len(x)-1.0):
    model += x[i]*(y - x[-1])**((i*2) + 1)
  return model

if __name__ == '__main__':

  args = docopt(__doc__)
  med_sep = plot_differences(args['<CATALOG_NAME>'],args['<IMAGE_NAME>'],args['<PLOT_NAME>'])

  print 'the median sep is',med_sep,' arcseconds'
