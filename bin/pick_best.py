# -*- coding: utf-8 -*-
import sys
from fitsio import FITS
import numpy as np
from astropy.io import fits as pf
import pickle
import os

"""
Pickle best solution from MCMC chain, ready for pipeline usage.

Usage:
  Pick_best.py <chain_name> <image> <outname>

input:

  chain_name          Name of the mcmc chain (output of emcee_catmatch.py)
  image               image solved to make mcmc chain
  outname             Where to save the best solution (do it somewhere that pipeline sees)


"""

def pickle_best(chain_name,image,outname,maxlength=1e6):

  chaindata, length = load_chain(chain_name,maxlength=maxlength)

  print 'Picking the best solution from the last ',int(maxlength),' links of a ',length,' long chain.'

  lp = chaindata['lp']
  samples = chaindata['x'][np.isfinite(lp) & (lp != 0.0)]
  lp = lp[np.isfinite(lp) & (lp != 0.0)]

  opt = samples[np.argmax(lp)]

  best = max(lp)

  print 'rms:', round((best/-2000)**0.5,2), 'arcseconds'

  with pf.open(image) as hdulist:
    XVAL = hdulist[0].header['NAXIS1']/2
    YVAL = hdulist[0].header['NAXIS2']/2
    TEL_RA = hdulist[0].header['TEL_RA']
    TEL_DEC = hdulist[0].header['TEL_DEC']

  if len(opt)==12:
    # 7th order fit
    dist_map = {'CRPIX1': 1.03259815e+03,'CRPIX2': 9.65505144e+02,'CD1_1': 1.41142333e-03,'CD2_2': 1.41109400e-03,'CD1_2': -1.89116218e-06,'CD2_1': 1.53342393e-06 ,'PV2_1': 1.0,'PV2_2':2.0,'PV2_3': 8.68515702e+00,'PV2_5': 2.70336203e+02,'PV2_7': 1.37726138e+04,'RA_s':-0.45807896397,'DEC_s':0.48575139999,'CTYPE1':'RA---ZPN','CTYPE2':'DEC--ZPN'}
    param_names = ['CRPIX1','CRPIX2','CD1_1','CD2_2','CD1_2','CD2_1','PV2_2','PV2_3','PV2_5','PV2_7','RA_s','DEC_s']
  elif len(opt)==11:
    # 7th order fit
    dist_map = {'CRPIX1': 1.03259815e+03,'CRPIX2': 9.65505144e+02,'CD1_1': 1.41142333e-03,'CD2_2': 1.41109400e-03,'CD1_2': -1.89116218e-06,'CD2_1': 1.53342393e-06 ,'PV2_1': 1.0,'PV2_2':0.0,'PV2_3': 8.68515702e+00,'PV2_5': 2.70336203e+02,'PV2_7': 1.37726138e+04,'RA_s':-0.45807896397,'DEC_s':0.48575139999,'CTYPE1':'RA---ZPN','CTYPE2':'DEC--ZPN'}
    param_names = ['CRPIX1','CRPIX2','CD1_1','CD2_2','CD1_2','CD2_1','PV2_3','PV2_5','PV2_7','RA_s','DEC_s']
  else:
    # 5th order fit
    dist_map = {'CRPIX1': 1.03186582e+03,'CRPIX2': 9.65390145e+02,'CD1_1': 1.41143441e-03,'CD2_2': 1.41109881e-03,'CD1_2': -1.91498220e-06,'CD2_1': 1.57312392e-06 ,'PV2_1': 1.0,'PV2_2':0,'PV2_3': 8.71721335e+00,'PV2_5': 2.18288567e+02,'PV2_7': 0,'RA_s':-0.45966596397,'DEC_s':0.48558489999,'CTYPE1':'RA---ZPN','CTYPE2':'DEC--ZPN'}
    param_names = ['CRPIX1','CRPIX2','CD1_1','CD2_2','CD1_2','CD2_1','PV2_2','PV2_3','PV2_5','RA_s','DEC_s']

  for i in range(0,len(opt)):
    dist_map[param_names[i]] = opt[i]

  dist_map['CRVAL1'] = TEL_RA + dist_map['RA_s']
  dist_map['CRVAL2'] = TEL_DEC + dist_map['DEC_s']

  print dist_map

  pickle.dump(dist_map,open(outname,'w'))

def load_chain(chain_name,maxlength=1e6):

  #keeping this seperate for now in case I want to implement burning or something.

  with FITS(chain_name,'r') as fits:
    length = fits['MCMC']._info['nrows']
    if length > maxlength:
      data_dict = fits['MCMC'][length-maxlength:length]
    else:
      data_dict = fits['MCMC'][0:length]
      
  return data_dict, length

if __name__ == "__main__":

  if len(sys.argv) != 4:
    print 'Usage: Pick_best.py <chain_name> <image> <outname>'
    quit()

  chain_name = os.path.abspath(sys.argv[1])

  image = os.path.abspath(sys.argv[2])

  outname = os.path.abspath(sys.argv[3])

  maxlength = 1e5

  pickle_best(chain_name,image,outname,maxlength=1e6)
