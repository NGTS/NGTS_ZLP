# -*- coding: utf-8 -*-
from NGTS_workpackage.catmatch import *
from emcee_tools import run_emcee
import numpy as np
import fitsio
import astropy.io.fits as pf
import sys


def main():

  casuin = sys.argv[1]
  mycatname = sys.argv[2]
  chain_name = sys.argv[3]
  catsrc = sys.argv[4]

  resume = False
  nwalkers = 1000
  nruns = 1e7
  nthreads = 16
  burns = 0
  start_size = 1e-2
  
  #if sys.argv[-1] == 'jun19':
    #casuin = 'june_19th.fits'
    #mycatname = 'new_19th_jun_cat.fits'
    #chain_name = '10_jun19.fits'
    #catsrc = '/storage/astro2/phrmat/NGTS_june5_astrometry/new_19th_catcache'

  #if sys.argv[-1] == 'jun11':
    #casuin = 'jun11.fits'
    #mycatname = 'jun11_cat.fits'
    #chain_name = '2_jun11.fits'
    #catsrc = '/storage/astro2/phrmat/NGTS_june5_astrometry/catcache_11th'

  #if sys.argv[-1] == 'jun26':
    #casuin = 'jun26.fits'
    #mycatname = 'jun26_cat.fits'
    #chain_name = '2_jun26.fits'
    #catsrc = '/storage/astro2/phrmat/NGTS_june5_astrometry/catcache'

  hdulist = fitsio.read_header(casuin)
  XVAL = hdulist['NAXIS1']/2
  YVAL = hdulist['NAXIS2']/2
  TEL_RA = hdulist['TEL_RA']
  TEL_DEC = hdulist['TEL_DEC']

  cat_names = []
  RA_lims = []
  DEC_lims = []
  for line in open(catsrc+'/index'):
    vals = line.strip('\n').split(' ')
    cat_names += [vals[0]]
    RA_lims += [[float(vals[2]),float(vals[3])]]
    DEC_lims += [[float(vals[4]),float(vals[5])]]

  cen_RA = np.array([(f[0] + f[1])/2.0 for f in RA_lims])
  cen_DEC = np.array([(f[0] + f[1])/2.0 for f in DEC_lims])

  sep = (((TEL_RA - cen_RA)*(np.cos(TEL_DEC*np.pi/180.0)))**2.0 + (TEL_DEC - cen_DEC)**2.0)**0.5

  cat_name = cat_names[np.argmin(sep)]

  with pf.open(catsrc+'/'+cat_name) as catd:
    catt = catd[1].data.copy()
  cat = {'ra':catt['ra'],'dec':catt['dec'],'Jmag':catt['Jmag']}

  with fitsio.FITS(mycatname) as mycatt:
    mycat = {'Aper_flux_3':mycatt[1]['Aper_flux_3'][:]}  
    my_X = mycatt[1]['x_coordinate'][:]
    my_Y = mycatt[1]['y_coordinate'][:]
  pix_coords = [[my_X[i],my_Y[i]] for i in range(0,len(my_X))]

  # 7th order fit
  dicty = {'CRPIX1': 1.03259815e+03,'CRPIX2': 9.65505144e+02,'CD1_1': 1.41142333e-03,'CD2_2': 1.41109400e-03,'CD1_2': -1.89116218e-06,'CD2_1': 1.53342393e-06 ,'PV2_1': 1.0,'PV2_3': 8.68515702e+00,'PV2_5': 2.70336203e+02,'PV2_7': 1.37726138e+04,'RA_s':-0.45807896397,'DEC_s':0.48575139999,'CTYPE1':'RA---ZPN','CTYPE2':'DEC--ZPN'}
  name_list = ['CRPIX1','CRPIX2','CD1_1','CD2_2','CD1_2','CD2_1','PV2_3','PV2_5','PV2_7','RA_s','DEC_s']

  # 5th order fit
#  dicty = {'CRPIX1': 1.03186582e+03,'CRPIX2': 9.65390145e+02,'CD1_1': 1.41143441e-03,'CD2_2': 1.41109881e-03,'CD1_2': -1.91498220e-06,'CD2_1': 1.57312392e-06 ,'PV2_1': 1.0,'PV2_3': 8.71721335e+00,'PV2_5': 2.18288567e+02,'PV2_7': 0,'RA_s':-0.45966596397,'DEC_s':0.48558489999,'CTYPE1':'RA---ZPN','CTYPE2':'DEC--ZPN'}
#  name_list = ['CRPIX1','CRPIX2','CD1_1','CD2_2','CD1_2','CD2_1','PV2_3','PV2_5','RA_s','DEC_s']

# june 11th fit
#  dicty = {'CD2_1': 1.5415410350820255e-06, 'CD2_2': 0.0014112047383328613, 'RA_s': -0.52982456599841854, 'CD1_2': -1.6950507735128331e-06, 'CD1_1': 0.0014112627531226202, 'CRVAL2': 49.59559537217558, 'CRPIX1': 1005.5871367339122, 'CRPIX2': 963.33376370419967, 'CRVAL1': 285.36598106573831, 'PV2_1': 1.0, 'PV2_2': 0.0, 'PV2_3': 9.3291677327369538, 'PV2_5': -542.20833721702797, 'PV2_7': 335143.48864371137, 'DEC_s': 0.40305908468262103, 'CTYPE2': 'DEC--ZPN', 'CTYPE1': 'RA---ZPN'}
#  name_list = ['CRPIX1','CRPIX2','CD1_1','CD2_2','CD1_2','CD2_1','PV2_3','PV2_5','PV2_7','RA_s','DEC_s']

# best fit for jun19
#  dicty = {'CD2_1': -9.3263106728144037e-07, 'CD2_2': 0.0013885458770754415, 'RA_s': -0.042574725586388078, 'CD1_2': 7.6821287599453613e-07, 'CD1_1': 0.0013884790887784609, 'CRVAL2': 49.099782962962607, 'CRPIX1': 1001.8384609589782, 'CRPIX2': 976.12370523564334, 'CRVAL1': 285.85331052744289, 'PV2_1': 1.0, 'PV2_2':0,'PV2_3': 8.3820912409837796, 'PV2_5': 690.45326416643377, 'PV2_7': -12173.91504945435, 'DEC_s': -0.092710638655683769, 'CTYPE2': 'DEC--ZPN', 'CTYPE1': 'RA---ZPN'}
#  name_list = ['CRPIX1','CRPIX2','CD1_1','CD2_2','CD1_2','CD2_1','PV2_3','PV2_5','PV2_7','RA_s','DEC_s']

# best fit for jun19, with the extra 2_2 parameter
#  dicty = {'CD2_1': -9.2503034157791214e-07, 'CD2_2': 0.0013884136198712995, 'RA_s': -0.042398901584286913, 'CD1_2': 7.622982784978669e-07, 'CD1_1': 0.001388342294464743, 'CRVAL2': 49.102018639433325, 'CRPIX1': 1001.9200387729178, 'CRPIX2': 977.74405281790348, 'CRVAL1': 285.85348635144499, 'PV2_1': 1.0, 'PV2_2': -0.0026221070341301313, 'PV2_3': 7.9629546187093201, 'PV2_5': 1523.5423889554663, 'PV2_7': -408922.22743904917, 'DEC_s': -0.090474962184963773, 'CTYPE2': 'DEC--ZPN', 'CTYPE1': 'RA---ZPN'}
#  name_list = ['CRPIX1','CRPIX2','CD1_1','CD2_2','CD1_2','CD2_1','PV2_3','PV2_5','PV2_7','RA_s','DEC_s']

#best fit for jun26

#  dicty = {'CD2_1': -3.1791182393163622e-09, 'CD2_2': 0.0013885328356863637, 'RA_s': 0.055516721576799946, 'CD1_2': -1.9130096940647768e-07, 'CD1_1': 0.0013883892608630331, 'CRVAL2': 49.413621212216562, 'CRPIX1': 1021.8925006108819, 'CRPIX2': 981.5660300722983, 'CRVAL1': 285.95141212097275, 'PV2_1': 1.0, 'PV2_2': 0.0, 'PV2_3': 8.0048907594320582, 'PV2_5': 1177.9976128674671, 'PV2_7': -158058.39490829041, 'DEC_s': 0.22113353324910145, 'CTYPE2': 'DEC--ZPN', 'CTYPE1': 'RA---ZPN'}
#  name_list = ['CRPIX1','CRPIX2','CD1_1','CD2_2','CD1_2','CD2_1','PV2_3','PV2_5','PV2_7','RA_s','DEC_s']

  dicty['CRVAL1'] = TEL_RA + dicty['RA_s']
  dicty['CRVAL2'] = TEL_DEC + dicty['DEC_s']

  cen = [[dicty['CRPIX1'],dicty['CRPIX2']]]

  old_world = load_wcs_from_keywords(hdulist,cen)

  TEL_RA = hdulist['TEL_RA']
  TEL_DEC = hdulist['TEL_DEC']

  dicty['RA_s'] = (old_world[0][0] - TEL_RA)
  dicty['DEC_s'] = (old_world[0][1] - TEL_DEC) 

  apply_correct(dicty,casuin,TEL_RA,TEL_DEC) 

  prior = []

#  for i in name_list:
#    prior += [dicty[i]]

  for i in name_list:
    prior += [hdulist[i]]
  prior[-1] = hdulist['CRVAL2'] - TEL_DEC
  prior[-2] = hdulist['CRVAL1'] - TEL_RA

  start_size = np.array([1e-3]*len(prior))

  for i in [2,3,4,5]:
    start_size[i] = 1e-3

  for i in [6,7,8,9]:
    start_size[i] = 1e-3

#  for i in range(0,len(prior)):
#    dicty[name_list[i]] = prior[i]

#  dicty['CRVAL1'] = TEL_RA + dicty['RA_s']
#  dicty['CRVAL2'] = TEL_DEC + dicty['DEC_s']

  args = [casuin,mycat,cat,XVAL,YVAL,TEL_RA,TEL_DEC,RA_lims,DEC_lims,my_X,my_Y,pix_coords,name_list,dicty]

#  x, success = leastsq(lmq_lnprob,prior,args=(casuin,mycat,cat,XVAL,YVAL,TEL_RA,TEL_DEC,RA_lims,DEC_lims,my_X,my_Y,pix_coords,dicty),epsfcn=1.0)
#  print x
#  quit()

#  run_emcee(prior,lnprob_initial,args,nwalkers,nruns,start_size,chain_name,burns,nthreads=nthreads,w=True,resume=resume)

  #rms = fit_shift_wcs_axis(prior,casuin,mycat,cat,XVAL,YVAL,TEL_RA,TEL_DEC,RA_lims,DEC_lims,my_X,my_Y,pix_coords)
  #os.system('python ../QCplots/vector_plot.py test_cat.fits test_saved.fits wtf.png')
  #quit()

#  dicty = {'CD2_1': -3.1791182393163622e-09, 'CD2_2': 0.0013885328356863637, 'RA_s': 0.055516721576799946, 'CD1_2': -1.9130096940647768e-07, 'CD1_1': 0.0013883892608630331, 'CRVAL2': 49.413621212216562, 'CRPIX1': 1021.8925006108819, 'CRPIX2': 981.5660300722983, 'CRVAL1': 285.95141212097275, 'PV2_1': 1.0, 'PV2_2': 0.0, 'PV2_3': 8.0048907594320582, 'PV2_5': 1177.9976128674671, 'PV2_7': -158058.39490829041, 'DEC_s': 0.22113353324910145, 'CTYPE2': 'DEC--ZPN', 'CTYPE1': 'RA---ZPN'}

#  name_list = ['CRPIX1','CRPIX2','CD1_1','CD2_2','CD1_2','CD2_1']

#  name_list = ['CRVAL1']

#  dicty['CRVAL1'] += 0.5

#  best_fit = lmq_fit(dicty.copy(),mycat,cat,RA_lims,DEC_lims,my_X,my_Y,TEL_RA,TEL_DEC,fitlist=name_list,factor=10.0,epsfcn=0.000001,shift=True)

#  print best_fit['CRVAL1'] - (dicty['CRVAL1'] - 0.5)

#  dicty = {'CD1_2': -2.0456617290539012e-07, 'CD1_1': 0.001384264414284212, 'CRVAL2': 49.435782429656442, 'CRPIX1': 1033.6952168544094, 'CRPIX2': 997.65007571830984, 'CRVAL1': 285.97633819019245, 'PV2_1': 1.0, 'PV2_2': 0.0, 'PV2_3': 8, 'PV2_5': 500, 'PV2_7': 0, 'DEC_s': 0.24329475068897921, 'CD2_1': 1.6349676741442095e-08, 'CTYPE2': 'DEC--ZPN', 'CTYPE1': 'RA---ZPN', 'CD2_2': 0.0013844216429580033, 'RA_s': 0.080442790796550118}

#  name_list = ['PV2_3','PV2_5','PV2_7']
#  best_fit = lmq_fit(dicty,mycat,cat,RA_lims,DEC_lims,my_X,my_Y,TEL_RA,TEL_DEC,fitlist=name_list,factor=1.0,epsfcn=0.01,shift=False)

#  name_list = ['CRPIX1','CRPIX2','CD1_1','CD2_2','CD1_2','CD2_1','PV2_3','PV2_5','PV2_7']
#  best_fit = lmq_fit(dicty,mycat,cat,RA_lims,DEC_lims,my_X,my_Y,TEL_RA,TEL_DEC,fitlist=name_list,factor=100.0,epsfcn=0.001,shift=False)

  run_emcee(prior,lnprob,args,nwalkers,nruns,start_size,chain_name,burns,nthreads=nthreads,w=True,resume=resume)

def lnprob(x,casuin,mycat,cat,XVAL,YVAL,TEL_RA,TEL_DEC,RA_lims,DEC_lims,my_X,my_Y,pix_coords,name_list,dicty):

  for i in range(0,len(x)):
    dicty[name_list[i]] = x[i]

  dicty['CRVAL1'] = TEL_RA + dicty['RA_s']
  dicty['CRVAL2'] = TEL_DEC + dicty['DEC_s']

  try:
    rms = fit_shift_wcs_axis(dicty,casuin,mycat,cat,XVAL,YVAL,TEL_RA,TEL_DEC,RA_lims,DEC_lims,my_X,my_Y,pix_coords)
  except:
    return -np.inf

  likelyhood = -2000*((np.median(rms))**2.0)

  lp = lnprior(dicty,rms)
  lnoutput = lp + likelyhood

  if not np.isfinite(lp):
      return -np.inf

  print lnoutput, np.median(rms)

  return lnoutput

def lnprior(dicty,rms):

  if np.all((len(rms) > 1000) and (0.0012 < dicty['CD1_1'] < 0.0025)  and (0.0012 < dicty['CD1_1'] < 0.0025) and ((abs(dicty['CD2_1']) < 1e-4)) and (abs(dicty['CD1_2']) < 1e-4) and (abs(dicty['RA_s']) < 1.0) and (abs(dicty['DEC_s'] < 1.0))):
      return 0.0
  return -np.inf

main()
