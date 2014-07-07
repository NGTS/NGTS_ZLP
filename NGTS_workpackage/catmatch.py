# -*- coding: utf-8 -*-
from astropy import wcs
import numpy as np
import fitsio

def shift_wcs_axis(dicty,mycat,cat,RA_lims,DEC_lims,my_X,my_Y,TEL_RA,TEL_DEC,iters=1):

  pix_coords = [[my_X[i],my_Y[i]] for i in range(0,len(my_X))]

  print dicty['DEC_s'], dicty['RA_s']
  for i in range(0,iters):
    dicty['CRVAL1'] = TEL_RA + dicty['RA_s']
    dicty['CRVAL2'] = TEL_DEC + dicty['DEC_s']
    world = load_wcs_from_keywords(dicty,pix_coords)

    xs,ys,RA_sep,DEC_sep,x_sep,y_sep,sep_list = calc_seps(mycat,cat,RA_lims,DEC_lims,world,my_X,my_Y,dicty)

    med_RA = np.median(RA_sep[sep_list > np.median(sep_list)])
    med_DEC = np.median(DEC_sep[sep_list > np.median(sep_list)])

    dicty['RA_s'] += med_RA
    dicty['DEC_s'] += med_DEC
    print dicty['DEC_s'], dicty['RA_s'], np.median(sep_list)

  print 'done!'

  return dicty

def lmq_fit(best_fit,mycat,cat,RA_lims,DEC_lims,my_X,my_Y,TEL_RA,TEL_DEC,fitlist=['RA_s','DEC_s','CD1_1','CD1_2','CD2_1','CD2_2']):
  import scipy.optimize as opt

  name_list = []
  priorl = []
  for key in best_fit:
    if any( key ==  np.array(fitlist)):
      priorl += [best_fit[key]]
      name_list += [key]

  pix_coords = [[my_X[i],my_Y[i]] for i in range(0,len(my_X))]

#  x, success = opt.leastsq(lmq_fit_model,priorl,args=(mycat,cat,RA_lims,DEC_lims,my_X,my_Y,pix_coords,TEL_RA,TEL_DEC,name_list,best_fit),factor=10.0,epsfcn=0.0000001)

  x, success = opt.leastsq(lmq_fit_model,priorl,args=(mycat,cat,RA_lims,DEC_lims,my_X,my_Y,pix_coords,TEL_RA,TEL_DEC,name_list,best_fit),factor=10.0,epsfcn=0.00001)



  for i in range(0,len(name_list)):
    best_fit[name_list[i]] = x[i]

  return best_fit

def lmq_fit_model(vals,mycat,cat,RA_lims,DEC_lims,my_X,my_Y,pix_coords,TEL_RA,TEL_DEC,name_list,dicty):

  for i in range(0,len(vals)):
    dicty[name_list[i]] = vals[i]

  if (abs(dicty['RA_s']) > 1.0):
    return [2e6]*len(vals)

  if (abs(dicty['DEC_s']) > 1.0):
    return [2e6]*len(vals)


  dicty['CRVAL1'] = TEL_RA + dicty['RA_s']
  dicty['CRVAL2'] = TEL_DEC + dicty['DEC_s']

  world = load_wcs_from_keywords(dicty,pix_coords)

  print dicty

  xs,ys,RA_sep,DEC_sep,x_sep,y_sep,sep_list = calc_seps(mycat,cat,RA_lims,DEC_lims,world,my_X,my_Y,dicty)

  goodness = [np.median(sep_list)*2000]*len(vals)

  print vals
  print np.median(sep_list)

  return goodness
  
def fit_shift_wcs_axis(dicty,casuin,mycat,cat,XVAL,YVAL,TEL_RA,TEL_DEC,RA_lims,DEC_lims,my_X,my_Y,pix_coords,thresh=100,reset=False,update=False):

  world = load_wcs_from_keywords(dicty,pix_coords)

  if update == True:
    apply_correct(dicty,casuin,TEL_RA,TEL_DEC)

  xs,ys,RA_sep,DEC_sep,x_sep,y_sep,sep_list = calc_seps(mycat,cat,RA_lims,DEC_lims,world,my_X,my_Y,dicty)

  return np.array(sep_list)

def apply_correct(dicty,casuin,TEL_RA,TEL_DEC):

  dicty['CRVAL1'] = TEL_RA + dicty['RA_s']
  dicty['CRVAL2'] = TEL_DEC + dicty['DEC_s']

  with fitsio.FITS(casuin,'rw') as fits:
    for key in dicty:
      fits[0].write_key(key,dicty[key])
  
def calc_seps(mycat,cat,RA_lims,DEC_lims,world,my_X,my_Y,dicty,in_test=2000):

  cat_RA_raw = cat['ra']
  cat_DEC_raw = cat['dec']

  zero = 21.5

  cat_Jmag = cat['Jmag']

  my_mag = zero - 2.512*np.log10(mycat['Aper_flux_3'])
    
  my_RA_raw = world[:,0]
  my_DEC_raw = world[:,1]

  sep_list = []
  DEC_sep = []
  RA_sep = []
  xs = []
  ys = []
  x_sep = []
  y_sep = []

  try:
    my_RA = my_RA_raw[(my_RA_raw > RA_lims[0][0]) & (my_RA_raw < RA_lims[0][1]) & (my_DEC_raw > DEC_lims[0][0]) & (my_DEC_raw < DEC_lims[0][1])]
    my_DEC = my_DEC_raw[(my_RA_raw > RA_lims[0][0]) & (my_RA_raw < RA_lims[0][1]) & (my_DEC_raw > DEC_lims[0][0]) & (my_DEC_raw < DEC_lims[0][1])]
    cat_RA = cat_RA_raw[(cat_RA_raw > min(my_RA)) & (cat_RA_raw < max(my_RA)) & (cat_DEC_raw > min(my_DEC)) & (cat_DEC_raw < max(my_DEC))]
    cat_DEC = cat_DEC_raw[(cat_RA_raw > min(my_RA)) & (cat_RA_raw < max(my_RA)) & (cat_DEC_raw > min(my_DEC)) & (cat_DEC_raw < max(my_DEC))]
  except:
    return xs,ys,RA_sep,DEC_sep,x_sep,y_sep,np.array([1000.0])

  my_X = my_X[(my_RA_raw > RA_lims[0][0]) & (my_RA_raw < RA_lims[0][1]) & (my_DEC_raw > DEC_lims[0][0]) & (my_DEC_raw < DEC_lims[0][1])]
  my_Y = my_Y[(my_RA_raw > RA_lims[0][0]) & (my_RA_raw < RA_lims[0][1]) & (my_DEC_raw > DEC_lims[0][0]) & (my_DEC_raw < DEC_lims[0][1])]

  my_brightest = np.argsort(my_mag[(my_RA_raw > RA_lims[0][0]) & (my_RA_raw < RA_lims[0][1]) & (my_DEC_raw > DEC_lims[0][0]) & (my_DEC_raw < DEC_lims[0][1])])[:in_test]
  c_b = np.argsort(cat_Jmag[(cat_RA_raw > min(my_RA)) & (cat_RA_raw < max(my_RA)) & (cat_DEC_raw > min(my_DEC)) & (cat_DEC_raw < max(my_DEC))])[:in_test]

  w = wcs.WCS(dicty)

  for i in my_brightest:
    RA = my_RA[i]
    DEC = my_DEC[i]
    sep = 3600*(((RA - cat_RA[c_b])*(np.cos(DEC*np.pi/180.0)))**2.0 + (DEC - cat_DEC[c_b])**2.0)**0.5
    index = np.argmin(sep)
    sep_list +=[sep[index]]
    RA_sep += [cat_RA[c_b][index] - RA]
    DEC_sep += [cat_DEC[c_b][index] - DEC]
    xs += [my_X[i]]
    ys += [my_Y[i]]
    cat_pix = w.wcs_world2pix([[cat_RA[c_b][index],cat_DEC[c_b][index]]], 1)
    x_sep += [cat_pix[0][0] - my_X[i]]
    y_sep += [cat_pix[0][1] - my_Y[i]]

  course_seps = np.array(sep_list)
  course_fit = np.median(course_seps)

  xs = np.array(xs)
  ys = np.array(ys)
  RA_sep = np.array(RA_sep)
  DEC_sep = np.array(DEC_sep)
  x_sep = np.array(x_sep)
  y_sep = np.array(y_sep)
  course_seps = np.array(course_seps)

  old_n = len(course_seps)
  n = old_n - 1
  while ((old_n-n) > 0):
    old_n = len(course_seps)
    stdev = np.std(course_seps)
    course_fit = np.median(course_seps)
#      print 'std',stdev,'med',course_fit,'n',len(course_seps)
    fs = [course_seps < (course_fit + 5.0*stdev)]
    xs = xs[fs]
    ys = ys[fs]
    RA_sep = RA_sep[fs]
    DEC_sep = DEC_sep[fs]
    x_sep = x_sep[fs]
    y_sep = y_sep[fs]
    course_seps = course_seps[fs]
    n = len(course_seps)
  return xs,ys,RA_sep,DEC_sep,x_sep,y_sep,course_seps
  
def old_fit_shift():

  for i in range(0,len(my_RA)):
    RA = my_RA[i]
    DEC = my_DEC[i]
    sep = 3600*(((RA - cat_RA)*(np.cos(DEC*np.pi/180.0)))**2.0 + (DEC - cat_DEC)**2.0)**0.5
    index = np.argmin(sep)
    sep_list +=[sep[index]]
    RA_sep += [cat_RA[index] - RA]
    DEC_sep += [cat_DEC[index] - DEC]
    xs += [my_X[i]]
    ys += [my_Y[i]]

  sep_list = np.array(sep_list)
  RA_sep = np.array(RA_sep)
  DEC_sep = np.array(DEC_sep)
  xs = np.array(xs)
  ys = np.array(ys)

  c = [sep_list < 3*5.0]

  print 'got right to the end!'

  return xs[c],ys[c],RA_sep[c],DEC_sep[c], sep_list[c]

def correct_catfile(catfile_name,image_name,nstars=False):

  with fitsio.FITS(catfile_name,'rw') as cat:
    im_header = fitsio.read_header(image_name)
    my_X = cat[1]['X_coordinate'][:]
    old_RA = cat[1]['RA'][:]
    my_Y = cat[1]['Y_coordinate'][:]
    pix_coords = [[my_X[i],my_Y[i]] for i in range(0,len(my_X))]
    world = load_wcs_from_keywords(im_header,pix_coords)
    RA = world[:,0]*np.pi/180
    DEC = world[:,1]*np.pi/180
    cat[1].write_column('RA',RA,clobber=True)
    cat[1].write_column('DEC',DEC,clobber=True)
  

def load_wcs_from_keywords(fheader,pixcrd):
# Load the WCS information from a fits header, and use it
# to convert pixel coordinates to world coordinates.
   # Load the FITS hdulist using astropy.io.fits

    #hdulist = fits.open(filename)
    #fheader = hdulist[0].header

    # Parse the WCS keywords in the primary HDU
    w = wcs.WCS(fheader)

    # Convert pixel coordinates to world coordinates
    # The second argument is "origin" -- in this case we're declaring we
    # have 1-based (Fortran-like) coordinates.
    world = w.wcs_pix2world(pixcrd, 1)

    return world
