# -*- coding: utf-8 -*-
from astropy.io import fits as pf
import os
from pylab import *
import scipy.optimize as opt
from util import load_wcs_from_file
import os

def shift_wcs_axis(casuin,mycatname,thresh=100):

  #this is the best solution for the catalog
  CRPIX1  =  -30.3589661445
  CRPIX2  =  -1.76800252831
  CRVAL1  =  -0.0134278702942
  CRVAL2  =  -0.0444093066511
  CD1_1   =  0.00138877871593
  CD2_2   =  0.00138875543313
  CD1_2   =  1.41054781057e-05
  CD2_1   =  -1.41656423353e-05
  PV2_1   =  0.999993897433
  PV2_3   =  8.11292725428
  PV2_5   =  901.974288037

  with pf.open(casuin) as hdulist:
    XVAL = hdulist[0].header['NAXIS1']/2
    YVAL = hdulist[0].header['NAXIS2']/2
    TEL_RA = hdulist[0].header['TEL_RA']
    TEL_DEC = hdulist[0].header['TEL_DEC']

  prior = [CRPIX1,CRPIX2,CRVAL1,CRVAL2,CD1_1,CD2_2,CD1_2,CD2_1,PV2_1,PV2_3,PV2_5]

  apply_correct(prior,casuin,XVAL,YVAL,TEL_RA,TEL_DEC)

  # this corrects efficiently for the 'offset' between the optical axis and the RA/DEC optical axis
  xs,ys,RA_sep,DEC_sep, sep_list = calc_seps(mycatname,casuin)
  prior[2] += median(RA_sep)
  prior[3] += median(DEC_sep)
  apply_correct(prior,casuin,XVAL,YVAL,TEL_RA,TEL_DEC)

def testbed():

  conf_dir = '../jul12/'



  #this is the best solution for the catalog
  CRPIX1  =  -30.3589661445
  CRPIX2  =  -1.76800252831
  CRVAL1  =  -0.0134278702942
  CRVAL2  =  -0.0444093066511
  CD1_1   =  0.00138877871593
  CD2_2   =  0.00138875543313
  CD1_2   =  1.41054781057e-05
  CD2_1   =  -1.41656423353e-05
  PV2_1   =  0.999993897433
  PV2_3   =  8.11292725428
  PV2_5   =  901.974288037


  nice = '/usr/bin/nice -n20 '

  casuin = 'IMAGE80620130715204739_cal.fits'

  casuin = 'catalog.fits'

#  casuin = '../combine_everything/IMAGE80620130606223430_cal.fits'

  mycatname = 'outputcat_dis.fits'

#  casuin = '../dithered_stack/outstack.fits'

#  casuin = '../combine_everything/IMAGE80620130802235108_cal.fits'
  casuin = '../combine_everything/IMAGE80620130711215053_cal.fits'
  casuin = '../combine_everything/IMAGE80620130801213346_cal.fits'

  extract = True
  fit = False
  thresh = 100

  prior_names = ['CRPIX1','CRPIX2','CRVAL1','CRVAL2','CD1_1','CD2_2','CD1_2','CD2_1','PV2_1','PV2_3','PV2_5']

#  prior = print_params(casuin)

  with pf.open(casuin) as hdulist:
    XVAL = hdulist[0].header['NAXIS1']/2
    YVAL = hdulist[0].header['NAXIS2']/2
    TEL_RA = hdulist[0].header['TEL_RA']
    TEL_DEC = hdulist[0].header['TEL_DEC']

  prior = [CRPIX1,CRPIX2,CRVAL1,CRVAL2,CD1_1,CD2_2,CD1_2,CD2_1,PV2_1,PV2_3,PV2_5]

  apply_correct(prior,casuin,XVAL,YVAL,TEL_RA,TEL_DEC)

  if extract == True:
    command = nice+'imcore '+casuin+' noconf outputcat_dis.fits 2 '+str(thresh)+' --noell'
    print command
    os.system(command)


  # this corrects efficiently for the 'offset' between the optical axis and the RA/DEC optical axis
  xs,ys,RA_sep,DEC_sep, sep_list = calc_seps(mycatname,casuin)
  prior[2] += median(RA_sep)
  prior[3] += median(DEC_sep)
  apply_correct(prior,casuin,XVAL,YVAL,TEL_RA,TEL_DEC)

  if extract == True:
    command = nice+'imcore '+casuin+' noconf outputcat_dis.fits 2 '+str(thresh)+' --noell'
    print command
    os.system(command)

#  prior = [CRPIX1,CRPIX2,CRVAL1,CRVAL2]


  if fit == True:
    final = opt.leastsq(find_dist_center, prior[0:2],args=(casuin,XVAL,YVAL,TEL_RA,TEL_DEC,prior[0],prior[1],prior[2],prior[3],thresh))
    prior = print_params(casuin)
    final = opt.leastsq(fit_distortion, prior[-2:],args=(casuin,XVAL,YVAL,TEL_RA,TEL_DEC,thresh))
  else:
    final = prior

  command = nice +'wcsfit '+casuin+' outputcat_dis.fits --site cds'
  print 'performing final wcs solution'
  print command
  os.system(command)

  oc = print_params(casuin)
  #the first 4 values are modifiers, not absoulte values

  if extract == True:
    command = nice+'imcore '+casuin+' noconf outputcat_dis.fits 2 '+str(thresh)+' --noell'
    print command
    os.system(command)

  filelist = 'stackfilelist'

#  filelist = '../combine_everything/filelist'

  plot_differences(mycatname,casuin)

  for l in open(filelist):
#    line = '../combine_everything/'+l.strip('\n').split('/')[-1].strip('.fits') + '_cal.fits'
    line = l.strip('\n')
    with pf.open(line.strip('\n')) as hdulist:
      XVAL = hdulist[0].header['NAXIS1']/2
      YVAL = hdulist[0].header['NAXIS2']/2
      TEL_RA = hdulist[0].header['TEL_RA']
      TEL_DEC = hdulist[0].header['TEL_DEC']
    print line.strip('\n')
    apply_correct(oc,line.strip('\n'),XVAL,YVAL,TEL_RA,TEL_DEC)

  create_source_cat(37,16,os.getcwd()+'/',conf_dir,True)

def print_params(casuin):
  oc = []

  names = ['CRPIX1','CRPIX2','CRVAL1','CRVAL2','CD1_1','CD2_2','CD1_2','CD2_1','PV2_1','PV2_3','PV2_5']

  with pf.open(casuin) as openfile:
    XLEN = openfile[0].header['NAXIS1']
    YLEN = openfile[0].header['NAXIS2']
    print '  CRPIX1',' = ', float(openfile[0].header['CRPIX1']) - XLEN/2
    oc += [float(openfile[0].header['CRPIX1']) - XLEN/2]
    print '  CRPIX2',' = ', float(openfile[0].header['CRPIX2']) - YLEN/2
    oc += [float(openfile[0].header['CRPIX2']) - YLEN/2]
    print '  CRVAL1',' = ', openfile[0].header['CRVAL1'] - openfile[0].header['TEL_RA']
    oc += [openfile[0].header['CRVAL1'] - openfile[0].header['TEL_RA']]
    print '  CRVAL2',' = ', openfile[0].header['CRVAL2'] - openfile[0].header['TEL_DEC']
    oc += [openfile[0].header['CRVAL2'] - openfile[0].header['TEL_DEC']]
    for name in names[4:]:
      print ' ',name,'  = ', openfile[0].header[name]
      oc += [openfile[0].header[name]]
  return oc

def find_dist_center(x,casuin,XVAL,YVAL,TEL_RA,TEL_DEC,x_prior,y_prior,CRVAL1_prior,CRVAL2_prior,thresh):
  # this should be the second step after the course fit, move only x and y. crvals are assumed to follow one to one. this should be good enough to center the distortion.

  nice = '/usr/bin/nice -n20 '

  apply_xy_correct(x,casuin,XVAL,YVAL,TEL_RA,TEL_DEC,x_prior,y_prior,CRVAL1_prior,CRVAL2_prior)
  
  sep = find_med_sep('outputcat_dis.fits',casuin)

  sep = mean(sep)

  print 'inputting: ',x
  print ' '

  oc = print_params(casuin)
  print 'Gave us mean sep: ',sep

  print ' '

  return [sep]*30

def apply_xy_correct(x,casuin,XVAL,YVAL,TEL_RA,TEL_DEC,x_prior,y_prior,CRVAL1_prior,CRVAL2_prior):

  x_shift = (x[0]-x_prior)
  y_shift = (x[1]-y_prior)

  plate_scale_y = 0.00138896050072
  plate_scale_x = 0.00138896050072/cos((TEL_DEC + y_shift*plate_scale_y + CRVAL2_prior)*pi/180.0)

  pf.setval(casuin,'CRPIX1',value=XVAL + x[0])
  pf.setval(casuin,'CRPIX2',value=YVAL + x[1])
  pf.setval(casuin,'CRVAL1',value=TEL_RA + x_shift*plate_scale_x + CRVAL1_prior)
  pf.setval(casuin,'CRVAL2',value=TEL_DEC + y_shift*plate_scale_y + CRVAL2_prior)

def apply_correct(x,casuin,XVAL,YVAL,TEL_RA,TEL_DEC):

  pf.setval(casuin,'CRPIX1',value=XVAL + x[0])
  pf.setval(casuin,'CRPIX2',value=YVAL + x[1])
  pf.setval(casuin,'CRVAL1',value=TEL_RA + x[2])
  pf.setval(casuin,'CRVAL2',value=TEL_DEC + x[3])
  pf.setval(casuin,'CD1_1',value=x[4])
  pf.setval(casuin,'CD2_2',value=x[5])
  pf.setval(casuin,'CD1_2',value=x[6])
  pf.setval(casuin,'CD2_1',value=x[7])
  if len(x) > 8:
    pf.setval(casuin,'PV2_1',value=x[8])
    pf.setval(casuin,'PV2_3',value=x[9])
    pf.setval(casuin,'PV2_5',value=x[10])

def apply_dist_correct(x,casuin,XVAL,YVAL,TEL_RA,TEL_DEC):

  pf.setval(casuin,'PV2_3',value=x[0])
  pf.setval(casuin,'PV2_5',value=x[1])

  
def plot_differences(mycatname,casuin):
  xs, ys, RA_sep, DEC_sep, sep_list = calc_seps(mycatname,casuin)

  print 'the median sep is',median(sep_list),' arcseconds'

  for i in range(0,len(xs)):
    plot([xs[i],(xs[i]+100000*RA_sep[i])],[ys[i],(ys[i]+100000*DEC_sep[i])],'k-')


  ylim([min(ys),max(ys)])
  xlim([min(xs),max(xs)])
  show()

def fit_distortion(x,casuin,XVAL,YVAL,TEL_RA,TEL_DEC,thresh):

  nice = '/usr/bin/nice -n20 '

  apply_dist_correct(x,casuin,XVAL,YVAL,TEL_RA,TEL_DEC)

  command = nice +'wcsfit '+casuin+' outputcat_dis.fits --site cds'
  print command
  os.system(command)
    
  sep = find_med_sep('outputcat_dis.fits',casuin)

  sep = mean(sep)

  print 'inputting: ',x
  print ' '

  oc = print_params(casuin)
  print 'Gave us mean sep: ',sep

  print ' '

  return [sep]*30

def calc_seps(mycatname,casuin):


  plate_scale = -3600.0/5.0

  mycat = pf.open(mycatname)

  cat_names = []
  RA_lims = []
  DEC_lims = []
  for line in open('catcache/index'):
    vals = line.strip('\n').split(' ')
    cat_names += [vals[0]]
    RA_lims += [[float(vals[2]),float(vals[3])]]
    DEC_lims += [[float(vals[4]),float(vals[5])]]

  n = 0

  cat_name = cat_names[n]
  cat = pf.open('catcache/'+cat_name)

  cat_RA_raw = cat[1].data['ra']
  cat_DEC_raw = cat[1].data['dec']

  zero = 21.5

  cat_Jmag = cat[1].data['Jmag']
  my_mag = zero - 2.512*log10(mycat[1].data['Aper_flux_3'])
    
  my_X = mycat[1].data['x_coordinate']
  my_Y = mycat[1].data['y_coordinate']

  pix_coords = [[my_X[i],my_Y[i]] for i in range(0,len(my_X))]

  world = load_wcs_from_file(casuin,pix_coords)

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
    my_RA = my_RA_raw[(my_RA_raw > RA_lims[n][0]) & (my_RA_raw < RA_lims[n][1]) & (my_DEC_raw > DEC_lims[n][0]) & (my_DEC_raw < DEC_lims[n][1])]
    my_DEC = my_DEC_raw[(my_RA_raw > RA_lims[n][0]) & (my_RA_raw < RA_lims[n][1]) & (my_DEC_raw > DEC_lims[n][0]) & (my_DEC_raw < DEC_lims[n][1])]
    cat_RA = cat_RA_raw[(cat_RA_raw > min(my_RA)) & (cat_RA_raw < max(my_RA)) & (cat_DEC_raw > min(my_DEC)) & (cat_DEC_raw < max(my_DEC))]
    cat_DEC = cat_DEC_raw[(cat_RA_raw > min(my_RA)) & (cat_RA_raw < max(my_RA)) & (cat_DEC_raw > min(my_DEC)) & (cat_DEC_raw < max(my_DEC))]
  except:
    return xs,ys,RA_sep,DEC_sep, array([10.0])

  my_X = my_X[(my_RA_raw > RA_lims[n][0]) & (my_RA_raw < RA_lims[n][1]) & (my_DEC_raw > DEC_lims[n][0]) & (my_DEC_raw < DEC_lims[n][1])]
  my_Y = my_Y[(my_RA_raw > RA_lims[n][0]) & (my_RA_raw < RA_lims[n][1]) & (my_DEC_raw > DEC_lims[n][0]) & (my_DEC_raw < DEC_lims[n][1])]

  my_brightest = argsort(my_mag[(my_RA_raw > RA_lims[n][0]) & (my_RA_raw < RA_lims[n][1]) & (my_DEC_raw > DEC_lims[n][0]) & (my_DEC_raw < DEC_lims[n][1])])[:20]

  c_b = argsort(cat_Jmag[(cat_RA_raw > min(my_RA)) & (cat_RA_raw < max(my_RA)) & (cat_DEC_raw > min(my_DEC)) & (cat_DEC_raw < max(my_DEC))])[:100]

  for i in my_brightest:
    RA = my_RA[i]
    DEC = my_DEC[i]
    sep = 3600*(((RA - cat_RA[c_b])*(cos(DEC*pi/180.0)))**2.0 + (DEC - cat_DEC[c_b])**2.0)**0.5
    index = argmin(sep)
    sep_list +=[sep[index]]
    RA_sep += [cat_RA[c_b][index] - RA]
    DEC_sep += [cat_DEC[c_b][index] - DEC]
    xs += [my_X[i]]
    ys += [my_Y[i]]
    x_sep += [RA_sep[-1]*plate_scale*cos(DEC*pi/180.0)]
    y_sep += [DEC_sep[-1]*plate_scale]

  clever = array([x for x in sort(sep_list)])
  
  course_fit = median(clever)

  if course_fit > 3:
    return xs,ys,RA_sep,DEC_sep, array(course_fit)
  
  for i in range(0,len(my_RA)):
    RA = my_RA[i]
    DEC = my_DEC[i]
    sep = 3600*(((RA - cat_RA)*(cos(DEC*pi/180.0)))**2.0 + (DEC - cat_DEC)**2.0)**0.5
    index = argmin(sep)
    sep_list +=[sep[index]]
    RA_sep += [cat_RA[index] - RA]
    DEC_sep += [cat_DEC[index] - DEC]
    xs += [my_X[i]]
    ys += [my_Y[i]]

  sep_list = array(sep_list)
  RA_sep = array(RA_sep)
  DEC_sep = array(DEC_sep)
  xs = array(xs)
  ys = array(ys)

  c = [sep_list < 3*5.0]

  return xs[c],ys[c],RA_sep[c],DEC_sep[c], sep_list[c]

def find_med_sep(mycatname,casuin):

  xs, ys, RA_sep, DEC_sep, sep_list = calc_seps(mycatname,casuin)

  med_sep = mean(sep_list)

  return med_sep
