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

def m_condense_data(filelist,nproc,appsize,verbose=False,outdir='./'):

  nfiles = 0
  for line in open(filelist):
    nfiles += 1

  starts, ends = thread_alloc(nfiles,nproc)

  processes = []
  for i in range(0,nproc):
    p = multiprocessing.Process(target=condense_data, args = (filelist,starts[i],ends[i],i+1,appsize,verbose))
    processes.append(p)
  [x.start() for x in processes]
  [x.join() for x in processes]

  filelist = array([ f for f in listdir(os.getcwd()) if 'output_' in f ])
  numberlist = array([int(f.split('_')[1].split('.')[0]) for f in filelist])
  ordered = numberlist.argsort() 
  stitch(filelist[ordered],outdir)

  for f in filelist:
    os.system('rm '+f)

def condense_data(filelist,minlen,maxlen,thread_no,appsize,verbose):

#Take all .phot outputs from casu imstack_ls and condense them into a single file with formatting suitible for reading by sysrem

  min_temp = 0
  max_clouds = 50
  max_fwhm = 10
  min_fwhm = 0.1
  max_shift = 100


  flux = []
  flux_err = []
  sky = []
  xpos = []
  ypos = []
  ALT = []
  AZ = []
  TEL_RA = []
  TEL_DEC = []
  time=[]
  centerra = []
  centerdec = []
  meanbias =[]
  stdbias = []
  skylevel =[]
  std_image = []
  T = []
  coolstat = []
  ADU_DEV = []
  fwhm = []
  SKY_MED = []
  CLOUDS = []
  SHIFT = []
  exposure = []
  seeing = []

  #the .copy() addition is necessary when dealing with long filelists - by default python list optimization keeps all files open otherwise,
  #leading to a crash from too many open files

  # an airmass correction term that I fitted myself to the data, might be best to remove this and let sysrem do it itself?
  k = 0.09624718

  pi =  3.14159265359

  npix = pi*appsize**2.0

  first_frame = True

  for i in range(minlen,maxlen):
    line = linecache.getline(filelist,i).strip('\n')
    status_checks = line.split(' ')[1:]
    image = line.split(' ')[0]
    if all([status == 'ok' for status in status_checks]):
      with pf.open(image+'.phot') as photdata:
        gain = photdata[1].header['GAINFACT']  
        try:
	  ambient = photdata[1].header['WXTEMP']
        except:
	  ambient = 30.0

        cloud_status = photdata[1].header['CLOUD_S']
        fwhm_frame = photdata[1].header['FWHM']
        frame_shift = photdata[1].header['SKY_MOVE']
	seeing_frame = photdata[1].header['SEEING']

        if ((cloud_status < max_clouds) & (fwhm_frame < max_fwhm) & (fwhm_frame > min_fwhm) & (frame_shift < max_shift)):
	  SHIFT += [frame_shift]
	  CLOUDS += [cloud_status]
	  SKY_MED += [photdata[1].header['SKYLEVEL']]
	  ALT +=[photdata[1].header['TEL_ALT']]
	  AZ +=[photdata[1].header['TEL_AZ']]
	  TEL_RA +=[photdata[1].header['TEL_RA']]
	  TEL_DEC +=[photdata[1].header['TEL_DEC']]
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
	  airmass = 1.0/cos((90.0-ALT[-1])*pi/180.0)
	  fluxcorrection = 10**(airmass*k/2.5)
	  sky += [photdata[1].data['Skylev'].copy()]
          frame_xpos = photdata[1].data['X_coordinate'].copy()
          frame_ypos = photdata[1].data['Y_coordinate'].copy()
	  xpos += [frame_xpos]
	  ypos += [frame_ypos]
	  utc = photdata[1].header['OBSSTART'].split('T')
	  yr, month, day = utc[0].split('-')
	  hr, min, sec = utc[1].split(':')
	  fwhm += [fwhm_frame]
	  seeing += [seeing_frame]
	  rawflux = photdata[1].data['Core_flux'].copy()
	  correctedflux = gain*rawflux*fluxcorrection/photdata[1].header['EXPOSURE']
	  flux += [correctedflux]
	  rel_err = 1.0/(rawflux*gain/sqrt(rawflux*gain + npix*sky[-1]*gain))
	  abs_err = rel_err*correctedflux
	  flux_err += [abs_err]
	  T +=[photdata[1].header['CCDTEMP']]
	  coolstat +=[photdata[1].header['COOLSTAT']]
          mjd = photdata[1].header['MJD']
	  time +=[[mjd]*len(flux[0])]
          if verbose == True:
	    print shape(time), line.split(' ')[0]+'.phot', thread_no
        else:
	  if verbose == True:  
            print 'Image ',image,' rejected for being too noisy! (sky SNR ',cloud_status,') (fwhm ',fwhm_frame,') (frame_shift ',frame_shift,')'

    else:
      print 'frame bad'
  # generate time of mid exposure array
  tarray = np.array(time)
  tmid = tarray[:,0] + ((array(exposure)/2.0)/(3600.0*24.0))


  #get some data that should be common to all frames
  with pf.open(image+'.phot') as photdata:
    RA = photdata[1].data['RA']
    DEC = photdata[1].data['DEC']

  zeros = tmid*0

  objid = np.arange(len(flux[0])) + 1

  fluxarray = np.array(flux).T
  flux_err_array = np.array(flux_err).T

  meanflux =[]
  for index in range (0,len(fluxarray[:,0])):
    meanflux += [np.mean(fluxarray[index,:])]

  meanarray = np.array(meanflux)

  npts = meanarray*0 + len(meanarray)

  c1 = pf.Column(name='OBJ_ID', format='26A', array=objid)
  c2 = pf.Column(name='FLUX_MEAN', format='1D', unit='Counts', array=meanarray)
  c3 = pf.Column(name='BLEND_FRACTION', format='1D', array=(meanarray)*0)
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

  hducatalogue=pf.new_table([c1,c2,c3,c4,c5,c6])

  hduimagelist=pf.new_table([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16])

  hduprime = pf.PrimaryHDU(np.array(flux).T)

  hduflux = pf.ImageHDU(np.array(flux).T)
  hdufluxerr = pf.ImageHDU(flux_err_array)
  hduxpos = pf.ImageHDU(np.array(xpos).T)
  hduypos = pf.ImageHDU(np.array(ypos).T)
  hdutime = pf.ImageHDU(np.array(time).T)
  hdunullq = pf.ImageHDU((np.array(time).T)*0 + 1)

  hdulist = pf.HDUList([hduprime] + [hducatalogue] + [hduimagelist] + [hdutime] + [hduflux] + [hdufluxerr] + [hdunullq] + [hduxpos] + [hduypos])

  hdulist[0].name = 'Primary'
  hdulist[1].name = 'CATALOGUE'
  hdulist[2].name = 'IMAGELIST'
  hdulist[3].name = 'HJD'
  hdulist[4].name = 'FLUX'
  hdulist[5].name = 'FLUXERR'
  hdulist[6].name = 'QUALITY'
  hdulist[7].name = 'CCDX'
  hdulist[8].name = 'CCDY'

  outname = 'output_'+str(thread_no)+'.fits'

  if len(exposure) > 0:
    hdulist.writeto(outname, clobber=True)

def stitch(filelist,outdir='./'):

  #combine the sub output files into a single master file, preserving the data types

  outname = outdir + '/output.fits'

  hdulist = []

  for filen in filelist:
    hdulist += [pf.open(filen)]

  headername ='FLUX'
  catalogue = hdulist[0][headername].data
  combine =[]
  for i in range(0,len(hdulist)):
    combine += list(hdulist[i][headername].data.T)
    print shape(combine)
  fluxmean = mean(combine, axis = 0)
  hduflux = pf.ImageHDU(array(combine).T)
  hduprime = pf.PrimaryHDU(array(combine).T)

  a = []
  headername ='IMAGELIST'
  for column in hdulist[0][headername].columns:
    combine =[]
    colname = column.name
    for i in range(0,len(hdulist)):
      combine += list(hdulist[i][headername].data[colname])
    a += [pf.Column(name=colname, format=column.format, array=combine)]

  c = []
  headername ='CATALOGUE'
  for column in hdulist[0][headername].columns:
    colname = column.name
    combine = list(hdulist[i][headername].data[colname])
    c += [pf.Column(name=colname, format=column.format, array=combine)]

  c[1].array = fluxmean

  headername_list = ['HJD','FLUXERR','QUALITY','CCDX','CCDY']
  dicty = {}

  for headername in headername_list:
    catalogue = hdulist[0][headername].data
    combine =[]
    for i in range(0,len(hdulist)):
      combine += list(hdulist[i][headername].data.T)
    dicty[headername] = pf.ImageHDU(array(combine).T)

  hduimagelist=pf.new_table(a)
  hducatalogue=pf.new_table(c)

  c1 = pf.Column(name='FLUX_MEAN', format='1D', unit='Counts', array=fluxmean)

  new_hdulist = pf.HDUList([hduprime] + [hducatalogue] + [hduimagelist] + [dicty['HJD']] + [hduflux] + [dicty['FLUXERR']] + [dicty['QUALITY']] + [dicty['CCDX']] + [dicty['CCDY']])

  new_hdulist[0].name = 'Primary'
  new_hdulist[1].name = 'CATALOGUE'
  new_hdulist[2].name = 'IMAGELIST'
  new_hdulist[3].name = 'HJD'
  new_hdulist[4].name = 'FLUX'
  new_hdulist[5].name = 'FLUXERR'
  new_hdulist[6].name = 'QUALITY'
  new_hdulist[7].name = 'CCDX'
  new_hdulist[8].name = 'CCDY'

  new_hdulist.writeto(outname, clobber=True)
