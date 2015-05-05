# -*- coding: utf-8 -*-

"""
Usage:
  super_sample.py [options] (<FILE_LIST>) (<OUT_NAME>)

Options:
  -h --help  Show help text
  --factor=FACTOR  What oversampling factor to use [default: 5]
  --size=SIZE  how large a region around each star to stack [default: 11]
  --stars=STARS  How many stars to stack in each quadrant [default: 100]
  --binning=BINNING  use a binning factor for long time series? [default: 1]

"""

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.io.fits as pf
import pickle
import threading
import scipy.optimize as opt
import os
import numpy as np

def super_sample(filelist,factor,size,stars,binning,tag):

  files = []
  for line in open(filelist,'r'):
    files += [line.strip('\n')]

  for file in files:
    call_find_fwhm(file,factor,size,stars,tag=tag)

 # labels = ['f_1','f_3','f_5','f_7','f_9']

#  for label in labels:
#    condense_data(label,tag=tag)

#  plot_everything(files[0:10],labels,binning,tag=tag)

def condense_data(label,tag=''):

  fwhm = []
  ellipse = []
  theta = []
  x = []
  y = []


  theta += pickle.load(open(tag+'theta'+label+'.p'))
  x += pickle.load(open(tag+'fwhm_x'+label+'.p'))
  y += pickle.load(open(tag+'fwhm_y'+label+'.p'))

  for i in range(0,len(x)):
    rad = theta[i]*np.pi/180.0
    theta[i] = theta[i]*np.pi/180.0
    theta[i] = theta[i] - 2*np.pi*(int(theta[i]/(2*np.pi)))
    print theta[i]

    if x[i] < y[i]:
      theta[i] -= np.pi/2.0

    if theta[i] < 0:
      theta[i] += 2*np.pi
    if theta[i] > np.pi:
      theta[i] -= np.pi

    theta[i] = theta[i]*180.0/np.pi

  #  print 'major?',abs(y[i]/np.cos(rad))

    ellipse_frame, average = find_minor_axis(rad,x[i],y[i])

#    A = 1.0

#    xp = arange(0,220)
#    yp = arange(0,220)

#    x0 = 110
#    y0 = 110

#    f = guassian2d(A,theta[i],x[i]*10,y[i]*10,x0,y0,xp,yp)
#    fig = plt.figure()
#    a = fig.add_subplot(1,1,1)
#    imshow(f)
#    a.set_ylim(0,220)

#    show()


    fwhm += [average]

    ellipse += [ellipse_frame]

  print 'the mean orientation of this star is:',mean(theta)

  pickle.dump(fwhm, open(tag+'_fwhm_'+str(label)+'.p','wb'))
  pickle.dump(ellipse, open(tag+'_ellipse_'+str(label)+'.p','wb'))
  pickle.dump(theta, open(tag+'_theta_'+str(label)+'.p','wb'))


def plot_everything(files,labels,binning,tag=''):

  i = 0

  fwhm = []
  ellipse = []
  theta = []

  for label in labels:
    fwhm += [pickle.load(open(tag+'_fwhm_'+str(label)+'.p'))]
    ellipse += [pickle.load(open(tag+'_ellipse_'+str(label)+'.p'))]
    theta += [pickle.load(open(tag+'_theta_'+str(label)+'.p'))]


  focus_position = []
  tel_alt = []
  mjd = []

  for line in files:
    with pf.open(line) as imdata:
      if 1.0*i/binning == int(i/binning):
	focus_position +=[imdata[0].header['FCSR_PHY']]
	tel_alt +=[imdata[0].header['TEL_ALT']]
	mjd +=[imdata[0].header['MJD']]
    i += 1
    if i > binning*(len(ellipse[0])-1): break

  subplot(2, 1, 1)
  for fwhms in fwhm:
    plot(mjd,fwhms,'o')
  ylabel('fwhm (pix)')
  xlabel('MJD')

  subplot(2, 1, 2)
  plot(mjd,tel_alt,'o')
  ylabel('Altitude')
  xlabel('MJD')
  ylim(0,100)

  savefig(tag+'_timeseries.png', bbox_inches=0)
  clf()

  for fwhms in fwhm:
    plot(tel_alt,fwhms,'o')
  ylabel('fwhm (pix)')
  xlabel('Altitude')
  xlim(0,100)
  savefig(tag+'_altitude.png', bbox_inches=0)
  clf()

  for thetas in theta:
    plot(mjd,thetas,'o') 
  ylabel('theta')
  xlabel('MJD')
  savefig(tag+'_theta.png', bbox_inches=0)
  clf()

  for ellipses in ellipse:
    plot(mjd,ellipses,'o') 
  ylabel('ellipticity')
  xlabel('MJD')
  savefig(tag+'_ellipse.png', bbox_inches=0)
  clf()

#for i in cross_sec_x:
#  plot(i)
#  show()

#im_ani.save('im.mp4',metadata={'artist':'Tom'})
#  i = 0
#  for fwhms in fwhm:
#    i += 1
#    plot(focus_position,array(fwhms) + i,'o') 
#  ylabel('fwhm (pix)')
#  xlabel('focus position (mm)')
#  savefig('plot_bin/'+tag+'_focuspos.png', bbox_inches=0)
#  clf()

def find_minor_axis(theta,sx,sy):

  # use math instead of numpy?

  a = ((np.cos(theta))**2.0)/(2*sx**2) + ((np.sin(theta))**2.0)/(2*sy**2)

  b = -np.sin(2.0*theta)/(4*sx**2.0) + sin(2.0*theta)/(4*sy**2.0)

  c = ((np.sin(theta))**2.0)/(2*sx**2) + ((np.cos(theta))**2.0)/(2*sy**2)

  w_1 = np.sqrt(np.log(2)) / sqrt(a + 2*b*tan(theta+np.pi/2.0) + c*(tan(theta+np.pi/2.0))**2.0)

  w_2 = np.sqrt(np.log(2)) / sqrt(a + 2*b*tan(theta) + c*(tan(theta))**2.0)

  if abs(theta) > 2.0*np.pi:
    theta = theta - 2.0*np.pi*int((theta / (2.0*np.pi)))

  if theta <0 :
    theta = theta + 2*np.pi


  w_1 = abs(w_1/np.sin(theta)) 
  w_2 = abs(w_2/np.cos(theta))

  if abs(w_1) > abs(w_2):
    w_major = w_1
    w_minor = w_2
  else:
    w_major = w_2
    w_minor = w_1

  ellipse = np.sqrt(1.0 - (w_minor / w_major)**2.0)

  average = (w_minor + w_major)/2.0

  xx = 0.1*np.arange(0,220)
  func = tan(theta+3*np.pi/2)
  yy = xx*func

  x0 = 11

  print ellipse, average

  return ellipse, average

def call_find_fwhm(file,factor,size,stars,tag=''):

  labels = ['f_1','f_3','f_5','f_7','f_9']

  fwhm_x = [[]]*len(labels)
  fwhm_y = [[]]*len(labels)
  fwhm = [[]]*len(labels)
  theta = [[]]*len(labels)

  zero_array = np.zeros((factor*size,factor*size))

  fwhm = {'f_1':[],'f_3':[],'f_5':[],'f_7':[],'f_9':[]}
  fwhm_x = {'f_1':[],'f_3':[],'f_5':[],'f_7':[],'f_9':[]}
  fwhm_y = {'f_1':[],'f_3':[],'f_5':[],'f_7':[],'f_9':[]}
  theta = {'f_1':[],'f_3':[],'f_5':[],'f_7':[],'f_9':[]}
  data = {'f_1':zero_array,'f_3':zero_array,'f_5':zero_array,'f_7':zero_array,'f_9':zero_array}
  lengths = {'f_1':True,'f_3':True,'f_5':True,'f_7':True,'f_9':True}

  f = plt.figure()

  label_no = 0

  for label in labels:
    fwhm_extract(file,factor,size,stars,tag)

    dat = pickle.load(open(tag+label+'.p','rb'))

    print dat

    data[label] += dat
    data[label] = data[label]/data[label].max()
    fwhm_x_frame, fwhm_y_frame, theta_frame = find_2dfwhm(data[label],factor,size)
    fwhm_x[label] += [fwhm_x_frame]
    fwhm_y[label] += [fwhm_y_frame]
    fwhm[label] += [(fwhm_x_frame + fwhm_y_frame)/2.0]
    theta[label] += [theta_frame]
    a = f.add_subplot(3,3,(label_no*2) + 1)
    ticks = factor*np.arange(size)
    a.set_yticks(ticks)
    a.set_yticklabels(np.arange(size))
    a.set_xticks(ticks)
    a.set_xticklabels(np.arange(size))
    reverse = (0,a.get_ylim()[1] + factor)
    a.set_ylim(reverse)
    cax = plt.imshow(data[label], interpolation='none')
    a.grid(True)
    center = factor*((size)/2.0)
    plt.plot(center,center,'gx')
    data[label] = zero_array
    lengths[label] = False
    label_no += 1
  plt.savefig(file.strip('.fits')+'_psf.png', bbox_inches=0)
  plt.clf()

  for label in labels:
    os.system('rm '+tag+label+'.p')

  cache = False

  if cache == True:
    for label in labels:
      pickle.dump(fwhm[label], open(tag+'fwhm'+label+'.p','wb'))
      pickle.dump(fwhm_x[label], open(tag+'fwhm_x'+label+'.p','wb'))
      pickle.dump(fwhm_y[label], open(tag+'fwhm_y'+label+'.p','wb'))
      pickle.dump(theta[label], open(tag+'theta'+label+'.p','wb'))

def find_2dfwhm(data,factor,size):
  x1 = [1.0,np.pi,0.6,0.4,size/2.0,size/2.0]
  x, success = opt.leastsq(guassian2d_fit, x1, args=(np.arange(1.0*factor*size)/factor,np.arange(1.0*factor*size)/factor,data))

  fwhm_x = 2.0*(np.sqrt(2.0*np.log(2.0)))*x[2]
  fwhm_y = 2.0*(np.sqrt(2.0*np.log(2.0)))*x[3]
  theta = x[1]*180.0/np.pi

  return fwhm_x, fwhm_y, theta

def guassian2d_fit(p,x,y,data):
  f = guassian2d(p[0],p[1],p[2],p[3],p[4],p[5],x,y)
  return (data - f).flatten()

def guassian2d(A,theta,sx,sy,x0,y0,x,y):
  xx, yy = np.meshgrid(x,y)

  a = ((np.cos(theta))**2.0)/(2*sx**2) + ((np.sin(theta))**2.0)/(2*sy**2)

  b = -np.sin(2.0*theta)/(4*sx**2.0) + np.sin(2.0*theta)/(4*sy**2.0)

  c = ((np.sin(theta))**2.0)/(2*sx**2) + ((np.cos(theta))**2.0)/(2*sy**2)

  f = A*np.exp(-(a*(x0-xx)**2 + 2*b*(x0- xx)*(y0 - yy) + c*(y0 - yy)**2.0))
  return f

def fwhm_extract(image_name,factor,size,stars,tag=''):

  with pf.open(image_name) as imdata:
    image = imdata[0].data

    size += 2

    with pf.open(image_name + '.phot') as photdata:
      mean_fluxes = photdata[1].data['Core3_flux']
      IQR = [(mean_fluxes < (np.median(mean_fluxes[mean_fluxes > np.median(mean_fluxes)]))) & (mean_fluxes > (np.median(mean_fluxes[mean_fluxes < np.median(mean_fluxes)])))]

      selection = [(np.argsort(mean_fluxes)[-(stars+100):-100])]

      xpos = photdata[1].data['X_coordinate']
      ypos = photdata[1].data['Y_coordinate']
      xpos = xpos[selection]
      ypos = ypos[selection]

#    with pf.open('../focus_test/test.cat') as photdata:
#      mean_fluxes = photdata[2].data['FLUX_APER']
#      IQR = [(mean_fluxes < (np.median(mean_fluxes[mean_fluxes > median(mean_fluxes)]))) & (mean_fluxes > (median(mean_fluxes[mean_fluxes < median(mean_fluxes)])))]
#      xpos = photdata[2].data['X_IMAGE'][mean_fluxes > 3*mean(mean_fluxes)]
#      ypos = photdata[2].data['Y_IMAGE'][mean_fluxes > 3*mean(mean_fluxes)]

    imdata.close()

    x_phase = np.array([abs(x-int(x)) for x in xpos])
    y_phase = np.array([abs(y-int(y)) for y in ypos])


#    n,bins,patches = hist(x_phase,10,normed=1)
#    show()
#    n,bins,patches = hist(y_phase,10,normed=1)
#    show()


#    top_left = [(xpos < mean(xpos)) & (ypos > mean(ypos))],xpos[(xpos < mean(xpos)) & (ypos > mean(ypos))]
#    top_right = [(xpos > mean(xpos)) & (ypos > mean(ypos))],xpos[(xpos > mean(xpos)) & (ypos > mean(ypos))]
#    bottom_left = [(xpos < mean(xpos)) & (ypos < mean(ypos))],xpos[(xpos < mean(xpos)) & (ypos < mean(ypos))]
#    bottom_right = [(xpos > mean(xpos)) & (ypos < mean(ypos))],xpos[(xpos > mean(xpos)) & (ypos < mean(ypos))]
#    condition_list = [top_left,top_right,bottom_left,bottom_right]
#    condition_name_list = ['top_left','top_right','bottom_left','bottom_right']

    mx = 2048
    my = 2048

    f_1 = [(xpos < mx/3) & (ypos > 2*my/3)]
    f_2 = [(xpos > mx/3) & (xpos < 2*mx/3) & (ypos < my/3)]
    f_3 = [(xpos > 2*mx/3) & (ypos > 2*my/3)]

    f_4 = [(xpos < mx/3) & (ypos < 2*my/3)]
    f_5 = [(xpos > mx/3) & (xpos < 2*mx/3) & (ypos < 2*my/3) & (ypos > mx/3)]
    f_6 = [(xpos > 2*mx/3) & (ypos < 2*my/3)]

    f_7 = [(xpos < mx/3) & (ypos < my/3)]
    f_8 = [(xpos > mx/3) & (xpos < 2*mx/3) & (ypos < my/3)]
    f_9 = [(xpos > 2*mx/3) & (ypos < my/3)]

    condition_list = [f_1,f_3,f_5,f_7,f_9]
    condition_name_list = ['f_1','f_3','f_5','f_7','f_9']

    for i in range(0,len(condition_list)):
      get_psf(ypos[condition_list[i][0]],xpos[condition_list[i][0]],image,size,factor,condition_name_list[i],tag)

#    threads = []
#    for i in range(0,len(condition_list)):
#      t = threading.Thread(target=get_psf, args = (ypos[condition_list[i][0]],xpos[condition_list[i][0]],image,size,factor,condition_name_list[i],tag))
#      threads.append(t)
#      threads = np.append(threads,t)
#    [x.start() for x in threads]
#    [x.join() for x in threads]


def get_psf(ypos,xpos,image,size,factor,condition_name,tag=''):

  a_size = size*factor
  stack = np.zeros((a_size,a_size))
  mini = 2*size
  maxi = 2048 -(2*size)
  condition = [(ypos > mini) & (ypos < maxi) & (xpos > mini) & (xpos < maxi)]

  xpos = np.array(xpos[condition])
  ypos = np.array(ypos[condition])

  for n in range(0,len(xpos)):
    coords = [ypos[n],xpos[n]]
    oversampled = return_sample_square(coords,image,size,factor)
    stack += oversampled
  # cleanup

  stack  = stack[factor:-factor,factor:-factor]
  stack = stack - np.median(stack)
  stack = stack / stack.max()

  outname = tag+condition_name+'.p'

  pickle.dump(stack, open(outname,'wb'))


def return_sample_square(c,image,size,factor):

  r = (size-1.0)/2.0
  raw = image[int(round(c[0])-r-1):int(round(c[0])+r),int(round(c[1])-r-1):int(round(c[1])+r)]

  yoffset = (round(c[1]) - c[1])*factor
  xoffset = (round(c[0]) - c[0])*factor

  oversampled = []
  for row in raw:
    oversampled_row = []
    for pixel in row:
      oversampled_row += [pixel]*factor
    oversampled += [oversampled_row]*factor
  oversampled = np.array(oversampled)
 
  # now recenter on the central pixel

  oversampled = recenter(oversampled,xoffset,yoffset)

  # now normalise the output

  maxval = oversampled.max() + 1.0

  oversampled = oversampled / maxval

  return oversampled

def recenter(oversampled,xshift,yshift):

  xshift = -int(round(xshift))
  yshift = -int(round(yshift))

  blanks = np.zeros((len(oversampled[:,]),abs(xshift))).T
  if xshift < 0:
    oversampled=np.append(blanks,oversampled[:xshift],axis=0)
  else:
    xshift = xshift
    blanks = np.zeros((len(oversampled[:,]),abs(xshift))).T
    oversampled=np.append(oversampled[xshift:],blanks,axis=0)

  blanks = np.zeros((len(oversampled[:,]),abs(yshift)))

  if yshift < 0:
    oversampled=np.append(blanks,oversampled[:,:yshift],axis=1)
  else:
    oversampled=np.append(oversampled[:,yshift:],blanks,axis=1)
  return oversampled

if __name__ == '__main__':

  description='''
  Produces diagnostic plot of the PSF by super-sampling and stacking the star images
  on the chip. Relies on good astrometry to function correctly.

  Currently requires a list of files and a casutools imcore output file. (note, not a imcore_list
  file)
  Searches for a 'catcache' directory in the same path as the image files, this is used for the
  comparison.
  Should modify so a direct path to the catanp.log file is given.
  '''

  parser = argparse.ArgumentParser(description=description)
  parser.add_argument('filelist')
  parser.add_argument('outname')
  parser.add_argument('--factor', type=int, required=False, default=5,
                      help='What oversampling factor to use [default: 5]')
  parser.add_argument('--size', type=int, required=False, default=11,
                      help='how large a region around each star to stack [default: 11]')
  parser.add_argument('--stars', type=int, required=False, default=100,
                      help='How many stars to stack in each quadrant [default: 100]')
  parser.add_argument('--binning', type=int, required=False, default=1,
                      help='use a binning factor for long time series? [default: 1]')

  args = parser.parse_args()

  super_sample(args.filelist,args.factor,args.size,args.stars,args.binning,args.outname)

# vim: sw=2
