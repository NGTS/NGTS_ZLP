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
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.io.fits as pf
import pickle
import threading
import scipy.optimize as opt
import os
import numpy as np
import multiprocessing.dummy as multithreading
import multiprocessing
from functools import partial

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import itertools

matplotlib.rc('text', usetex=False)

def super_sample(filelist,inputcat,factor,size,stars,binning,tag,side,nproc=4):

    p = multiprocessing.Pool(nproc)

    # make curry

    files = []
    files_model = ''
    files_residuals = ''
    files_psf = ''
    files_surface = ''

    for line in open(filelist,'r'):
        files += [line.rstrip('\n')]
        files_psf += files[-1].rstrip('.fits')+'_psf.png,'
        files_model += files[-1].rstrip('.fits')+'_model.png,'
        files_residuals += files[-1].rstrip('.fits')+'_residuals.png,'
        files_surface += files[-1].rstrip('.fits')+'_surface.png,'
    
    data_points = len(files)

    f_5x = []
    f_5y = []
    theta = []
    mjd = []

    curry = []
    for i in range(0,data_points):
        curry += [[files[i],inputcat,factor,size,stars,tag,side]]
        with pf.open(files[i]) as imdata:
            mjd +=[imdata[0].header['MJD']]

    output = p.map(uncurry_call_find_fwhm,curry)

    for o in output:
        f_5x += o[0]['f_5']
        f_5y += o[1]['f_5']
        theta += o[2]['f_5']

    plt.plot(mjd,f_5x,'bo')
    plt.plot(mjd,f_5y,'ro')

    plt.savefig(tag+'_xy.png', bbox_inches=0)

    plt.close()

    plt.plot(mjd,theta,'ro')

    plt.savefig(tag+'_theta.png', bbox_inches=0)

    plt.close()

    fps = 5

    outputf = tag+'_psf.avi'
    command = 'mencoder mf://'+files_psf.rstrip(',')+' -mf w=800:h=600:fps='+str(fps)+':type=png -ovc raw -oac copy -o '+outputf
    print command
    os.system(command)

    outputf = tag+'_residuals.avi'
    command = 'mencoder mf://'+files_residuals.rstrip(',')+' -mf w=800:h=600:fps='+str(fps)+':type=png -ovc raw -oac copy -o '+outputf
    print command
    os.system(command)

    outputf = tag+'_model.avi'
    command = 'mencoder mf://'+files_model.rstrip(',')+' -mf w=800:h=600:fps='+str(fps)+':type=png -ovc raw -oac copy -o '+outputf
    print command
    os.system(command)

    outputf = tag+'_surface.avi'
    command = 'mencoder mf://'+files_surface.rstrip(',')+' -mf w=800:h=600:fps='+str(fps)+':type=png -ovc raw -oac copy -o '+outputf
    print command
    os.system(command)
  
def uncurry_call_find_fwhm(c):
    return call_find_fwhm(c[0],c[1],c[2],c[3],c[4],tag=c[5],side=c[6])


def call_find_fwhm(file,inputcat,factor,size,stars,tag='',side=3,label_filt=False):

    if(side % 2 == 0):
      print 'Side must be an odd number!'
      quit()

    tag = file.rstrip('.fits')

    if label_filt != False:
      labels = label_filt
      frame_center = int(label_filt[len(label_filt)/2].split('_')[-1])
    else:
      labels = []
      for i in range(0,side**2):
        labels += ['f_'+str(i+1)]
      frame_center = int((side**2 + 1)/2)


    zero_array = np.zeros((factor*size,factor*size))

    fwhm = {}
    fwhm_a = {}
    fwhm_b = {}
    theta = {}
    data = {}
    lengths = {}

    for label in labels:
      fwhm[label] = []
      fwhm_a[label] = []
      fwhm_b[label] = []
      theta[label] = []
      data[label] = zero_array.copy()
      lengths[label] = True

    f1 = plt.figure()
    f2 = plt.figure()
    f3 = plt.figure()

    stacks = fwhm_extract(file,inputcat,factor,size,stars,labels,side,tag)

    for label in labels:

        data[label] = stacks[label]/stacks[label].max()
        fwhm_a_frame, fwhm_b_frame, theta_frame, residuals, model = find_2dfwhm(data[label],factor,size)
        print fwhm_a_frame, fwhm_b_frame, theta_frame, label
        
        fwhm_a[label] += [fwhm_a_frame]
        fwhm_b[label] += [fwhm_b_frame]
        fwhm[label] += [(fwhm_a_frame + fwhm_b_frame)/2.0]
        theta[label] += [theta_frame]

        label_no = int(label.split('_')[-1])
        a1 = f1.add_subplot(side,side,label_no)
        ticks = factor*np.arange(size)
        a1.set_yticks(ticks)
        a1.set_yticklabels(np.arange(size))
        a1.set_xticks(ticks)
        a1.set_xticklabels(np.arange(size))
        reverse = (0,a1.get_ylim()[1] + factor)
        a1.set_ylim(reverse)
        cax = a1.imshow(data[label], interpolation='none',cmap='afmhot')
        a1.grid(True)
        center = factor*((size)/2.0)
        a1.plot(center,center,'gx')
        lengths[label] = False
        a1.get_xaxis().set_ticklabels([])
        a1.get_yaxis().set_ticklabels([])

        # ONLY COMPUTE THIS ONCE...
        if label_no == frame_center:
            stub = os.path.basename(file)
            av_fwhm_a = round(fwhm_a_frame,2)
            av_fwhm_b = round(fwhm_b_frame,2)
            eccentricity = round(np.sqrt(1.0-((av_fwhm_b/av_fwhm_a)**2)),2)
            av_theta = round(theta_frame,0)



        a2 = f2.add_subplot(side,side,label_no)
        ticks = factor*np.arange(size)
        a2.set_yticks(ticks)
        a2.set_yticklabels(np.arange(size))
        a2.set_xticks(ticks)
        a2.set_xticklabels(np.arange(size))
        reverse = (0,a2.get_ylim()[1] + factor)
        a2.set_ylim(reverse)
        cax = a2.imshow(abs(residuals), interpolation='none',cmap='afmhot',vmin=min(data[label].flatten()),vmax=max(data[label].flatten()))
        a2.grid(True)
        center = factor*((size)/2.0)
        a2.plot(center,center,'gx')
        lengths[label] = False
        a2.get_xaxis().set_ticklabels([])
        a2.get_yaxis().set_ticklabels([])


        a3 = f3.add_subplot(side,side,label_no)
        ticks = factor*np.arange(size)
        a3.set_yticks(ticks)
        a3.set_yticklabels(np.arange(size))
        a3.set_xticks(ticks)
        a3.set_xticklabels(np.arange(size))
        reverse = (0,a3.get_ylim()[1] + factor)
        a3.set_ylim(reverse)
        cax = a3.imshow(model, interpolation='none',cmap='afmhot',vmin=min(data[label].flatten()),vmax=max(data[label].flatten()))
        a3.grid(True)
        center = factor*((size)/2.0)
        a3.plot(center,center,'gx')
        data[label] = zero_array.copy()
        lengths[label] = False
        a3.get_xaxis().set_ticklabels([])
        a3.get_yaxis().set_ticklabels([])

    f1.suptitle('Center PSF: {0:.2f}        e: {1:.2f}      theta: {2:.0f}'.format(av_fwhm_a, eccentricity, av_theta),size=15)
    f2.suptitle('Center PSF: {0:.2f}        e: {1:.2f}      theta: {2:.0f}'.format(av_fwhm_a, eccentricity, av_theta),size=15)
    f3.suptitle('Center PSF: {0:.2f}        e: {1:.2f}      theta: {2:.0f}'.format(av_fwhm_a, eccentricity, av_theta),size=15)

    f1.text(0.5, 0.05, '{0}'.format(stub), ha="center", va="center",size=15)
    f2.text(0.5, 0.05, '{0}'.format(stub), ha="center", va="center",size=15)
    f3.text(0.5, 0.05, '{0}'.format(stub), ha="center", va="center",size=15)

    f1.savefig(file.rstrip('.fits')+'_psf.png', bbox_inches=0)
    f2.savefig(file.rstrip('.fits')+'_residuals.png', bbox_inches=0)
    f3.savefig(file.rstrip('.fits')+'_model.png', bbox_inches=0)

    plt.close()
    plt.close()
    plt.close()

    fit_focus_surface(labels,fwhm_a,fwhm_b,side,tag)

    cache = False
    if cache == True:
        for label in labels:
            pickle.dump(fwhm[label], open(tag+'fwhm'+label+'.p','wb'))
            pickle.dump(fwhm_a[label], open(tag+'fwhm_a'+label+'.p','wb'))
            pickle.dump(fwhm_b[label], open(tag+'fwhm_b'+label+'.p','wb'))
            pickle.dump(theta[label], open(tag+'theta'+label+'.p','wb'))

    return fwhm_a, fwhm_b, theta

def find_2dfwhm(data,factor,size):
    x1 = [1.0,np.pi,0.6,0.4,size/2.0,size/2.0]


    x1 = [1.0,np.pi,0.6,0.1,size/2.0,size/2.0]

    x = np.arange(1.0*factor*size)/factor
    y = np.arange(1.0*factor*size)/factor

    p, success = opt.leastsq(gaussian2d_fit, x1, args=(x,y,data))

    model = gaussian2d(p[0],p[1],p[2],p[2]+abs(p[3]),p[4],p[5],x,y)

    fwhm_b = 2.0*(np.sqrt(2.0*np.log(2.0)))*p[2]
    fwhm_a = 2.0*(np.sqrt(2.0*np.log(2.0)))*(p[2] + abs(p[3]))
    theta = p[1]*180.0/np.pi

    #while theta > 360:
        #theta -= 180

    #while theta < 0:
        #theta += 180

    #plt.imshow(model, interpolation='none')
    #plt.show()

    residuals = model -data

    return fwhm_a, fwhm_b, theta, residuals, model

def gaussian2d_fit(p,x,y,data):

    f = gaussian2d(p[0],p[1],p[2],p[2]+abs(p[3]),p[4],p[5],x,y)

    residuals = (data - f).flatten()/np.sqrt(abs(data.flatten()))

    residuals[abs(residuals) == np.inf] = 0.0

    return residuals

def gaussian2d(A,theta,sx,sy,x0,y0,x,y):

    # http://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function

    xx, yy = np.meshgrid(x,y)

    a = ((np.cos(theta))**2.0)/(2*sx**2) + ((np.sin(theta))**2.0)/(2*sy**2)

    b = -np.sin(2.0*theta)/(4*sx**2.0) + np.sin(2.0*theta)/(4*sy**2.0)

    c = ((np.sin(theta))**2.0)/(2*sx**2) + ((np.cos(theta))**2.0)/(2*sy**2)

    f = A*np.exp(-(a*(x0-xx)**2 + 2*b*(x0- xx)*(y0 - yy) + c*(y0 - yy)**2.0))
    return f

def fwhm_extract(image_name,inputcat,factor,size,stars,condition_name_list,side,tag=''):

    with pf.open(image_name) as imdata:
        image = imdata[0].data

        size += 2

        with pf.open(image_name + '.phot') as photdata:


            with pf.open(inputcat) as incat:
                mean_fluxes = incat[1].data['isophotal_flux']
                IQR = [(mean_fluxes < (np.median(mean_fluxes[mean_fluxes > np.median(mean_fluxes)]))) & (mean_fluxes > (np.median(mean_fluxes[mean_fluxes < np.median(mean_fluxes)])))]
                selection = [(np.argsort(mean_fluxes)[-(stars+100):-100])]

            xpos = photdata[1].data['X_coordinate']
            ypos = photdata[1].data['Y_coordinate']
            xpos = xpos[selection]
            ypos = ypos[selection]

        imdata.close()

        x_phase = np.array([abs(x-int(x)) for x in xpos])
        y_phase = np.array([abs(y-int(y)) for y in ypos])

        imdata.close()

        mx = 2048
        my = 2048
  
        condition_list = []

        i = 0
        for y in range(0,side)[::-1]:
          ymin = y*my/side
          ymax = (y+1)*my/side
          for x in range(0,side):
            xmin = x*mx/side
            xmax = (x+1)*mx/side
            condition = [(xpos > xmin) & (xpos < xmax) & (ypos < ymax) & (ypos > ymin)]
            condition_list += [condition]

        dat = {}

        for x in condition_name_list:
            i = int(x.split('_')[-1])-1
            stack = get_psf(ypos[condition_list[i][0]],xpos[condition_list[i][0]],image,size,factor,x,tag)
            dat[x] = stack

        return dat


def get_psf(ypos,xpos,image,size,factor,condition_name,tag=''):

    a_size = size*factor
    stack = np.zeros((a_size,a_size))
    mini = 2*size
    maxi = 2048 -(2*size)
    condition = [(ypos > mini) & (ypos < maxi) & (xpos > mini) & (xpos < maxi)]

    xpos = np.array(xpos[condition])
    ypos = np.array(ypos[condition])

    for n in range(0, len(xpos)):
        coords = [ypos[n], xpos[n]]
        oversampled = return_sample_square(coords, image, size, factor)
        stack += oversampled
    # cleanup

    stack  = stack[factor:-factor,factor:-factor]
    stack = stack - np.median(stack)
    stack = stack / stack.max()

    return stack


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

def fit_focus_surface(labels,fwhm_a,fwhm_b,side,tag):

  mx = 2048.0
  my = 2048.0

  condition_list = []

  i = 0

  xcens = []
  ycens = []

  for y in range(0,side)[::-1]:
    ycen = (y+0.5)*my/side
    for x in range(0,side):
      xcen = (x + 0.5)*mx/side
      xcens += [xcen]
      ycens += [ycen]

  xcens = np.array(xcens)
  ycens = np.array(ycens)

  Z = []
  for label in labels:
    psf_area = [fwhm_a[label][0]*fwhm_b[label][0]*np.pi]
    #Z += psf_area
    Z += [fwhm_b[label][0]]

  Z = np.array(Z)

  #order = np.argsort(ycens.copy())
  #xcens = xcens[order]
  #ycens = ycens[order]
  #Z = Z[order]

  X, Y = np.meshgrid(xcens, ycens)

  m = polyfit2d(xcens,ycens,Z,order=2)

  nx = 20
  ny = 20

  xx, yy = np.meshgrid(np.linspace(0, mx, nx), 
                        np.linspace(0, my, ny))
  zz = polyval2d(xx, yy, m)

  plt.imshow(zz, extent=(0, my, 0, mx), cmap=cm.afmhot)
  plt.scatter(xcens, ycens, c=Z, cmap=cm.afmhot)
  plt.savefig(tag+'_surface.png', bbox_inches=0)
  plt.close()

  pickle.dump(m, open(tag+'_fitparams_2.p','wb'))

  #fig = plt.figure()
  #ax = fig.gca(projection='3d')

  #surf = ax.plot_surface(xx, yy, zz, rstride=1, cstride=1, cmap=cm.coolwarm,
          #linewidth=0, antialiased=False)

  #ax.set_zlim(-1.01, 1.01)
  #ax.zaxis.set_major_locator(LinearLocator(10))
  #ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

  #fig.colorbar(surf, shrink=0.5, aspect=5)

  #plt.show()

def polyfit2d(x, y, z, order=3):
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, _, _, _ = np.linalg.lstsq(G, z)
    return m

def polyval2d(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z

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
    parser.add_argument('inputcat')
    parser.add_argument('outname')
    parser.add_argument('--factor', type=int, required=False, default=5,
                                            help='What oversampling factor to use [default: 5]')
    parser.add_argument('--size', type=int, required=False, default=11,
                                            help='how large a region around each star to stack [default: 11]')
    parser.add_argument('--stars', type=int, required=False, default=1000,
                                            help='How many stars to stack in each quadrant [default: 100]')
    parser.add_argument('--binning', type=int, required=False, default=1,
                                            help='use a binning factor for long time series? [default: 1]')
    parser.add_argument('--nproc', type=int, required=False, default=4,
                                            help='use multiprocessing? [default: 4]')
    parser.add_argument('--side', type=int, required=False, default=3,
                                            help='How many times should the image be split on each axis? (Must be odd) [default: 3]')

    args = parser.parse_args()

    super_sample(args.filelist,args.inputcat,args.factor,args.size,args.stars,args.binning,args.outname,args.side,args.nproc)

# vim: sw=2
