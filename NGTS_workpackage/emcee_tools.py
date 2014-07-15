# -*- coding: utf-8 -*-
import numpy as np
from emcee import EnsembleSampler
from fitsio import write as fwrite
from fitsio import read as fread
from fitsio import FITS as fFITS

# this is a simple front end for the emcee package for easy startup, featuring incremental saving and multithreading
# / mpi support.

def run_emcee(x,lnprob,args,nwalkers,nruns,fudge,chain_name,burns,pool=None,nthreads=1,namearray=[],resume=False,w=False):

  ndim = len(x)

  p0 = []

  if resume == True:
    p0, ndone = resume_file(chain_name,ndim,nwalkers)
    nruns  -= ndone
    n = (ndone + burns)/nwalkers
  else:
    for i in range(0,nwalkers):
      shuffle = (10**(fudge*(np.random.rand(ndim) - 0.5)))
      p0 += [list(shuffle*x)]
    initiate_file(chain_name,ndim,blob_list=namearray,w=w)
    n = 0


  iterations = int(nruns/nwalkers)

  if pool != None:
    sampler = EnsembleSampler(nwalkers,ndim,lnprob,args=args,pool=pool)
  else:
    sampler = EnsembleSampler(nwalkers,ndim,lnprob,args=args,threads=nthreads)

  for result in sampler.sample(p0, iterations=iterations, storechain=False):
      n += 1
      if (n > burns/nwalkers):
	position = result[0]
	logl = result[1]
	with fFITS(chain_name,'rw') as fits:
	  for k in range(position.shape[0]):
	    output = {'lp':np.array([logl[k]]),'x':np.array([position[k]])}
	    for i in range(0,len(namearray)):
	      blob = result[3][k][i]
	      output[namearray[i]] = np.array([blob])
	    if np.isfinite(logl[k]):
	      fits['MCMC'].append(output)
  pool.close()

def resume_file(chain_name,ndim,nwalkers):

  with fFITS(chain_name,'r') as fits:
    length = fits['MCMC']._info['nrows']
    p0 = fits['MCMC']['x'][length-nwalkers:length]

  return p0, length

def initiate_file(chain_name,ndim,blob_list=[],w=False):

  type_array = [('x','f8',(ndim)),('lp','f4')]
  for i in range(0,len(blob_list)):
    type_array += [(blob_list[i],'f4')]

  initiate = np.zeros(0, dtype=type_array)
  if w == False:
    print 'Not allowed to overwrite'
    quit()
  fwrite(chain_name,initiate,clobber=w,extname='MCMC')