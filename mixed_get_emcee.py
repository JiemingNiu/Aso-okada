#!/usr/bin/env python
import numpy as np
import pickle
import sys
from forward_dc3d_aso import *
from timeit import default_timer as timer
from scipy.optimize import differential_evolution
from multiprocessing import Pool
from PyMAP import *
import emcee
import gc

evt = sys.argv[1]
###############################################
###############################################
###############################################
tic = timer()
###############################################
###############################################
###############################################
stations = ['ASHV','ASIV','ASNV','ASTV']
channels = ['BU','LE','LN','BE','BN']
nstation = len(stations)
if 'C5' in evt: nchannel = 5
else: nchannel = 3
###############################################
###############################################
###############################################
ichns = {'BU':0,'LE':1,'LN':2,'BE':3,'BN':4}
jstas = {'ASHV':0,'ASIV':1,'ASNV':2,'ASTV':3,'HND':4,'TKD':5,'TMC':6}
###############################################
###############################################
###############################################
infile = 'Offset_ref_%s' % evt
###############################################
###############################################
###############################################
data = np.zeros([nstation*nchannel,2])
with open(infile,'r') as f:
	for line in f.readlines():
		if line[0] == '#': continue
		lines = line.split()
		sta = lines[0].split('.')[1]
		chn = lines[0].split('.')[2]
		amp = float(lines[1])
		damp = float(lines[2])
		ichn = ichns[chn]
		jsta = jstas[sta]
		n1 = jsta*nchannel+ichn
		data[n1,0] = amp
		data[n1,1] = damp
###############################
###############################
###############################
def get_bound(parameterisation):
	if parameterisation == 'mt':
		bounds = [(-10,10),(-10,10),(0,10),(0,2*np.pi),(0,np.pi/2),(-5,5)]
		force = False
		nx = 4
	else:
		bounds = [(-10,10),(-10,10),(0,10),(0,2*np.pi),(0,np.pi/2),(-5,5),(-5,5)]
		force = True
		nx = 7
	return (nx, force, bounds)
###############################
###############################
###############################
def log_probability(theta, data2d=None, bounds=None, force=True, nx=7):
	p_geo = theta[:-nx]
	p_slip = theta[-nx:]
	lp = ln_prior_bound(p_geo,bounds)
	if not np.isfinite(lp): return -np.inf
	return lp - point2posterior_mix(p_geo, data2d, force, nstation, nchannel, input_xm=True, _xm_=p_slip)
###############################
###############################
###############################
np.random.seed(3145)
xseeds = np.random.randint(9999,size=1000)
###############################
###############################
###############################
pool = Pool()
with open('initial_mixed_%s.txt' % evt,'r') as fin:
	for line in fin.readlines():
		lines = line.split(',')
		parameterisation = lines[0]
		pid = int(lines[1])
		nx, force, bounds = get_bound(parameterisation)
		x_mean = []
		for i in range(2,len(lines)):
			x_mean.append(float(lines[i]))
		x_mean = np.array(x_mean)
		x_std = np.zeros_like(x_mean)
		x_std[:-nx] += 1.0
		x_std[-nx:] += np.abs(x_mean[-nx:])
		print " -- Iteration [%s]: %d " % (parameterisation, pid)
		print " -- Parameter mean:", x_mean
		print " -- Parameter std:", x_std
		ndim = len(x_mean)
		nwalk = 100
		nstep = 1E5
		initials = np.zeros([nwalk, ndim])
		for i in range(ndim):
			np.random.seed(xseeds[i+pid*ndim])
			initials[:,i] = np.random.normal(x_mean[i], x_std[i], size=nwalk)
		###############################
		###############################
		###############################
		sampler = emcee.EnsembleSampler(nwalk, ndim, log_probability, kwargs={'data2d':data, 'bounds':bounds, 'force':force, 'nx':nx}, pool=pool)
		sampler.run_mcmc(initials, nstep, progress=True)
		###############################
		samples = sampler.get_chain(thin=50)
		###############################
		###############################
		with open('mcmc_mixed_output_%s_PID%d_%s.pkl' % (evt,pid,parameterisation), 'w') as f:
			pickle.dump([data, x_mean, initials, samples, bounds], f)
		###############################
toc = timer()
print 'Elapsed time [s]:', toc-tic
print('###############################')
gc.collect()
