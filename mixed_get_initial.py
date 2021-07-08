#!/usr/bin/env python
import numpy as np
import pickle
import sys
from forward_dc3d_aso import *
from timeit import default_timer as timer
from scipy.optimize import differential_evolution
from multiprocessing import Pool
from PyMAP import *
from mpi4py import MPI
import emcee
import gc

comm = MPI.COMM_WORLD
iprc = comm.Get_rank()
nprc = comm.Get_size()

###############################################
###############################################
###############################################
evt = sys.argv[1]
parameterisation = sys.argv[2]
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
if 'noi' not in evt: infile = 'Offset_ref_%s' % evt
else: infile = 'syn_offset_%s.txt' % evt
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
Niter = 100 
sbin = int(np.ceil(float(Niter)/float(nprc)))
bins = []
for i in range(nprc-1):
	bins.append(np.arange(sbin*i, sbin*i+sbin))
bins.append(np.arange(sbin*(nprc-1), Niter))
###############################
###############################
###############################
if parameterisation == 'mt':
	bounds = [(-10,10),(-10,10),(0,10),(0,2*np.pi),(0,np.pi/2),(-10,10)]
	force = False
	nx = 4
else:
	bounds = [(-10,10),(-10,10),(0,10),(0,2*np.pi),(0,np.pi/2),(-10,10),(-10,10)]
	force = True
	nx = 7
ndim = len(bounds)
np.random.seed(3145)
xseeds = np.random.randint(9999,size=1000)
###############################
###############################
###############################
nbin = len(bins[iprc])
initial = np.zeros([nbin,ndim+nx])
for kiter in range(nbin):
	iseed = bins[iprc][kiter]
	de_solution = differential_evolution(point2posterior_mix, bounds, args=(data,force,nstation,nchannel), maxiter=2000, popsize=200, mutation=1.5, recombination=0.5, polish=True, disp=False, seed=xseeds[iseed])
	###############################
	xm, syn = point2posterior_mix(de_solution['x'], data, force, nstation, nchannel, output_syn=True)
	###############################
	x_mean = np.append(de_solution['x'], xm)
	###############################
	initial[kiter,:] = x_mean.copy()
###############################
###############################
###############################
if iprc > 0:
	comm.Send(initial, dest=0, tag=9999)
else:
	initials = np.zeros([Niter,ndim+nx])
	initials[bins[0],:] = initial.copy()
	for i in range(1, nprc):
		initial = np.zeros([len(bins[i]),ndim+nx])
		comm.Recv(initial, source=i, tag=9999)
		initials[bins[i],:] = initial.copy()
	for i in range(Niter):
		initial_str = "%s,%d" % (parameterisation, i)
		x_mean = initials[i,:]
		for i in range(len(x_mean)):
			initial_str += ",%.3e" % x_mean[i]
		print(initial_str)
