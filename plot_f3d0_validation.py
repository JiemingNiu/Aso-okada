#!/usr/bin/env python
import matplotlib
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.size'] = 12 
matplotlib.rcParams['axes.linewidth'] = 0.5
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['ytick.right'] = True
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['xtick.major.size'] = 4
matplotlib.rcParams['xtick.minor.size'] = 1
matplotlib.rcParams['ytick.major.size'] = 4
matplotlib.rcParams['ytick.minor.size'] = 1
matplotlib.rcParams['mathtext.fontset']='custom'
matplotlib.rcParams['mathtext.rm']='Arial'
matplotlib.rcParams['mathtext.it']='Arial:italic'
matplotlib.rcParams['mathtext.bf']='Arial:bold'
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
import matplotlib.patches as patches
import numpy as np

validations = np.genfromtxt('f3d0_validation_log')

szz = np.unique(validations[:,0])
nszz = len(szz)
labels = ['Normalised vertical displacement [SI]','Normalised radial displacement [SI]','Normalised radial tilt [1000$\\times$SI]']
for i in range(nszz):
	fig = plt.figure(figsize=(8,4))
	axgrid = gridspec.GridSpec(1,3,bottom=0.15,left=0.1,right=0.97,top=0.9,wspace=0.5,hspace=0.25)
	#
	flag = np.where(validations[:,0]==szz[i])[0]
	tmp = validations[flag,:]
	#
	dist = tmp[:,1]/1000.
	mxd = np.max(tmp[:,2])
	for j in range(3):
		ax = plt.subplot(axgrid[0,j])
		dc3d0 = tmp[:,2+j]/mxd
		f3d0 = tmp[:,6+j]/mxd
		if j == 2:
			dc3d0 *= 1000
			f3d0 *= 1000
		ax.plot(dist, dc3d0, '-', color='0', linewidth=1.5, label='dc3d0',clip_on=False)
		ax.plot(dist, f3d0, '--', color='r', linewidth=1.5, dashes=(2,1), label='f3d0',clip_on=False)
		ax.set_ylabel(labels[j])
		ax.set_xlim([0,20])
		ax.set_xlabel('Distance [km]')
		ax.grid(linewidth=0.5,linestyle=':') 
		if j < 2: ax.set_ylim([0,1.1])
		else: ax.set_ylim([-1.1,0])
		if j == 1: ax.set_title('Source depth:%d km' % (szz[i]/1000),fontsize=12)
		if j == 0: ax.legend(frameon=0, ncol=1, loc=1, fontsize=12)
	plt.savefig('fig_validation_Z%d.jpg' % (szz[i]), dpi=300)
	plt.close()
	print i, nszz
	
