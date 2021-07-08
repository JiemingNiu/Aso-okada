#!/usr/bin/env python
from forward_dc3d_aso import *

###############################
###############################
###############################
def point2posterior_mix(x,d,force,nstation,nchannel,input_xm=False,_xm_=None,output_syn=False):
	indices = np.arange(nstation*nchannel)
	i_bu = np.where(indices%nchannel==0)[0]
	i_be = np.where((indices%nchannel==3)|(indices%nchannel==4))[0]
	i_le = np.where((indices%nchannel==1)|(indices%nchannel==2))[0]
	n1 = len(i_bu)
	n2 = len(i_be)
	n3 = len(i_le)
	p_geo = x[:5]
	# + currently discarded
	p_sigma_bu = 1
	p_sigma_be = 1
	p_sigma_le = 1
	# - currently discarded
	p_lambda = 10**x[5]
	if force: p_lambda_f = 10**x[6]
	G = ok_point2syn(p_geo,nstation=nstation,nchannel=nchannel)
	if force:
		Gf = sf_point2syn(p_geo[:3],nstation=nstation,nchannel=nchannel)
		G = np.hstack((G,Gf))
	n, p = G.shape
	L = d[:,0].copy()
	D = d[:,1].copy()
	#######################
	Lnorm = np.abs(L)
	L /= Lnorm
	D /= Lnorm
	G /= np.tile(np.reshape(Lnorm, (n,1)), (1,p))
	Pvec = 10.*D[0]/D
	#######################
	Pvec[i_bu] /= p_sigma_bu
	Pvec[i_be] /= p_sigma_be
	Pvec[i_le] /= p_sigma_le
	P = np.diag(Pvec)
	Gp = np.dot(P, G)
	Lp = np.dot(P, L)
	if force:
		p1 = 4
		p2 = 3
		I = np.identity(p)
		I[0,0] /= p_lambda
		I[1,1] /= p_lambda
		I[2,2] /= p_lambda
		I[3,3] /= p_lambda
		I[4,4] /= p_lambda_f
		I[5,5] /= p_lambda_f
		I[6,6] /= p_lambda_f
		Gv = np.vstack((Gp, I))
		Lv = np.append(Lp, np.zeros(p))
	else:
		Gv = np.vstack((Gp, np.identity(p)/p_lambda))
		Lv = np.append(Lp, np.zeros(p))
	if input_xm: xm = _xm_
	else: xm = np.linalg.lstsq(Gv, Lv, rcond=-1)[0]
	v = np.dot(Gv, xm) - Lv
	fs = np.sum(v**2)
	if force:
		ln_posterior_positive = 2 * n1 * np.log(p_sigma_bu) + 2 * n2 * np.log(p_sigma_be) + 2 * n3 * np.log(p_sigma_le) + 2 * p1 * np.log(p_lambda) + 2 * p2 * np.log(p_lambda_f) + fs + np.log(np.linalg.det(np.dot(Gv.T,Gv)))
	else:
		ln_posterior_positive = 2 * n1 * np.log(p_sigma_bu) + 2 * n2 * np.log(p_sigma_be) + 2 * n3 * np.log(p_sigma_le) + 2 * p * np.log(p_lambda) + fs + np.log(np.linalg.det(np.dot(Gv.T,Gv)))
	if output_syn:
		G = ok_point2syn(p_geo,nstation=nstation,nchannel=nchannel)
		if force:
			Gf = sf_point2syn(p_geo[:3],nstation=nstation,nchannel=nchannel)
			G = np.hstack((G,Gf))
		syn = np.dot(G, xm)
		return (xm, syn)
	else: return ln_posterior_positive 
###############################
###############################
###############################
def ln_prior_bound(theta, bounds):
	for i in range(len(theta)):
		if theta[i] < bounds[i][0] or theta[i] > bounds[i][1]:
			prior = -np.inf
			break
		else:
			prior = 0.0
	return prior
