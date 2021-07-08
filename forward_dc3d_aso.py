#!/usr/bin/env python
import numpy as np
import okada92
import f3d0
from ll2xy import *

###############################################
###############################################
###############################################
vp = 1500.0
vs = 800.0
rho = 1700.0
mu = vs**2. * rho
lame1st = vp**2. * rho - 2 * mu
alpha = 1 - (vs/vp)**2.0
mss = 1.0E18
mds = 1.0E18
mic = 1.0E18
mex = 1.0E18
dss = 1.0
dds = 1.0
dts = 1.0
fp = 1.0E15
###############################################
###############################################
###############################################
def ok_fault2syn(parameters, nstation=4, nchannel=5):
	sxx, syy, szz, strike, dip, L_1, L_2, W_1, W_2 = parameters
	###############################################
	sxx *= 1000.0
	syy *= 1000.0
	szz *= 1000.0
	strike *= 180.0 / np.pi
	dip *= 180.0 / np.pi
	L_1 *= -1000.0
	L_2 *= 1000.0
	W_1 *= -1000.0
	W_2 *= 1000.0
	###############################################
	syn = np.zeros([nstation*nchannel,4])
	for ista in range(nstation):
		rlat = station_loc[ista][1]
		rlon = station_loc[ista][2]
		rzz1 = station_loc[ista][3]
		rzz2 = station_loc[ista][4]
		rxx, ryy = ll2xy(rlat, rlon)
		rxxn = rxx - sxx
		ryyn = ryy - syy
		#
		n1 = ista * nchannel
		#
		phi = np.deg2rad(strike)
		rxxp = np.sin(phi)*rxxn + np.cos(phi)*ryyn
		ryyp = -np.cos(phi)*rxxn + np.sin(phi)*ryyn
		depth = rzz1 + szz
		rzzp = rzz2 - rzz1
		#
		uxs, uys, uzs, uxxs, uyxs, uzxs, uxys, uyys, uzys, uxzs, uyzs, uzzs, irets = okada92.dc3d(alpha, rxxp, ryyp, 0.0, depth, dip, L_1, L_2, W_1, W_2, dss, 0.0, 0.0)
		uxd, uyd, uzd, uxxd, uyxd, uzxd, uxyd, uyyd, uzyd, uxzd, uyzd, uzzd, iretd = okada92.dc3d(alpha, rxxp, ryyp, 0.0, depth, dip, L_1, L_2, W_1, W_2, 0.0, dds, 0.0)
		uxc, uyc, uzc, uxxc, uyxc, uzxc, uxyc, uyyc, uzyc, uxzc, uyzc, uzzc, iretc = okada92.dc3d(alpha, rxxp, ryyp, 0.0, depth, dip, L_1, L_2, W_1, W_2, 0.0, 0.0, dts)
		#####
		uxe, uye, uze, uxxe, uyxe, uzxe, uxye, uyye, uzye, uxze, uyze, uzze, irete = okada92.dc3d0(alpha, rxxp, ryyp, 0.0, depth, dip, 0, 0, 0, mex/mu)
		###
		syn[n1,0] += uzs
		syn[n1,1] += uzd
		syn[n1,2] += uzc
		syn[n1,3] += uze
		#
		if nchannel > 3:
			syn[n1+3,0] += uxs*np.sin(phi) - uys*np.cos(phi)
			syn[n1+3,1] += uxd*np.sin(phi) - uyd*np.cos(phi)
			syn[n1+3,2] += uxc*np.sin(phi) - uyc*np.cos(phi)
			syn[n1+3,3] += uxe*np.sin(phi) - uye*np.cos(phi)
	
			syn[n1+4,0] += uys*np.sin(phi) + uxs*np.cos(phi)
			syn[n1+4,1] += uyd*np.sin(phi) + uxd*np.cos(phi)
			syn[n1+4,2] += uyc*np.sin(phi) + uxc*np.cos(phi)
			syn[n1+4,3] += uye*np.sin(phi) + uxe*np.cos(phi)
		
		uxs, uys, uzs, uxxs, uyxs, uzxs, uxys, uyys, uzys, uxzs, uyzs, uzzs, irets = okada92.dc3d(alpha, rxxp, ryyp, rzzp, depth, dip, L_1, L_2, W_1, W_2, dss, 0.0, 0.0)
		uxd, uyd, uzd, uxxd, uyxd, uzxd, uxyd, uyyd, uzyd, uxzd, uyzd, uzzd, iretd = okada92.dc3d(alpha, rxxp, ryyp, rzzp, depth, dip, L_1, L_2, W_1, W_2, 0.0, dds, 0.0)
		uxc, uyc, uzc, uxxc, uyxc, uzxc, uxyc, uyyc, uzyc, uxzc, uyzc, uzzc, iretc = okada92.dc3d(alpha, rxxp, ryyp, rzzp, depth, dip, L_1, L_2, W_1, W_2, 0.0, 0.0, dts)
		uxe, uye, uze, uxxe, uyxe, uzxe, uxye, uyye, uzye, uxze, uyze, uzze, irete = okada92.dc3d0(alpha, rxxp, ryyp, rzzp, depth, dip, 0, 0, 0, mex/mu)
		#
		duzsdx = uzxs*np.sin(phi) - uzys*np.cos(phi)
		duzsdy = uzys*np.sin(phi) + uzxs*np.cos(phi)
		#
		duzddx = uzxd*np.sin(phi) - uzyd*np.cos(phi)
		duzddy = uzyd*np.sin(phi) + uzxd*np.cos(phi)
		#
		duzcdx = uzxc*np.sin(phi) - uzyc*np.cos(phi)
		duzcdy = uzyc*np.sin(phi) + uzxc*np.cos(phi)
		#
		duzedx = uzxe*np.sin(phi) - uzye*np.cos(phi)
		duzedy = uzye*np.sin(phi) + uzxe*np.cos(phi)
		#
		syn[n1+1,0] += -duzsdx
		syn[n1+2,0] += -duzsdy
		#
		syn[n1+1,1] += -duzddx
		syn[n1+2,1] += -duzddy
		#
		syn[n1+1,2] += -duzcdx
		syn[n1+2,2] += -duzcdy
		#
		syn[n1+1,3] += -duzedx
		syn[n1+2,3] += -duzedy
		#
	return syn
###############################################
###############################################
###############################################
def ok_point2syn(parameters, nstation=4, nchannel=5):
	###############################################
	sxx, syy, szz, strike, dip = parameters
	###############################################
	sxx *= 1000.0
	syy *= 1000.0
	szz *= 1000.0
	strike *= 180.0 / np.pi
	dip *= 180.0 / np.pi
	###############################################
	syn = np.zeros([nstation*nchannel,4])
	for ista in range(nstation):
		rlat = station_loc[ista][1]
		rlon = station_loc[ista][2]
		rzz1 = station_loc[ista][3]
		rzz2 = station_loc[ista][4]
		rxx, ryy = ll2xy(rlat, rlon)
		rxxn = rxx - sxx
		ryyn = ryy - syy
		#
		n1 = ista * nchannel
		#
		phi = np.deg2rad(strike)
		rxxp = np.sin(phi)*rxxn + np.cos(phi)*ryyn
		ryyp = -np.cos(phi)*rxxn + np.sin(phi)*ryyn
		depth = rzz1 + szz
		rzzp = rzz2 - rzz1
		uxs, uys, uzs, uxxs, uyxs, uzxs, uxys, uyys, uzys, uxzs, uyzs, uzzs, irets = okada92.dc3d0(alpha, rxxp, ryyp, 0.0, depth, dip, mss/mu, 0, 0, 0)
		uxd, uyd, uzd, uxxd, uyxd, uzxd, uxyd, uyyd, uzyd, uxzd, uyzd, uzzd, iretd = okada92.dc3d0(alpha, rxxp, ryyp, 0.0, depth, dip, 0, mds/mu, 0, 0)
		uxc, uyc, uzc, uxxc, uyxc, uzxc, uxyc, uyyc, uzyc, uxzc, uyzc, uzzc, iretc = okada92.dc3d0(alpha, rxxp, ryyp, 0.0, depth, dip, 0, 0, mic/lame1st, 0)
		uxe, uye, uze, uxxe, uyxe, uzxe, uxye, uyye, uzye, uxze, uyze, uzze, irete = okada92.dc3d0(alpha, rxxp, ryyp, 0.0, depth, dip, 0, 0, 0, mex/mu)
		#
		syn[n1,0] += uzs
		syn[n1,1] += uzd
		syn[n1,2] += uzc
		syn[n1,3] += uze
		#
		if nchannel > 3:
			syn[n1+3,0] += uxs*np.sin(phi) - uys*np.cos(phi)
			syn[n1+3,1] += uxd*np.sin(phi) - uyd*np.cos(phi)
			syn[n1+3,2] += uxc*np.sin(phi) - uyc*np.cos(phi)
			syn[n1+3,3] += uxe*np.sin(phi) - uye*np.cos(phi)
	
			syn[n1+4,0] += uys*np.sin(phi) + uxs*np.cos(phi)
			syn[n1+4,1] += uyd*np.sin(phi) + uxd*np.cos(phi)
			syn[n1+4,2] += uyc*np.sin(phi) + uxc*np.cos(phi)
			syn[n1+4,3] += uye*np.sin(phi) + uxe*np.cos(phi)

		uxs, uys, uzs, uxxs, uyxs, uzxs, uxys, uyys, uzys, uxzs, uyzs, uzzs, irets = okada92.dc3d0(alpha, rxxp, ryyp, rzzp, depth, dip, mss/mu, 0, 0, 0)
		uxd, uyd, uzd, uxxd, uyxd, uzxd, uxyd, uyyd, uzyd, uxzd, uyzd, uzzd, iretd = okada92.dc3d0(alpha, rxxp, ryyp, rzzp, depth, dip, 0, mds/mu, 0, 0)
		uxc, uyc, uzc, uxxc, uyxc, uzxc, uxyc, uyyc, uzyc, uxzc, uyzc, uzzc, iretc = okada92.dc3d0(alpha, rxxp, ryyp, rzzp, depth, dip, 0, 0, mic/lame1st, 0)
		uxe, uye, uze, uxxe, uyxe, uzxe, uxye, uyye, uzye, uxze, uyze, uzze, irete = okada92.dc3d0(alpha, rxxp, ryyp, rzzp, depth, dip, 0, 0, 0, mex/mu)
		#
		duzsdx = uzxs*np.sin(phi) - uzys*np.cos(phi)
		duzsdy = uzys*np.sin(phi) + uzxs*np.cos(phi)
		#
		duzddx = uzxd*np.sin(phi) - uzyd*np.cos(phi)
		duzddy = uzyd*np.sin(phi) + uzxd*np.cos(phi)
		#
		duzcdx = uzxc*np.sin(phi) - uzyc*np.cos(phi)
		duzcdy = uzyc*np.sin(phi) + uzxc*np.cos(phi)
		#
		duzedx = uzxe*np.sin(phi) - uzye*np.cos(phi)
		duzedy = uzye*np.sin(phi) + uzxe*np.cos(phi)
		#
		syn[n1+1,0] += -duzsdx
		syn[n1+2,0] += -duzsdy
		#
		syn[n1+1,1] += -duzddx
		syn[n1+2,1] += -duzddy
		#
		syn[n1+1,2] += -duzcdx
		syn[n1+2,2] += -duzcdy
		#
		syn[n1+1,3] += -duzedx
		syn[n1+2,3] += -duzedy
	return syn
###############################################
###############################################
###############################################
def sf_point2syn(parameters, nstation=4, nchannel=5):
	###############################################
	sxx, syy, szz = parameters
	###############################################
	sxx *= 1000.0
	syy *= 1000.0
	szz *= 1000.0
	###############################################
	norm_fp = fp / (8.0*np.pi*mu)
	###############################################
	syn = np.zeros([nstation*nchannel,3])
	for ista in range(nstation):
		rlat = station_loc[ista][1]
		rlon = station_loc[ista][2]
		rzz1 = station_loc[ista][3]
		rzz2 = station_loc[ista][4]
		rxx, ryy = ll2xy(rlat, rlon)
		rxxp = rxx - sxx
		ryyp = ryy - syy
		n1 = ista * nchannel
		#
		depth = rzz1 + szz
		rzzp = rzz2 - rzz1
		ux1, uy1, uz1, uzx1, uzy1 = f3d0.f3d0(alpha, rxxp, ryyp, 0.0, depth, norm_fp, 1.0, 0.0, 0.0)
		ux2, uy2, uz2, uzx2, uzy2 = f3d0.f3d0(alpha, rxxp, ryyp, 0.0, depth, norm_fp, 0.0, 1.0, 0.0)
		ux3, uy3, uz3, uzx3, uzy3 = f3d0.f3d0(alpha, rxxp, ryyp, 0.0, depth, norm_fp, 0.0, 0.0, 1.0)
		#
		syn[n1,   0] = uz1
		syn[n1,   1] = uz2
		syn[n1,   2] = uz3
		#
		if nchannel > 3:
			syn[n1+3, 0] = ux1
			syn[n1+3, 1] = ux2
			syn[n1+3, 2] = ux3
			syn[n1+4, 0] = uy1
			syn[n1+4, 1] = uy2
			syn[n1+4, 2] = uy3
		ux1, uy1, uz1, uzx1, uzy1 = f3d0.f3d0(alpha, rxxp, ryyp, rzzp, depth, norm_fp, 1.0, 0.0, 0.0)
		ux2, uy2, uz2, uzx2, uzy2 = f3d0.f3d0(alpha, rxxp, ryyp, rzzp, depth, norm_fp, 0.0, 1.0, 0.0)
		ux3, uy3, uz3, uzx3, uzy3 = f3d0.f3d0(alpha, rxxp, ryyp, rzzp, depth, norm_fp, 0.0, 0.0, 1.0)
		#
		syn[n1+1, 0] = -uzx1
		syn[n1+2, 0] = -uzy1
		#
		syn[n1+1, 1] = -uzx2
		syn[n1+2, 1] = -uzy2
		#
		syn[n1+1, 2] = -uzx3
		syn[n1+2, 2] = -uzy3
		#
	return syn 
