#!/usr/bin/env python
import numpy as np
from ll2xy import *
import f3d0
import okada92

VP=1500.0
VS=800.0
RHO=1700.0
MU = VS * VS * RHO
LAMBDA = VP * VP * RHO  - 2 * MU
ALPHA = (LAMBDA+MU)/(LAMBDA+2.0*MU)

F   = 1E18 / (8 * np.pi * MU)
MEX = 1E18 / MU

rxx = np.arange(0,20001.0,100.0)
szz = np.arange(1000.0,10001.0,1000.0)

rxxs, szzs = np.meshgrid(rxx, szz)
rxxs = rxxs.flatten()
szzs = szzs.flatten()

npair = len(rxxs)
#dh = 0.001

for ipair in range(npair):
	rxx = rxxs[ipair]
	szz = szzs[ipair]
	#####
	### by DC3D0 for reference
        #####
	ex_ux, ex_uy, ex_uz, ex_uxx, ex_uyxe, ex_uzx, ex_uxy, ex_uyy, ex_uzy, ex_uxz, ex_uyz, ex_uzz, irete = okada92.dc3d0(ALPHA, rxx, 0.0, 0.0, szz, 0.0, 0, 0, 0, MEX)

	####
	### by f3d0 for comparison
	####
	#dhx = 0.001 * rxx
	#dhz = 0.001 * szz 
	dhx = 1E-2 
	dhz = 1E-2
	fx1_ux, fx1_uy, fx1_uz, fx1_uzx, fx1_uzy = f3d0.f3d0(ALPHA, rxx+dhx, 0.0,  0.0, szz,     F, 1.0, 0.0, 0.0)	
	fx2_ux, fx2_uy, fx2_uz, fx2_uzx, fx2_uzy = f3d0.f3d0(ALPHA, rxx-dhx, 0.0,  0.0, szz,     F, 1.0, 0.0, 0.0)	
	####
	fy1_ux, fy1_uy, fy1_uz, fy1_uzx, fy1_uzy = f3d0.f3d0(ALPHA, rxx,     dhx,  0.0, szz,     F, 0.0, 1.0, 0.0)	
	fy2_ux, fy2_uy, fy2_uz, fy2_uzx, fy2_uzy = f3d0.f3d0(ALPHA, rxx,     -dhx, 0.0, szz,     F, 0.0, 1.0, 0.0)	
	####
	fz1_ux, fz1_uy, fz1_uz, fz1_uzx, fz1_uzy = f3d0.f3d0(ALPHA, rxx,     0.0,  0.0, szz+dhz, F, 0.0, 0.0, 1.0)	
	fz2_ux, fz2_uy, fz2_uz, fz2_uzx, fz2_uzy = f3d0.f3d0(ALPHA, rxx,     0.0,  0.0, szz-dhz, F, 0.0, 0.0, 1.0)	
	####
	fex_ux  = (fx2_ux  - fx1_ux  + fy2_ux  - fy1_ux  + fz2_ux  - fz1_ux)/2./dhx
	fex_uy  = (fx2_uy  - fx1_uy  + fy2_uy  - fy1_uy  + fz2_uy  - fz1_uy)/2./dhx
	fex_uz  = (fx2_uz  - fx1_uz  + fy2_uz  - fy1_uz  + fz2_uz  - fz1_uz)/2./dhz
	fex_uzx = (fx2_uzx - fx1_uzx + fy2_uzx - fy1_uzx + fz2_uzx - fz1_uzx)/2./dhx
	fex_uzy = (fx2_uzy - fx1_uzy + fy2_uzy - fy1_uzy + fz2_uzy - fz1_uzy)/2./dhx
	#####
	print "%7d%7d%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e" % (szz, rxx, ex_uz, ex_ux, ex_uzx, ex_uzy, fex_uz, fex_ux, fex_uzx, fex_uzy)

