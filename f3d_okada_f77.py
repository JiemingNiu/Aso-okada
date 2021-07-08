#!/usr/bin/env python
import numpy as np
import okada92
from sympy import *
from sympy.functions.special.tensor_functions import KroneckerDelta

################################
# define variables
################################
# here F = F / (8*pi*mu)
alpha, lam, mu = symbols('ALPHA, LAMBDA, MU')
F = symbols('FP')
x1, x2, x3 = symbols('RX, RY, RZ')
y1, y2, y3 = symbols('SX, SY, SZ')

###############################
# define distance funcitons
###############################
def R(x1_, x2_, x3_, y1_, y2_, y3_):
	R1	=  x1_ - y1_
	R2	=  x2_ - y2_
	R3 	= -x3_ - y3_
	R	= sqrt(R1**2 + R2**2 + R3**2)
	return (R1, R2, R3, R)

##############################
# define displacment functions
# uj_i: ith component of displacement from jth component of force
#############################
def uA(j, i, x1_, x2_, x3_, y1_, y2_, y3_):
	Rs = R(x1_, x2_, x3_, y1_, y2_, y3_)
	return F * ((2-alpha)*KroneckerDelta(i,j)/Rs[3] + alpha * Rs[i] * Rs[j] / Rs[3]**3)	

def uB(j, i, x1_, x2_, x3_, y1_, y2_, y3_):
	Rs = R(x1_, x2_, x3_, y1_, y2_, y3_)
	A = KroneckerDelta(i,j)/(Rs[2]+Rs[3])
	B = (Rs[i]*KroneckerDelta(j,2)-Rs[j]*KroneckerDelta(i,2)*(1-KroneckerDelta(j,2)))/(Rs[3]*(Rs[3]+Rs[2]))
	C = Rs[i]*Rs[j]*(1-KroneckerDelta(i,2))*(1-KroneckerDelta(j,2))/(Rs[3]*(Rs[3]+Rs[2])**2)
	
	return 2 * F * (KroneckerDelta(i,j)/Rs[3] + Rs[i] * Rs[j] / Rs[3]**3 + \
		(1-alpha)/alpha * (A+B-C))

def uC(j, i, x1_, x2_, x3_, y1_, y2_, y3_):
	Rs = R(x1_, x2_, x3_, y1_, y2_, y3_)
	return 2 * F * (1-2*KroneckerDelta(i,2)) * \
		((2-alpha)*(Rs[i]*KroneckerDelta(j,2)-Rs[j]*KroneckerDelta(i,2))/Rs[3]**3 + \
		alpha * y3_ * (KroneckerDelta(i,j)/Rs[3]**3 - 3 * Rs[i] * Rs[j] / Rs[3]**5))

def okada92_eq1(j, i):
	return uA(j, i, x1, x2, -x3, y1, y2, y3) - uA(j, i, x1, x2, x3, y1, y2, y3) + \
		uB(j, i, x1, x2, x3, y1, y2, y3) + x3 * uC(j, i, x1, x2, x3, y1, y2, y3) 
#############################
#############################
#############################
print("# test function input ")
#############################
#############################
#############################
ex_w     = diff(okada92_eq1(0,2), y1) + diff(okada92_eq1(1,2), y2) + diff(okada92_eq1(2,2),y3)
ex_dwdx  = diff(ex_w, x1)
ex_dwdy  = diff(ex_w, x2)
fex_w    = lambdify((F, alpha, x1, x2, x3, y1, y2, y3), ex_w)
fex_dwdx = lambdify((F, alpha, x1, x2, x3, y1, y2, y3), ex_dwdx)
fex_dwdy = lambdify((F, alpha, x1, x2, x3, y1, y2, y3), ex_dwdy)
#############################
#############################
#############################
VP=1500.0
VS=800.0
RHO=1700.0
_MU_ = VS * VS * RHO
_LAMBDA_ = VP * VP * RHO  - 2 * _MU_
_ALPHA_ = (_LAMBDA_+_MU_)/(_LAMBDA_+2.0*_MU_)
_F_ = 1E18 / (8 * np.pi * _MU_)
_MEX_ = 1E18 / _MU_
#############################
#############################
#############################
rxx = np.arange(0,10001.0,1000.0)
szz = np.array([5000.0])

rxxs, szzs = np.meshgrid(rxx, szz)
rxxs = rxxs.flatten()
szzs = szzs.flatten()

npair = len(rxxs)

for ipair in range(npair):
	rxx = rxxs[ipair]
	szz = szzs[ipair]
	fex_uz  = fex_w   (_F_, _ALPHA_, rxx, 0.0, 0.0, 0.0, 0.0, -szz)
	fex_uzx = fex_dwdx(_F_, _ALPHA_, rxx, 0.0, 0.0, 0.0, 0.0, -szz)
	fex_uzy = fex_dwdy(_F_, _ALPHA_, rxx, 0.0, 0.0, 0.0, 0.0, -szz)
	ex_ux, ex_uy, ex_uz, ex_uxx, ex_uyxe, ex_uzx, ex_uxy, ex_uyy, ex_uzy, ex_uxz, ex_uyz, ex_uzz, iree = okada92.dc3d0(_ALPHA_, rxx, 0.0, 0.0, szz, 0.0, 0, 0, 0, _MEX_)
	print "%7d%7d%16.8e%16.8e%16.8e%16.8e%16.8e%16.8e" % (szz, rxx, ex_uz, fex_uz, ex_uzx, fex_uzx, ex_uzy, fex_uzy)
#############################
#############################
#############################
print("# test finished")
#############################
#############################
#############################
print("      SUBROUTINE  F3D0(ALPHA,RX,RY,RZ,DEPTH,FP,FX,\n\
     *               FY,FZ,UX,UY,UZ,UZX,UZY)\n\
      IMPLICIT NONE\n\
      DOUBLE PRECISION ALPHA,RX,RY,RZ,DEPTH,SZ,FP,\n\
     *                 FX,FY,FZ,UX,UY,UZ,UZX,UZY\n\
\n\
Cf2py intent(in) alpha,sx,sy,sz,depth,fp,fx,fy,fz\n\
Cf2py intent(out) ux,uy,uz,uzx,uzy\n\
\n\
      DOUBLE PRECISION U1_1, U2_1, U3_1, U1_2, U2_2,\n\
     *                 U3_2, U1_3, U2_3, U3_3\n\
      DOUBLE PRECISION U1_3_px1, U1_3_px2\n\
      DOUBLE PRECISION U2_3_px1, U2_3_px2\n\
      DOUBLE PRECISION U3_3_px1, U3_3_px2\n\
\n\
      SZ = -DEPTH\n\
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\
      ! code generated by sympy!\n\
      ! here FP = F/(8*pi*mu)\n\
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
y1 = 0.0
y2 = 0.0
for j in range(3):
	for i in range(3):
		uz  = okada92_eq1(j,i)
		uzx = diff(uz, x1)
		uzy = diff(uz, x2) 

		print(fcode(uz, assign_to='U%d_%d' % (j+1, i+1)))
		print ""
		if i == 2:
			print(fcode(uzx, assign_to='U%d_%d_px1' % (j+1, i+1)))
			print ""
			print(fcode(uzy, assign_to='U%d_%d_px2' % (j+1, i+1)))
			print ""
print ("      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\
      UX = U1_1 * FX + U2_1 * FY + U3_1 * FZ\n\
      UY = U1_2 * FX + U2_2 * FY + U3_2 * FZ\n\
      UZ = U1_3 * FX + U2_3 * FY + U3_3 * FZ\n\
      UZX = U1_3_px1 * FX + U2_3_px1 * FY + U3_3_px1 * FZ\n\
      UZY = U1_3_px2 * FX + U2_3_px2 * FY + U3_3_px2 * FZ\n\
      END")
