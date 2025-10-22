
import os
import pdb
import sys
import math
import numpy as np

try :
	# works when calling script with python3 script_file
	dir_path = os.path.dirname(os.path.realpath(__file__))
except:
	# works when calling script with exec from python console
	dir_path = os.getcwd()

import unduwave as uw
import unduwave.constants as uc
from unduwave.unduwave_incl import *

wave = uw.wave()
ebeam=wave._ebeam_paras
ebeam.get_std_bessy_III_paras()
ebeam.beam_en.set( ebeam.beam_en.get()*0.001 )
ebeam.update_values()

undu_paras = wave._undu_paras # getting parameter object

period_length = 0.02 # m
nperiods = 250
b0_y = 1

aundu=uw.anas.undulator(
	bEff=b0_y,
	periodLength=period_length,
	numPeriods=nperiods,
	lengthEndPeriodsRelative=0.0,
	ebeam=ebeam,
	thetaObservation=0.0,
	unduK=None
	)

aspec=uw.ana_eb.aspectrum(ebeam=ebeam,undulator=aundu)

num_integration_pnts = 100
r = 6 # m
dX_scrn = 0.01 # m
dY_scrn = 0.01 # m

xs = np.linspace(-dX_scrn/2.0,dX_scrn/2.0,num_integration_pnts)
ys = np.linspace(-dY_scrn/2.0,dY_scrn/2.0,num_integration_pnts)
zs = [r]
coord_tensor = uw.ana_eb.spherical_coord_from_xyz(xs,ys,zs)

shape_tensor = ( len(xs), len(ys) )
distro_flux = np.zeros(shape_tensor)

# for ind_x, x in enumerate(xs) : 
# 	for ind_y, y in enumerate(ys) : 
# 		distro_flux[ind_x,ind_y] = \
# 			aspec.flux_strong_undu(
# 				x=x,
# 				y=y,
# 				z=10,
# 				omega=425.612*uc.q_el/uc.hbar,
# 				order=1,
# 				)

# fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
# fig.suptitle("Angular Distribution", fontsize=14)
# ax = plt.gca()
# # ax.contour(theta_vals, phi_vals, distro_flux)
# ax = fig.add_subplot(111, projection='3d')


# X, Y = np.meshgrid(xs,ys)
# surf = ax.plot_surface(X, Y, distro_flux, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)

# ax.set_xlabel('x', fontsize=12)
# ax.set_ylabel('y', fontsize=12)
# # ax.set_zlabel('Flux', fontsize=12)
# plt.savefig(f"flux_distro.png"   , bbox_inches='tight')
# plt.show()

# pdb.set_trace()

theta_vals = np.linspace(-1/ebeam.gammaFactor.get(),1/ebeam.gammaFactor.get(),100)
phi_vals = np.linspace(0,2*math.pi,100)

[ zg_t, xg_t, yg_t ] = uw.ana_eb.make_coords( 
	fun = aspec.fun_transv, 
	coords1 = theta_vals, 
	coords2 = phi_vals 
	)

[ zg_ur, xg_ur, yg_ur ] = uw.ana_eb.make_coords( 
	fun = aspec.fun_transv_ultra_rel, 
	coords1 = theta_vals, 
	coords2 = phi_vals 
	)

mlab.mesh(xg_t,yg_t,zg_t)
# Simple plot.
mlab.xlabel("x")
mlab.ylabel("y")
mlab.zlabel("z")

mlab.figure()

mlab.mesh(xg_ur,yg_ur,zg_ur)
# Simple plot.
mlab.xlabel("x")
mlab.ylabel("y")
mlab.zlabel("z")

mlab.show()


pdb.set_trace()
