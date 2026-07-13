
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

sys.path.insert(0, '../../../')

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

anaUndu=uw.undulator.undulatorCharacterization(
	bEffY=b0_y,
	numPeriods=nperiods,
	periodLength=period_length,
	ebeam=ebeam,
	)

aspec=uw.ana_eb.aspectrum(ebeam=ebeam,undulator=anaUndu)
aspec.angular_power_distro_order(
	r=6,
	dX_scrn=0.04,
	dY_scrn=0.04,
	order=1,
	num_integration_pnts=201,
	omega=anaUndu.properFrequencyStrong.get(),
	sigmaOff=False,
	piOff=False,
	)
m0=anaUndu.firstHarmEnergyStrong.get()
aspec.flux_factor_comparisson(order=1)
aspec.onAxisAngularSpectralFluxDensity(
	ens=np.linspace(m0*0.5,m0*5.1,2000),
	orders=[1,3,5],
	freqs=None
	)

order=1
X,Y,angularPowerDistro,avrg_total_power=aspec.angularSpectralFluxDistro_order(
	r=6,
	dX_scrn=0.1,
	dY_scrn=0.1,
	order=order,
	num_integration_pnts=201,
	omega=anaUndu.properFrequencyStrong.get(),
	sigmaOff=False,
	piOff=False,
	)

fig = plt.figure(figsize=(13*uw.cm_inch, 6.5*uw.cm_inch), dpi=150)
fig.suptitle(f"Power Distribution \n Avrg Pow: {avrg_total_power:2f} W", fontsize=14)
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, angularPowerDistro, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
ax.set_xlabel('x [m]', fontsize=12)
ax.set_ylabel('y [m]', fontsize=12)
ax.set_zlabel('Power [W]', fontsize=12)
plt.savefig(dir_path+f"/power_distro_order_{order}.png"   , bbox_inches='tight')
plt.show()

pdb.set_trace()

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
