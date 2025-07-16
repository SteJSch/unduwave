
import os
import pdb
import sys
import math
import numpy as np
sys.path.insert(0, '../../../')

try :
	# works when calling script with python3 script_file
	dir_path = os.path.dirname(os.path.realpath(__file__))
except:
	# works when calling script with exec from python console
	dir_path = os.getcwd()

import unduwave as uw

field_folder = f'/'
res_folder = dir_path+'/res/'

"""
Getting wave
"""
wave = uw.wave(wave_mode='undu_endp') # Simple undulator model with endpoles
wave_prog_paras = wave._wave_prog_paras

"""
Setting Undulator Parameter

We give the K-parameters and the B-Amplitudes are calculated from there
"""

K_para_y = np.linspace(0.5,2.34,10) # Possible K-Values in horizontal (z) mode
K_para_y=K_para_y[0] 

period_length = 0.021 # m
nperiods = 238

wave_num_k = 2 * math.pi / ( period_length )
b0_y = K_para_y*uw.m_el * uw.v_c * wave_num_k/uw.q_el

undu_paras = wave._undu_paras # getting parameter object
undu_paras.b0Halbasy.set(b0_y)
undu_paras.ahwpolHalbasy.set(2*nperiods+1) # we count the number of B-field peaks here - one extra for the end-fields (odd)
undu_paras.xlHalbasy.set(period_length)

"""
Setting Program Parameter
"""
wave_prog_paras.res_folder.set(res_folder)
wave_prog_paras.spec_calc.set(True)
wave_prog_paras.nthreads.set(6)

"""
Setting Spectrometer Parameter
"""

spectrometer_paras = wave._spectrometer_paras
spectrometer_paras.freq_num.set(100) 
spectrometer_paras.freq_low.set(100) # lower spectrum energie in eV
spectrometer_paras.freq_high.set(500) # highest spectrum energie in eV
spectrometer_paras.undu.set(1)

"""
Setting Screen Parameter
"""
screen_paras = wave._screen_paras
screen_paras.pinh_nz.set(30) 
screen_paras.pinh_ny.set(30)
screen_paras.pinh_w.set(40) # pinhole width mm
screen_paras.pinh_h.set(40) # pinhole height mm

"""
Setting Beam Parameter
"""
ebeam_paras = wave._ebeam_paras
ebeam_paras = ebeam_paras.get_std_bessy_III_paras()

"""
Run
"""

wave.run()

"""
Get Results and Plot
"""

results = wave.get_results() # Return Result Object
nfig=0

"""
Plot Trajectory
"""
traj_x = results.get_result(which='traj_x')
traj_y = results.get_result(which='traj_y')
traj_z = results.get_result(which='traj_z')

traj_y.plot_over(x_quant=traj_x,nfig=nfig,nosave=True)
nfig=traj_z.plot_over(x_quant=traj_x,nfig=nfig,title='Trajectory')
nfig=traj_z.plot_parametric_3d(x_quant=traj_x,y_quant=traj_y,title='Trajectory',nfig=nfig)

"""
Plot B-Field
"""
By = results.get_result(which='By')
By.plot_over(x_quant=traj_x,nfig=nfig,nosave=True)
Bz = results.get_result(which='Bz')
nfig=Bz.plot_over(x_quant=traj_x,nfig=nfig,title='B-Field')

"""
Plot Power-Distribution
"""
power_z = results.get_result(which='power_z')
power_y = results.get_result(which='power_y')
power_distro = results.get_result(which='power_distribution')
nfig=power_distro.plot_over_3d(x_quant=power_y,y_quant=power_z,file_name=None,nosave=False,nfig=nfig)

"""
Plot Flux
"""
en_flux = results.get_result(which='en_flux')
flux = results.get_result(which='flux')
nfig=flux.plot_over(x_quant=en_flux,file_name=None,nosave=False,nfig=nfig,loglog=True)

"""
Plot Brilliance
"""
en_brill = results.get_result(which='en_brill')
brill0 = results.get_result(which='brill0')
brill0e = results.get_result(which='brill0e')
brill0f = results.get_result(which='brill0f')
brill0ef = results.get_result(which='brill0ef')
brill0.plot_over(x_quant=en_brill,nfig=nfig,nosave=True,loglog=True)
brill0e.plot_over(x_quant=en_brill,nfig=nfig,nosave=True,loglog=True)
brill0f.plot_over(x_quant=en_brill,nfig=nfig,nosave=True,loglog=True)
nfig=brill0ef.plot_over(x_quant=en_brill,nfig=nfig,loglog=True)

"""
Plot FLux density distributions
"""
en_fd = results.get_result(which='en_fd')
flux_density_onaxis = results.get_result(which='flux_density')
nfig=flux_density_onaxis.plot_over(x_quant=en_fd,nfig=nfig,loglog=True)

flux_dens_distr_ens_loaded = results.find_load_flux_density_distribution(energies=[460,504])
fd_y = results.get_result(which='fd_y')
fd_z = results.get_result(which='fd_z')
for en in flux_dens_distr_ens_loaded :
	fd = results.get_result(which=f'flux_density_distribution_{en:.2f}')
	nfig=fd.plot_over_3d(x_quant=fd_y,y_quant=fd_z,file_name=None,nosave=False,nfig=nfig)

pdb.set_trace()
