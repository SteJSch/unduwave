
import os
import pdb
import sys
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
wave = uw.wave(wave_mode='undu_ellip')
wave_prog_paras = wave._prog_paras

"""
Setting Program Parameter
"""
wave_prog_paras.res_folder.set(res_folder)
wave_prog_paras.calc_spectrum.set(True)
wave_prog_paras.nthreads.set(6)

"""
Setting Spectrometer Parameter
"""

spectrometer_paras = wave._spectrometer_paras
spectrometer_paras.spectrum_n_energies.set(100)
spectrometer_paras.spectrum_min_energy.set(100)
spectrometer_paras.spectrum_max_energy.set(500)
spectrometer_paras.spectrum_undu_mode.set(1)

"""
Setting Screen Parameter
"""
screen_paras = wave._screen_paras
screen_paras.screen_segm_hor.set(30) 
screen_paras.screen_segm_vert.set(30)
screen_paras.screen_extent_hor.set(40) # pinhole width mm
screen_paras.screen_extent_vert.set(40) # pinhole height mm

"""
Setting Undulator Parameter
"""

undu_paras = wave._undu_paras # getting parameter object
undu_paras.elliptUnduB0Y.set(1.18)
undu_paras.elliptUnduB0Z.set(0.0)
undu_paras.elliptUnduNumPeriods.set(10)
undu_paras.elliptUnduPerLength.set(0.020)
undu_paras.elliptUnduPerShift.set(0.5)

"""
Setting Beam Parameter
"""
ebeam_paras = wave._ebeam_paras
ebeam_paras.get_std_bessy_II_paras()

"""
Run
"""

wave.run()

"""
Get Results and Plot
"""

results = wave.get_results()
nfig=0

traj_x = results.get_result(which='traj_x')
traj_y = results.get_result(which='traj_y')
traj_z = results.get_result(which='traj_z')

traj_y.plot_over(x_quant=traj_x,nfig=nfig,nosave=True)
nfig=traj_z.plot_over(x_quant=traj_x,nfig=nfig)

By = results.get_result(which='By')
By.plot_over(x_quant=traj_x,nfig=nfig,nosave=True)
Bz = results.get_result(which='Bz')
nfig=Bz.plot_over(x_quant=traj_x,nfig=nfig)

power_z = results.get_result(which='power_z')
power_y = results.get_result(which='power_y')
power_distro = results.get_result(which='power_distribution')
nfig=power_distro.plot_over_3d(x_quant=power_y,y_quant=power_z,file_name=None,nosave=False,nfig=nfig)

en_flux = results.get_result(which='en_flux')
flux = results.get_result(which='flux')
nfig=flux.plot_over(x_quant=en_flux,file_name=None,nosave=False,nfig=nfig,loglog=True)

en_brill = results.get_result(which='en_brill')
brill0 = results.get_result(which='brill0')
brill0e = results.get_result(which='brill0e')
brill0f = results.get_result(which='brill0f')
brill0ef = results.get_result(which='brill0ef')
brill0.plot_over(x_quant=en_brill,nfig=nfig,nosave=True,loglog=True)
brill0e.plot_over(x_quant=en_brill,nfig=nfig,nosave=True,loglog=True)
brill0f.plot_over(x_quant=en_brill,nfig=nfig,nosave=True,loglog=True)
nfig=brill0ef.plot_over(x_quant=en_brill,nfig=nfig,loglog=True)

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
