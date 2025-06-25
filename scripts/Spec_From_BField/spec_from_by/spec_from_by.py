
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

wave = uw.wave(undu_mode='By')
wave_prog_paras = wave._wave_prog_paras
ebeam_paras = wave._ebeam_paras

"""
Setting Program Parameter
"""

wave_prog_paras.res_folder.set(res_folder)
wave_prog_paras.field_files.set( [ 'field.dat' ] )# The magnetic field files to be used in the simulation
wave_prog_paras.field_folder.set(dir_path+field_folder)# Field Folder
wave_prog_paras.spec_calc.set(True)
wave_prog_paras.nthreads.set(6)

"""
Setting Spectrometer Parameter
"""

spectrometer_paras = wave._spectrometer_paras
spectrometer_paras.freq_num.set(100)
spectrometer_paras.freq_low.set(100)
spectrometer_paras.freq_high.set(500)
spectrometer_paras.undu.set(0)

"""
Setting Screen Parameter
"""
screen_paras = wave._screen_paras
screen_paras.pinh_nz.set(30)
screen_paras.pinh_ny.set(30)
screen_paras.pinh_w.set(40)
screen_paras.pinh_h.set(40)

"""
Setting Undulator Parameter
"""

undu_paras = wave._undu_paras # getting parameter object
undu_paras.b0y.set(1.18)
undu_paras.b0z.set(0.0)
undu_paras.nper.set(10)
undu_paras.perl_x.set(0.051)
undu_paras.ell_shift.set(0.5)

"""
Setting Beam Parameter
"""

ebeam_paras.beam_en.set(1.722) # [GeV]
ebeam_paras.current.set(0.3) # [A]
ebeam_paras.bsigz.set(275e-6) # 
ebeam_paras.bsigzp.set(28.1e-6) #
ebeam_paras.bsigy.set(22.5e-6) # 
ebeam_paras.bsigyp.set(6.8e-6) # 
ebeam_paras.espread.set(1e-3) # 
ebeam_paras.emitt_h.set(7.7e-9)
ebeam_paras.emitt_v.set(15.4e-11)

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

flux_dens_distr_ens_loaded = results.find_load_flux_density_distribution(energies=[414])
fd_y = results.get_result(which='fd_y')
fd_z = results.get_result(which='fd_z')
for en in flux_dens_distr_ens_loaded :
	fd = results.get_result(which=f'flux_density_distribution_{en:.2f}')
	nfig=fd.plot_over_3d(x_quant=fd_y,y_quant=fd_z,file_name=None,nosave=False,nfig=nfig)

pdb.set_trace()
