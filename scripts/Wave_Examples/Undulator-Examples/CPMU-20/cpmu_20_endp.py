import unduwave as uw
from unduwave.unduwave_incl import *
from unduwave import undulatorComponents
from unduwave import undu_blocks
from unduwave import magneticObjectGeometries

try :
	# works when calling script with python3 script_file
	dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
except:
	# works when calling script with exec from python console
	dir_path = Path(os.getcwd())

res_folder = dir_path/Path('res/')

"""
Getting wave
"""
wave = uw.wave(wave_mode='undu_endp') # Simple undulator model with endpoles
wave_prog_paras = wave._prog_paras

"""
Setting Undulator Parameter

We give the K-parameters and the B-Amplitudes are calculated from there
"""

K_para_y=0.5

period_length = 0.02 # m
nperiods = 3

undu_paras = wave._undu_paras # getting parameter object
undu_paras.unduParameterKY.set(K_para_y)
undu_paras.numPeriods.set(nperiods) # number full periods
undu_paras.periodLength.set(period_length)

"""
Setting Program Parameter
"""
wave_prog_paras.res_folder.set(res_folder)
wave_prog_paras.calc_spectrum.set(False)
wave_prog_paras.nthreads.set(6)

"""
Setting Spectrometer Parameter
"""

spectrometer_paras = wave._spectrometer_paras
spectrometer_paras.spectrum_n_energies.set(10)
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
Setting Beam Parameter
"""
ebeam_paras = wave._ebeam_paras
ebeam_paras = ebeam_paras.get_std_bessy_III_paras()

"""
Run
"""
print("run")
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

pdb.set_trace()
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
nfig=flux.plot_over(x_quant=en_flux,nfig=nfig,loglog=True)

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
