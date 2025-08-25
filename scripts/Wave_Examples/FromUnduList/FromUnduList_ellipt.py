import pandas as pd
import os
import pdb
import sys

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

"""
Setting Program Parameter
"""
wave_prog_paras = wave._prog_paras
wave_prog_paras.res_folder.set(res_folder)
wave_prog_paras.calc_spectrum.set(True)
wave_prog_paras.nthreads.set(6)

"""
Setting Beam Parameter
"""
ebeam_paras = wave._ebeam_paras
ebeam_paras.get_std_bessy_II_paras()

bessyIIlist=uw.loadBessyIIundulatorList()

bfieldEnds=uw.bfield.bfield(
		unitsXB=[0.001,1.0],
	)

bfieldEnds.create_harm_field(
		period_length=0.03,
		amplitude=1,
		numPer=20, 
		phase_shift = 0, 
		num_pnts_per_period = 100, 
		colx = 'x', 
		coly = 'By',
		) 

bfieldEnds=uw.bfield.create_field_with_ends(
		nperiods=20,
		num_pnts_per_period=100,
		colx = 'x', 
		coly = 'By',
		unitsXB=[1.0,1.0],
		shift=0,
		)

nfig=0
bfieldEnds.by.plot_over(x_quant=bfieldEnds.xvals,nfig=nfig,nosave=True)

pdb.set_trace()
"""
Setting Undulator Parameter
"""
undu_paras = wave._undu_paras # getting parameter object
undu_paras.set_paras_from_bessyII_undu_list(unduName='UE112')
undu_paras.bEffY.set(1.0)
undu_paras.bEffZ.set(1.0)
"""
Horizontal mode. 
The bessyIIlist['beffy [T]'] value is used only.
"""
undu_paras.shift.set(0.0)
"""
Inclined Mode (anti-parallel shift)
The bessyIIlist['beffy [T]'] value is used for both,
the effective bfield in y- and z-direction
"""
undu_paras.shift.set(-0.5)
"""
Circular mode (right-handed). 
The bessyIIlist['beffy [T]'] and bessyIIlist['beffz [T]'] values are used.
"""
undu_paras.shift.set(0.25)
"""
Vertical mode. 
The bessyIIlist['beffz [T]'] value is used only.
"""
undu_paras.shift.set(0.5)
"""
Circular mode (left-handed). 
The bessyIIlist['beffy [T]'] and bessyIIlist['beffz [T]'] values are used.
"""
undu_paras.shift.set(0.75)

"""
Setting Spectrometer Parameter
"""

spectrometer_paras = wave._spectrometer_paras
spectrometer_paras.spectrum_n_energies.set(200)
spectrometer_paras.spectrum_min_energy.set(1.0)
spectrometer_paras.spectrum_max_energy.set(100.0)
spectrometer_paras.spectrum_undu_mode.set(1)

"""
Setting Screen Parameter
"""
screen_paras = wave._screen_paras
screen_paras.screen_segm_hor.set(20) 
screen_paras.screen_segm_vert.set(20)
screen_paras.screen_extent_hor.set(40) # pinhole width mm
screen_paras.screen_extent_vert.set(40) # pinhole height mm

"""
Run
"""

wave.run()

"""
Get Results and Plot
"""

results = wave.get_results()

summary=results._summary
print(f"The summary:\n{pd.DataFrame([summary])}")

nfig=0

traj_x = results.get_result(which='traj_x')
traj_y = results.get_result(which='traj_y')
traj_z = results.get_result(which='traj_z')

traj_y.plot_over(x_quant=traj_x,nfig=nfig,nosave=True)
nfig=traj_z.plot_over(x_quant=traj_x,nfig=nfig)

nfig=traj_y.plot_over(
	x_quant=traj_z,
	nfig=nfig,
	file_name=f'trajectory_projection.png',
	plot=False,
	clear=True,
	title=f'Trajectory Projection',
	# xlim=[-1.5e-5,1.5e-5],
	# ylim=[-1.5e-5,1.5e-5],
	)

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

flux_dens_distr_ens_loaded = results.find_load_flux_density_distribution(energies=[7])
fd_y = results.get_result(which='fd_y')
fd_z = results.get_result(which='fd_z')
for en in flux_dens_distr_ens_loaded :
	fd = results.get_result(which=f'flux_density_distribution_{en:.2f}')
	nfig=fd.plot_over_3d(
		x_quant=fd_y,
		y_quant=fd_z,
		nosave=False,
		nfig=nfig,
		plot=False,
		file_name=None,
		title=f'Flux Density - Hirst Harmonic\n',
		zLabelAdd='\n\n',
		clear=True,
		)

	s1 = results.get_result(which=f's1_distribution_{en:.2f}')
	nfig=s1.plot_over_3d(
		x_quant=fd_y,
		y_quant=fd_z,
		nosave=False,
		nfig=nfig,
		plot=False,
		title=f'S1 - Hirst Harmonic',
		zLabelAdd='\n\n',
		clear=True,
		)

	s2 = results.get_result(which=f's2_distribution_{en:.2f}')
	nfig=s2.plot_over_3d(
		x_quant=fd_y,
		y_quant=fd_z,
		nosave=False,
		nfig=nfig,
		plot=False,
		title=f'S2 - Hirst Harmonic',
		zLabelAdd='\n\n',
		clear=True,
		)

	s3 = results.get_result(which=f's3_distribution_{en:.2f}')
	nfig=s3.plot_over_3d(
		x_quant=fd_y,
		y_quant=fd_z,
		nosave=False,
		nfig=nfig,
		plot=False,
		title=f'S3 - Hirst Harmonic',
		zLabelAdd='\n\n',
		clear=True,
		)

s0 = results.get_result(which='stokes_para_0')
en_s0 = results.get_result(which='en_s0')
nfig = s0.plot_over(
	x_quant=en_s0,
	nfig=nfig,
	nosave=False,
	clear=True,
	)
s1 = results.get_result(which='stokes_para_1')
en_s1 = results.get_result(which='en_s1')
nfig = s1.plot_over(
	x_quant=en_s1,
	nfig=nfig,
	nosave=False,
	clear=True,
	)
s2 = results.get_result(which='stokes_para_2')
en_s2 = results.get_result(which='en_s2')
nfig = s2.plot_over(
	x_quant=en_s2,
	nfig=nfig,
	nosave=False,
	clear=True,
	)
s3 = results.get_result(which='stokes_para_3')
en_s3 = results.get_result(which='en_s3')
nfig = s3.plot_over(
	x_quant=en_s3,
	nfig=nfig,
	nosave=False,
	clear=True,
	)

pdb.set_trace()
