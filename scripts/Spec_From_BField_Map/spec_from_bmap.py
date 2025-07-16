
import os
import pdb
import sys
sys.path.insert(0, '../../')

try :
	# works when calling script with python3 script_file
	dir_path = os.path.dirname(os.path.realpath(__file__))
except:
	# works when calling script with exec from python console
	dir_path = os.getcwd()

import unduwave as uw

res_folder=dir_path+'/res/' # Result-Folder
field_dir=dir_path+'/field_map/' # Field-Map Folder
field_file='field_map.map' # Field Map File

wave = uw.wave(wave_mode='bmap') # Get a wave object for bmap-calculation

waveparas = wave._wave_prog_paras
waveparas.ihisascii.set(111) # extract all result-data
waveparas.field_folder.set(field_dir)# Field Folder
waveparas.field_files.set( [ field_file ] )# The magnetic field files to be used in the simulation
waveparas.pics_folder.set(f'../res/Pics/') # set the picture-result folder
waveparas.nthreads.set(6) # set the number of cpu-cores
waveparas.wave_res_copy_behaviour.set('copy_essentials') #'copy_all'
waveparas.res_folder.set(res_folder) # where to store the wave results

ebeam_paras=wave._ebeam_paras
ebeam_paras=ebeam_paras.get_std_bessy_II_paras()

"""
Setting Screen Parameter
"""
screen_paras = wave._screen_paras
screen_paras.pinh_x.set(10)
screen_paras.pinh_nz.set(15)
screen_paras.pinh_ny.set(15)
screen_paras.pinh_w.set(40)
screen_paras.pinh_h.set(40)
spectrometer_paras = wave._spectrometer_paras
spectrometer_paras.freq_num.set(100)
spectrometer_paras.freq_low.set(100)
spectrometer_paras.freq_high.set(1000)
spectrometer_paras.undu.set(1)

wave.run()

results = wave.get_results()
nfig=0

add = '' # to distinguish separate simulation results

traj_x = results.get_result(which='traj_x')
traj_y = results.get_result(which='traj_y')
traj_z = results.get_result(which='traj_z')

traj_y.plot_over(
	x_quant=traj_x,
	nfig=nfig,
	nosave=True,
	clear=True,
	dataFile=f'/{res_folder}/traj_yx{add}.txt'
	)
nfig=traj_z.plot_over(
	x_quant=traj_x,
	nfig=nfig,
	file_name=f'trajectory{add}.png',
	plot=False,
	title=f'Trajectory\n',
	dataFile=f'/{res_folder}/traj_zx{add}.txt'
	)

By = results.get_result(which='By')
By.plot_over(
	x_quant=traj_x,
	nfig=nfig,
	nosave=True,
	clear=True,
	dataFile=f'/{res_folder}/byfield{add}.txt'
	)
Bz = results.get_result(which='Bz')
nfig=Bz.plot_over(
	x_quant=traj_x,
	nfig=nfig,
	file_name=f'bfield{add}.png',
	plot=False,
	title=f'Magn. Induction',
	dataFile=f'/{res_folder}/bzfield{add}.txt'
	)

power_z = results.get_result(which='power_z')
power_y = results.get_result(which='power_y')
power_distro = results.get_result(which='power_distribution')
nfig=power_distro.plot_over_3d(
	x_quant=power_y,
	y_quant=power_z,
	file_name=f'power_distro{add}.png',
	nosave=False,
	nfig=nfig,
	plot=False,
	title=f'Power Distribution',
	clear=True,
	dataFile=f'/{res_folder}/power_distro{add}'
	)

pdb.set_trace()

