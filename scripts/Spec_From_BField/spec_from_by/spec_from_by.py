
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
traj_x = results.get_result(which='traj_x')
By = results.get_result(which='By')
By.plot_over(x_quant=traj_x)

power_z = results.get_result(which='power_z')
power_y = results.get_result(which='power_y')
power_distro = results.get_result(which='power_distribution')
power_distro.plot_over_3d(x_quant=power_y,y_quant=power_z,file_name=None,nosave=False)

pdb.set_trace()
