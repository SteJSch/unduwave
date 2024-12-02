
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

field_folder = f'/'
res_folder = dir_path+'/res/'

"""
Getting wave
"""

wave = uw.wave(undu_mode='undu_ellip')
wave_prog_paras = wave._wave_prog_paras
ebeam_paras = wave._ebeam_paras

"""
Setting Program Parameter
"""

wave_prog_paras.res_folder.set(res_folder)
wave_prog_paras.spec_calc.set(True)

spectrometer_paras = wave._spectrometer_paras
spectrometer_paras.freq_num.set(10)

"""
Setting Undulator Parameter
"""

undu_paras = wave._undu_paras
undu_paras.b0y.set(1.5)
undu_paras.b0z.set(1.0)
undu_paras.nper.set(10)
undu_paras.perl_x.set(0.051)
undu_paras.ell_shift.set(0.25)

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
Bz = results.get_result(which='Bz')
Bz.plot_over(x_quant=traj_x)

power_z = results.get_result(which='power_z')
power_y = results.get_result(which='power_y')
power_distro = results.get_result(which='power_distribution')
power_distro.plot_over_3d(x_quant=power_y,y_quant=power_z,file_name=None,nosave=False)

pdb.set_trace()
