import unittest

import unduwave as uw
from unduwave.unduwave_incl import *
from unduwave import undu_blocks

"""
Get the directory of this file
"""
try :
	# works when calling script with python3 script_file
	dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
except:
	# works when calling script with exec from python console
	dir_path = Path(os.getcwd())

wave = uw.wave() # Get a wave object
wave_prog_paras = wave._prog_paras # Get the wave program parameters
wave_prog_paras.res_folder.set(dir_path/Path('res')) # Set the result folder
wave_prog_paras.calc_spectrum.set(0)

wave.run() # Run wave

results = wave.get_results()
nfig=0

traj_x = results.get_result(which='traj_x')
traj_y = results.get_result(which='traj_y')
traj_z = results.get_result(which='traj_z')

traj_y.plot_over(x_quant=traj_x,nfig=nfig,nosave=True)
nfig=traj_z.plot_over(x_quant=traj_x,nfig=nfig)

nfig=traj_y.plot_over(x_quant=traj_z,nfig=nfig)

By = results.get_result(which='By')
By.plot_over(x_quant=traj_x,nfig=nfig,nosave=True)
Bz = results.get_result(which='Bz')
nfig=Bz.plot_over(x_quant=traj_x,nfig=nfig)

pdb.set_trace()

