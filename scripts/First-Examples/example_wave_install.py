"""
Easy example script that just runs wave
"""
import pdb
import sys
import os

"""
Get the directory of this file
"""
try :
	# works when calling script with python3 script_file
	dir_path = os.path.dirname(os.path.realpath(__file__))
except:
	# works when calling script with exec from python console
	dir_path = os.getcwd()

import unduwave as uw # Import unduwave

wave = uw.wave() # Get a wave object
wave_prog_paras = wave._wave_prog_paras # Get the wave program parameters
wave_prog_paras.res_folder.set(dir_path+'/res/') # Set the result folder

wave.run() # Run wave

pdb.set_trace()

