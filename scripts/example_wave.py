
import pdb
import sys
sys.path.insert(0, '../')

try :
	# works when calling script with python3 script_file
	dir_path = os.path.dirname(os.path.realpath(__file__))
except:
	# works when calling script with exec from python console
	dir_path = os.getcwd()

import unduwave as uw

"""
A change 2
"""
wave = uw.wave()
wave_prog_paras = wave._wave_prog_paras
wave_prog_paras.res_folder.set(dir_path+'/res/')

wave.run()

pdb.set_trace()

