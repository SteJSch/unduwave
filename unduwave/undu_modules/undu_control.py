"""
Functionality for running the Undumag simulation
"""

from unduwave.unduwave_incl import *
import unduwave.helpers.file_folder_helpers as f_h

try :
	# works when calling script with python3 script_file
	dir_path = os.path.dirname(os.path.realpath(__file__))
except:
	# works when calling script with exec from python console
	dir_path = os.getcwd()

class undu_control():
	def __init__(self, undu_api):
		"""
		An object for controling undumag, allowing to run simulations

		:param undu_api undu_api: external wave_api class
		:param current_folder: folder to which you want to jump back after wave was run
		:type current_folder: str or None
		"""
		self._undu_api = undu_api
		self._undu_folder = self._undu_api._prog_paras.undumag_curr_folder.get()

	def run( self ) : 
		"""
		Run an undumag simulation
		"""
		os.chdir(self._undu_folder / 'stage/' )

		exe = self._undu_folder / "bin" / ("undumag_win.exe" if os.name == "nt" else "undumag.exe")
		subprocess.run([str(exe)], check=True)

		os.chdir(ROOT_DIR)

		# if os.name == 'nt' :
		# 	os.system("../bin/undumag_win.exe")        
		# else:
		# 	os.system("../bin/undumag.exe")        

