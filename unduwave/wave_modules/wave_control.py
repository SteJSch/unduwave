from unduwave.unduwave_incl import *
import unduwave.helpers.file_folder_helpers as f_h

try :
	# works when calling script with python3 script_file
	dir_path = os.path.dirname(os.path.realpath(__file__))
except:
	# works when calling script with exec from python console
	dir_path = os.getcwd()

class wave_control():
	"""
	Internal API for the WAVE program
	"""  
	def __init__(self, wave_api,current_folder = None):
		"""
		Initialize the internal API
		wave_api : external wave_api class
		current_folder : folder to which you want to jump back after wave was run
		"""	  
		self._wave_api = wave_api
		self.wave_folder = self._wave_api._prog_paras.wave_prog_folder()

	def run(self):
		"""Run Wave from the self.wave_folder.

		If given, change the directory back to self.current_folder
		"""
		"""
		update paras
		"""

		wave_folder = Path(self._wave_api._prog_paras.wave_curr_folder())
		self._wave_api._undu_paras.update_values()
		self._wave_api._ebeam_paras.update_values()
		os.chdir(wave_folder/'stage/' )
		if os.name == 'nt' :
			os.system("../bin/wave_win.exe")
		else:
			os.system("../bin/wave.exe")        
		os.chdir(ROOT_DIR)
