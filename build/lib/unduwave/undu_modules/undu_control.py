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
	def __init__(self, undu_api,current_folder = None):
		"""
		An object for controling undumag, allowing to run simulations

		:param undu_api undu_api: external wave_api class
		:param current_folder: folder to which you want to jump back after wave was run
		:type current_folder: str or None
		"""	  
		self._undu_api = undu_api
		self._current_folder = current_folder
		self._undu_folder = self._undu_api._prog_paras.undumag_prog_folder.get()

	def run( self ) : 
		"""
		Run an undumag simulation
		"""
		os.chdir(self._undu_folder + 'stage/' )
		if os.name == 'nt' :
			with open(f_h.convert_path_to_win(dir_path+'/../../External-Software/where_is_cygwin_installation.txt'), 'r') as o_f:
				cygwinfile = o_f.readlines()
			wherecyg = f_h.convert_path_to_win(cygwinfile[0].replace("'","").strip())
			whereundupy = f_h.convert_path_to_win(cygwinfile[1].replace("'","").strip())
			subprocess.call(f"{wherecyg}bin\\bash.exe --login -c 'cd {whereundupy}unduwave/External-Software/Undumag/stage; ../bin/undumag.exe'")        
		else:
			os.system("../bin/undumag.exe")        
		if not (self._current_folder is None):
			os.chdir(self._current_folder )

