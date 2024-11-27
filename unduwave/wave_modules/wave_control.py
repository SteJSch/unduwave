from unduwave.unduwave_incl import *
import unduwave.helpers.file_folder_helpers as f_h

try :
	# works when calling script with python3 script_file
	dir_path = os.path.dirname(os.path.realpath(__file__))
except:
	# works when calling script with exec from python console
	dir_path = os.getcwd()

class wave_control():
	"""API for the WAVE program

	Args:
		current_folder (str): _description_
		wave_folder (str): _description_
	"""  
	def __init__(self, wave_api,current_folder = None):
	  
		self._wave_api = wave_api
		self.current_folder = current_folder
		self.wave_folder    = self._wave_api._wave_prog_paras.wave_prog_folder.get()

	def run(self):
		"""Run Wave from the self.wave_folder.

		If given, change the directory back to self.current_folder
		"""        
		os.chdir(self.wave_folder + 'stage/' )
		if os.name == 'nt' :
			with open(f_h.convert_path_to_win(dir_path+'/../../External-Software/where_is_cygwin_installation.txt'), 'r') as o_f:
				cygwinfile = f_h.convert_path_to_win(o_f.readlines())
			wherecyg = cygwinfile[0].replace("'","").strip()
			whereundupy = f_h.convert_path_to_win(cygwinfile[1].replace("'","").strip())
			pdb.set_trace()
			subprocess.call(f"{wherecyg}bin\\bash.exe --login -c 'cd {whereundupy}unduwave/External-Software/WAVE/stage; ../bin/wave.exe'")        
		else:
			os.system("../bin/wave.exe")        
		if not (self.current_folder is None):
			os.chdir(self.current_folder )