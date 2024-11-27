from unduwave.unduwave_incl import *

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
			pdb.set_trace()
			with open(folder+"wave.out", 'r') as o_f:
				# read an store all lines into list
				wave_out_file = o_f.readlines()
			subprocess.call("C:\\cygwin64\\bin\\bash.exe --login -c 'cd unduwave/External-Software/WAVE/stage; ../bin/wave.exe'")        
		else:
			os.system("../bin/wave.exe")        
		if not (self.current_folder is None):
			os.chdir(self.current_folder )