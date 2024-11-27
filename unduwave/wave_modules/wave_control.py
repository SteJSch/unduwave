from unduwave.unduwave_incl import *

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
			subprocess.call('"../bin/wave.exe"')        
		else:
			os.system("../bin/wave.exe")        
		if not (self.current_folder is None):
			os.chdir(self.current_folder )