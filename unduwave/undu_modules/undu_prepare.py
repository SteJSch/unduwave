from unduwave.unduwave_incl import *
from unduwave.undu_modules.undu_parameters import *

class undu_prepare():
	"""_summary_

		Args:
			undu_api (undu_api): Standard Parameters used for the simulation
	""" 
	def __init__(self, undu_api):
		"""
		Ini wave-prepare class for setting everything up for calcs
		"""
		self._undu_api = undu_api	
		self._all_paras = [ 
			self._undu_api._prog_paras 
			]

	def create_fresh_clc(self) :
		"""
		Writes/Wrote the clc file from some parameter list - but this is deprecated for now
		"""
		undu_paras = self._undu_api._prog_paras
		undu_folder= undu_paras.undumag_prog_folder.get()
		if undu_paras.undu_mode.get() == 'from_clc_file' : 
			if len(undu_paras.copy_clc_folder.get()) > 0 : 
				clcFile=undu_paras.copy_clc_folder.get()+undu_paras.in_file_clc.get()
				os.system( 'cp ' + clcFile + ' ' + undu_folder + 'stage/undumag.clc' )
		elif undu_paras.undu_mode.get() == 'from_undu_magns' :
			clc_lines=undu_paras.in_file_clc_lines.get()
			with open( undu_folder + 'stage/undumag.clc', 'w') as o_f:
				for ind, line in enumerate(clc_lines) :
					if line.find('$PerLen = ')>=0 :
						line=f'$PerLen = {undu_paras.periodLength.get()*1e3:.4f}'
					o_f.write(line)

	def create_fresh_nam( self ) :
		"""
		Updates the nam file from parameters.

		Loads the input file set in undu_paras, updates properties 
		based on other undu_paras properties,and copies the 
		resulting file to the Undu program folder.
		"""
		undu_paras = self._undu_api._prog_paras

		undu_folder   = self._undu_api._prog_paras.undumag_prog_folder.get()
		inp_folder    = self._undu_api._prog_paras.in_file_folder.get()
		configFile_in = self._undu_api._prog_paras.in_file_nam.get()
		# open the configuration file
		undu_in_file = []
		with open(inp_folder+configFile_in, 'r') as o_f:
			# read an store all lines into list
			undu_in_file = o_f.readlines()

		for ind, line in enumerate(undu_in_file) :

			if line.find('=') >= 0 :
				for para_list in self._all_paras :					
					for para in para_list.children():
						if not (para.get_in_name() is None) : 
							if line.find(f' {para.get_in_name()}=') >= 0 :
								stuff1 = line.split(f'{para.get_in_name()}=')
								stuff2 = stuff1[-1].split('!')
								try:
									undu_in_file[ind] = stuff1[0]+f'{para.get_in_name()}='+f'{para.get()*para.get_fac()}'+' !' + stuff2[-1]
								except:
									pdb.set_trace()
		with open( undu_folder + 'stage/undumag.nam', 'w') as o_f:
			for ind, line in enumerate(undu_in_file) :
				o_f.write(line)

	def create_undu_input(self):
		"""
		Creates all the files needed as input for Undumag.
		"""

		self.create_fresh_nam()
		self.create_fresh_clc()

