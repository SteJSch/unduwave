"""
Undumag api definitions
"""
from unduwave.unduwave_incl import *

# import unduwave.undu_parameters as undu_parameters

from unduwave.undu_modules.undu_parameters import *
from unduwave.undu_modules.undu_prepare import *
from unduwave.undu_modules.undu_control import *
from unduwave.undu_modules.undu_postprocess import *
from unduwave.undu_modules.undu_results import *

class undu_api :
	"""
	Undumag API-class for controlling basic undumag-functionality. 
	Holds the basic parameter classes.
	"""
	def __init__(self,undu_mode='from_clc_file') :
		"""
		Initialize the Undumag parameters

		:param str undu_mode: can be one of the following: 'from_clc_file', 'from_undu_magns'
		"""
		self._prog_paras = undu_prog_parameters()
		self._prog_paras.get_std_paras(undu_mode=undu_mode)

	def run(self) :
		"""
		Runs undumag with the given settings, prepares and postprocesses data
		"""
		prep = undu_prepare(undu_api=self)
		prep.create_undu_input()

		script_folder = os.getcwd()
		undu_instance = undu_control(undu_api=self,current_folder=script_folder)
		undu_instance.run()
		post= undu_postprocess(undu_api=self)
		post.copy_results()
		post.cleanup()

	def load_clc_raw(self) :
		clc_raw=self._prog_paras.in_file_folder.get()+self._prog_paras.in_file_clc_raw.get()
		with open(clc_raw, 'r') as o_f:
			load_clc = o_f.readlines()
		self._prog_paras.in_file_clc_lines.set(load_clc)

	def set_force_calc(self,object) :
		center_coord, maxs, mins = object.get_max_extent()
		self._prog_paras.calc_force.set(1)
		self._prog_paras.force_box_center_x.set(center_coord._x)
		self._prog_paras.force_box_center_y.set(center_coord._y)
		self._prog_paras.force_box_center_z.set(center_coord._z)
		self._prog_paras.force_box_len_x.set(maxs[0]-mins[0])
		self._prog_paras.force_box_len_y.set(maxs[1]-mins[1])
		self._prog_paras.force_box_len_z.set(maxs[2]-mins[2])

	def get_results(self) :
		"""
		Returns the results from a given simulation as result-object.
		"""
		results = undu_results(undu_api=self)
		results.load_from_res_folder()
		return results
