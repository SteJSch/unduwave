"""
Undumag api definitions
"""

from unduwave.unduwave_incl import *

from unduwave.undu_modules.undu_parameters import *
from unduwave.undu_modules.undu_prepare import *
from unduwave.undu_modules.undu_control import *
from unduwave.undu_modules.undu_postprocess import *
from unduwave.undu_modules.undu_results import *

class undu_api :
	def __init__(self,undu_mode='from_clc_file') :
		"""
		Undumag API-class for controlling basic undumag-functionality. 

		:param str undu_mode: Determines the basic mode by which a magnet structure is
		constructed. Can be one of the following: \n
				'from_clc_file': The magnet structure is read from an Undumag .clc file. \n
				'from_undu_magns': A model of the magnet structure is constructed using 
				the undu_magns classes. \n
		"""
		self._prog_paras = undu_prog_parameters()
		self._prog_paras.get_std_paras(undu_mode=undu_mode)

	def run(self) :
		"""
		Runs undumag with the given settings, prepares and postprocesses data.
		"""
		prep = undu_prepare(undu_api=self) # create prepare class
		prep.create_undu_input() # creating and copying undumag input

		script_folder = os.getcwd()
		# create control class
		undu_instance = undu_control(undu_api=self,current_folder=script_folder)
		undu_instance.run() # run the simulation
		post= undu_postprocess(undu_api=self) # create class for postprocessing
		post.copy_results()
		post.cleanup()

	def load_clc_raw(self) :
		"""
		Loads a raw-clc file template from the folder set in prog_paras.in_file_folder+prog_paras.in_file_clc_raw.
		After loading prog_paras.in_file_clc_lines holds the list representation of the file in form of a list of strings.
		"""
		clc_raw=self._prog_paras.in_file_folder.get()+self._prog_paras.in_file_clc_raw.get()
		with open(clc_raw, 'r') as o_f:
			load_clc = o_f.readlines()
		self._prog_paras.in_file_clc_lines.set(load_clc)

	def set_force_calc(self,object) :
		"""
		Prepares the parameter list for calculation of force on the undu_magnets object object.

		:param undu_magnets object: An undu_magnets object (i.e. a collection of magnets) on which the force is to be calculated.
		"""
		# the objects.get_max_extent returns the minimal, rectangular volume info 
		# containing all of the objects in the object-list.
		center_coord, maxs, mins = object.get_max_extent()
		self._prog_paras.calc_force.set(1) # switching force calc on
		# setting the center coordinates of where to calc the force
		self._prog_paras.force_box_center_x.set(center_coord._x)
		self._prog_paras.force_box_center_y.set(center_coord._y)
		self._prog_paras.force_box_center_z.set(center_coord._z)
		# setting the box lengths
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
