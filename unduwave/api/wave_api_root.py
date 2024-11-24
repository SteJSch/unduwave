"""
The basic api
"""
from unduwave.unduwave_incl import *
from unduwave.wave_modules.wave_parameters import *
from unduwave.wave_modules.wave_prepare import *
from unduwave.wave_modules.wave_control import *
from unduwave.wave_modules.wave_postprocess import *

class wave_api :
	def __init__(self,undu_mode='undu_easy') :
		"""
		"""

		self._wave_prog_paras = wave_prog_parameters()
		self._wave_prog_paras.get_std_paras(undu_mode=undu_mode)
		self._ebeam_paras = ebeam_parameters()
		self._ebeam_paras.get_std_paras()
		self._screen_paras = screen_parameters()
		self._screen_paras.get_std_paras()
		self._spectrometer_paras = spectrometer_paras()
		self._spectrometer_paras.get_std_paras()
		self._bfield_paras = bfield_paras()
		self._bfield_paras.get_std_paras()
		self._undu_paras = undu_paras()
		self._undu_paras.get_std_paras()

	def set_bessy_II_elliptical_undu(self,nperiods) :
		"""
		Sets standard settings for bessy II and some helical undulator
		"""
		pass

	def run(self) :
		"""
		Runs wave with the given settings, prepares and postprocesses data
		"""
		
		prep = wave_prepare(wave_api=self)
		prep.create_wave_input()
		prep.prepare_b_files_for_wave()
		script_folder = os.getcwd()
		wave_instance = wave_control(wave_api=self,current_folder=script_folder)
		wave_instance.run()
		post= wave_postprocess(wave_api=self)
		post.copy_results()
		post.cleanup()
		pdb.set_trace()
		results = wave_results(post)
		return results
