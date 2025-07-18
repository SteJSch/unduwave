"""
Wave api definitions
"""
from unduwave.unduwave_incl import *
from unduwave.wave_modules.wave_parameters import *
from unduwave.wave_modules.wave_prepare import *
from unduwave.wave_modules.wave_control import *
from unduwave.wave_modules.wave_postprocess import *
from unduwave.wave_modules.wave_results import *

class wave_api :

	_prog_paras = wave_prog_parameters()
	_ebeam_paras = ebeam_parameters()
	_screen_paras = screen_parameters()
	_spectrometer_paras = spectrometer_paras()
	_bfield_paras = bfield_paras()
	_undu_paras = undu_paras()

	def __init__(self,wave_mode='undu_endp') :
		"""
		Wave API-class for controlling basic wave-functionality. 
		Holds the basic parameter classes.

		:param str wave_mode: can be one of the following.
			'By','Byz',	'Bxyz',	'undu_ellip', 'undu_endp' :
		"""

		self._prog_paras.get_std_paras(wave_mode=wave_mode)
		self._ebeam_paras.get_std_paras()
		self._screen_paras.get_std_paras()
		self._spectrometer_paras.get_std_paras()
		self._bfield_paras.get_std_paras()
		self._undu_paras.get_std_paras(wave_mode=wave_mode,ebeam=self._ebeam_paras)

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

	def get_results(self) :
		"""
		Returns the results from a given simulation as result-object.
		"""
		results = wave_results(wave_api=self)
		results.load_from_res_folder()
		return results
