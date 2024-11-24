"""
The basic api
"""
from unduwave.unduwave_incl import *

class wave_api :
	def __init__(self) :
		"""
		Modes:
		- from B-Field, different b-fields, y/z
		- wiggler / undulator / helical undulator
		"""

		"""
		self._ebeam_paras = StdEbeamParas()
		self._screen_paras = StdScreenParas()
		self._spectrometer_paras = StdSpectrometerParas()
		self._bfield_paras = StdBFieldParas()
		self._undu_paras = StdUnduParas()
		"""

	def set_bessy_II_elliptical_undu(self,nperiods) :
		"""
		Sets standard settings for bessy II and some helical undulator
		"""
		pass

	def say_hello(self) :
		print("Hello")

	def run(self) :
		"""
		Runs wave with the given settings, prepares and postprocesses data
		"""

		prep = prepare_wave(wave_instance=self)
		prep.create_wave_in()
		prep.prepare_bfields()
		wave_instance = wave_api(prep)
		wave_instance.run(paras)
		post= postprocess_waverun(wave_instance)
		post.copy_results()
		post.clean_up()
		results = wave_results(post)
		return results

def hello() : 
	print("hello")