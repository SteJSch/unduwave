import unittest

import unduwave as uw
from unduwave.unduwave_incl import *
from unduwave import undu_blocks

try :
	# works when calling script with python3 script_file
	dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
except:
	# works when calling script with exec from python console
	dir_path = Path(os.getcwd())

class wave_settings_test(unittest.TestCase) :

	def test_wave_settings(self) :

		self.res_folder = dir_path/'res/'

		"""
		Getting wave
		"""
		wave = uw.wave(wave_mode='undu_ellip')
		wave_prog_paras = wave._prog_paras

		"""
		Setting Program Parameter
		"""
		wave_prog_paras.res_folder.set(self.res_folder)
		wave_prog_paras.calc_spectrum.set(True)
		wave_prog_paras.nthreads.set(6)

		"""
		Setting Spectrometer Parameter
		"""

		spectrometer_paras = wave._spectrometer_paras
		spectrometer_paras.spectrum_n_energies.set(10)
		spectrometer_paras.spectrum_min_energy.set(400)
		spectrometer_paras.spectrum_max_energy.set(500)
		spectrometer_paras.spectrum_undu_mode.set(1)

		"""
		Setting Screen Parameter
		"""
		screen_paras = wave._screen_paras
		screen_paras.screen_segm_hor.set(11) 
		screen_paras.screen_segm_vert.set(11)
		screen_paras.screen_extent_hor.set(40) # pinhole width mm
		screen_paras.screen_extent_vert.set(40) # pinhole height mm

		"""
		Setting Undulator Parameter
		"""

		undu_paras = wave._undu_paras # getting parameter object
		undu_paras.bEffY.set(1.18)
		undu_paras.bEffZ.set(0.5)
		undu_paras.elliptUnduNumPeriods.set(10)
		undu_paras.elliptUnduPerLength.set(0.020)
		undu_paras.elliptUnduPerShift.set(0.0)

		"""
		Setting Beam Parameter
		"""
		ebeam_paras = wave._ebeam_paras
		ebeam_paras.get_std_bessy_II_paras()

		"""
		Run
		"""

		wave.run()

	def setUp(self):
		self.resource = "Resource allocated"
		print("Setting up test resources...")

	def tearDown(self):
		print("Cleaning up resource...")

		ROOT_DIR_TEST = Path(__file__).resolve().parent
		del_dirs=[
			ROOT_DIR_TEST/'__pycache__',
			ROOT_DIR_TEST/self.res_folder,
			]

		for del_dir in del_dirs :
			if os.path.exists(del_dir) :		
				shutil.rmtree(del_dir)

if __name__ == '__main__':
	unittest.main()
