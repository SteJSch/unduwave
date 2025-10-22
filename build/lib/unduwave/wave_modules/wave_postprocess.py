from unduwave.unduwave_incl import *
from unduwave.wave_modules.wave_parameters import *
import unduwave.helpers.file_folder_helpers as f_h

class wave_postprocess:
	"""
	A class for postprocessing WAVE files.
	
	Args:
		wave_api (StdParameters): An instance of the StdParameters class.
	"""
	
	def __init__(self, wave_api):
		"""
		Takes external api and create postprocess class
		"""
		self._wave_api = wave_api

	def copy_results(self):
		"""
		Cleans the wave-stage folder and copies the desired files to their location,
		deletes non-desired files, and zips the results based on the wave_res_copy_behaviour setting.
		"""
		#'copy_all', 'copy_none' - only writes res_summary, 'copy_essentials' 
		wave_folder    = self._wave_api._prog_paras.wave_prog_folder.get()
		res_folder     = self._wave_api._prog_paras.res_folder.get()
		res_wave       = res_folder + self._wave_api._prog_paras.wave_data_res_folder.get()
		copy_behav     = self._wave_api._prog_paras.wave_res_copy_behaviour.get()
		zip_res_folder = self._wave_api._prog_paras.zip_res_folder.get()
		files_dont_del = []

		if copy_behav == 'copy_all' :
			wave_res_extract   = self._wave_api._prog_paras.wave_ending_extract.get()
			wave_res_copy      = self._wave_api._prog_paras.wave_ending_copy.get()
			wave_files_no_copy = self._wave_api._prog_paras.no_copy.get()
			files_del = []
		elif copy_behav == 'copy_del_none' :
			wave_res_extract   = []
			wave_res_copy      = []
			wave_files_no_copy = []
			files_del          = []
		elif copy_behav == 'copy_essentials' :
			wave_res_extract   = self._wave_api._prog_paras.wave_files_essentials.get()
			wave_res_copy      = self._wave_api._prog_paras.wave_ending_copy.get()
			wave_files_no_copy = []
			files_del          = self._wave_api._prog_paras.wave_ending_extract.get()
			files_dont_del     = self._wave_api._prog_paras.no_copy.get()

		if not os.path.exists(res_wave):
			os.makedirs(res_wave)

		files_moved = f_h.mv_cp_files(hints = wave_res_extract, exptns = wave_files_no_copy\
				, folder_in = wave_folder + 'stage/', folder_out = res_wave, move = True, add_string = '')
		files_copied = f_h.mv_cp_files(hints = wave_res_copy, exptns = wave_files_no_copy\
				, folder_in = wave_folder + 'stage/', folder_out = res_wave, move = False, add_string = '')
		f_h.del_files(hints = files_del, exptns = files_dont_del, folder = wave_folder + 'stage/')
		# if zip_res_folder : 
		#     f_h.zip_files_in_folder(folder_to_pack = res_wave)
		return files_moved + files_copied
	
	def cleanup(self):
		"""
		Cleans up the WAVE run by removing the 'WAVE.mhb' file if it exists in the specified folder.

		Args:
			wave_folder (str): The folder containing the WAVE run files.
		"""
		wave_folder    = self._wave_api._prog_paras.wave_prog_folder.get()
		wave_mbh_file = os.path.join(wave_folder, 'stage', 'WAVE.mhb')
		if os.path.exists(wave_mbh_file):
			os.remove(wave_mbh_file)
