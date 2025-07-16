from unduwave.unduwave_incl import *
from unduwave.wave_modules.wave_parameters import *
import unduwave.helpers.file_folder_helpers as f_h

class undu_postprocess:
	"""
	A class for postprocessing Undumag files.
	
	Args:
		undu_api 
	"""
	
	def __init__(self, undu_api):
		"""
		Takes external api and create postprocess class
		"""
		self._undu_api = undu_api

	def copy_results(self):
		"""
		Cleans the undumag-stage folder and copies the desired files to their location,
		deletes non-desired files, and zips the results based on the undu_res_copy_behaviour setting.
		"""
		#'copy_all', 'copy_none' - only writes res_summary, 'copy_essentials' 
		undu_folder    = self._undu_api._prog_paras.undumag_prog_folder.get()
		res_folder     = self._undu_api._prog_paras.res_folder.get()
		res_undu       = res_folder + self._undu_api._prog_paras.undu_data_res_folder.get()
		copy_behav     = self._undu_api._prog_paras.undu_res_copy_behaviour.get()
		zip_res_folder = self._undu_api._prog_paras.zip_res_folder.get()
		files_dont_del = []

		if copy_behav == 'copy_all' :
			undu_res_extract   = self._undu_api._prog_paras.undu_ending_extract.get()
			undu_res_copy      = self._undu_api._prog_paras.undu_ending_copy.get()
			undu_files_no_copy = self._undu_api._prog_paras.no_copy.get()
			files_del = []
		elif copy_behav == 'copy_del_none' :
			undu_res_extract   = []
			undu_res_copy      = []
			undu_files_no_copy = []
			files_del          = []
		elif copy_behav == 'copy_essentials' :
			undu_res_extract   = self._undu_api._prog_paras.undu_files_essentials.get()
			undu_res_copy      = self._undu_api._prog_paras.undu_ending_copy.get()
			undu_files_no_copy = []
			files_del          = self._undu_api._prog_paras.undu_ending_extract.get()
			files_dont_del     = self._undu_api._prog_paras.no_copy.get()

		if not os.path.exists(res_undu):
			os.makedirs(res_undu)
		files_moved = f_h.mv_cp_files(hints = undu_res_extract, exptns = undu_files_no_copy\
				, folder_in = undu_folder + 'stage/', folder_out = res_undu, move = True, add_string = '')
		files_copied = f_h.mv_cp_files(hints = undu_res_copy, exptns = undu_files_no_copy\
				, folder_in = undu_folder + 'stage/', folder_out = res_undu, move = False, add_string = '')
		f_h.del_files(hints = files_del, exptns = files_dont_del, folder = undu_folder + 'stage/')
		# if zip_res_folder : 
		#     f_h.zip_files_in_folder(folder_to_pack = res_undu)
		return files_moved + files_copied
	
	def cleanup(self):
		"""
		Cleans up the Undumag run

		"""
		pass
