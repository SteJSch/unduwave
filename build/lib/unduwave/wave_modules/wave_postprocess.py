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
		wave_folder    = self._wave_api._wave_prog_paras.wave_prog_folder.get()
		res_folder     = self._wave_api._wave_prog_paras.res_folder.get()
		res_wave       = res_folder + self._wave_api._wave_prog_paras.wave_data_res_folder.get()
		copy_behav     = self._wave_api._wave_prog_paras.wave_res_copy_behaviour.get()
		zip_res_folder = self._wave_api._wave_prog_paras.zip_res_folder.get()
		files_dont_del = []

		if copy_behav == 'copy_all' :
			wave_res_extract   = self._wave_api._wave_prog_paras.wave_ending_extract.get()
			wave_res_copy      = self._wave_api._wave_prog_paras.wave_ending_copy.get()
			wave_files_no_copy = self._wave_api._wave_prog_paras.no_copy.get()
			files_del = []
		elif copy_behav == 'copy_del_none' :
			wave_res_extract   = []
			wave_res_copy      = []
			wave_files_no_copy = []
			files_del          = []
		elif copy_behav == 'copy_essentials' :
			wave_res_extract   = self._wave_api._wave_prog_paras.wave_files_essentials.get()
			wave_res_copy      = self._wave_api._wave_prog_paras.wave_ending_copy.get()
			wave_files_no_copy = []
			files_del          = self._wave_api._wave_prog_paras.wave_ending_extract.get()
			files_dont_del     = self._wave_api._wave_prog_paras.no_copy.get()

		if not os.path.exists(res_wave):
			os.makedirs(res_wave)

		files_moved = f_h.mv_cp_files(hints = wave_res_extract, exptns = wave_files_no_copy\
				, folder_in = wave_folder + 'stage/', folder_out = res_wave, move = True, add_string = '')
		files_copied = f_h.mv_cp_files(hints = wave_res_copy, exptns = wave_files_no_copy\
				, folder_in = wave_folder + 'stage/', folder_out = res_wave, move = False, add_string = '')
		f_h.del_files(hints = files_del, exptns = files_dont_del, folder = wave_folder + 'stage/')
		summary = self.extract_summary(folder = res_wave)
		self._wave_api._summary = summary
		# if zip_res_folder : 
		#     f_h.zip_files_in_folder(folder_to_pack = res_wave)
		return files_moved + files_copied
	
	def extract_summary(self, folder=None):
		"""
		Extracts summary information from a WAVE run's files in the specified folder
		and stores the results in the file self.wave_paras.res_summary_file within the folder.

		Args:
			folder (str): The folder containing the WAVE run files.
		"""
		if folder is None :
			res_folder     = self._wave_api._wave_prog_paras.res_folder.get()
			folder       = res_folder + self._wave_api._wave_prog_paras.wave_data_res_folder.get()

		wave_out_file = []
		with open(folder+"wave.out", 'r') as o_f:
			# read an store all lines into list
			wave_out_file = o_f.readlines()

		summary = {}
		gamma = 0
		open_angle = 0
		# power, fund freq, cone opening, pinhole width at observaion point
		for ind, line in enumerate(wave_out_file) : 
			if line.find('(X,Y,Z), width, height:') >= 0 :
				line_vals = wave_out_file[ind+2]
				vals = []
				for elem in line_vals.split(' ') :
					if (len(elem) > 0) and not (elem == '\n') and not ( elem == ' ')  :
						vals.append(float(elem))
				summary.update( { 'pinhole_x [m]' : vals[0] } )
			if line.find('K, lx, N of effective periodical device:') >= 0 :
				line_vals = wave_out_file[ind+1]
				vals = []
				for elem in line_vals.split(' ') :
					if (len(elem) > 0) and not (elem == '\n') and not ( elem == ' ')  :
						vals.append(float(elem))
				summary.update( { 'Undu_Para' : vals[0] } )
			if line.find('Deflection parameter K, 1. harmonical [eV] (main poles):') >= 0 :
				line_vals = wave_out_file[ind+1]
				parts = line_vals.split(' ')
				vals=[]
				for part in parts:
					if len(part.strip()) > 0 :
						vals.append(float(part.strip()))
				harm = vals[1]
				summary.update( { 'omega0' : harm } )
			if line.find('first harmonical [eV] (estimate)') >= 0 :
				line_vals = wave_out_file[ind+1]
				parts = line_vals.split('->')
				harm = float(parts[-1])
				summary.update( { 'omega0' : harm } )
			if line.find('1. harmonical [nm] and [eV]') >= 0 :
				line_vals = wave_out_file[ind]
				parts = line_vals.split('1. harmonical [nm] and [eV], omega [1/s]: ')
				harms=parts[-1].split(' ')
				ind_f = 0
				for el in harms :
					if len(el) > 0 :
						try:
							num = float(el.strip())
							ind_f = ind_f + 1
						except:
							continue
						if ind_f == 2 :
							harm = num
							break
				summary.update( { 'omega0' : harm } )
			if line.find('Initial energy [GeV] and gamma:') >= 0 :
				part_last = line.split('Initial energy [GeV] and gamma:')[-1]
				vals = []
				for elem in part_last.split(' ') :
					if (len(elem) > 0) and not (elem == '\n') and not ( elem == ' ')  :
						vals.append(float(elem))
				gamma = vals[-1]
				open_angle = 1/gamma
			if line.find('BYmax,BYmin:') >= 0 :
				part_last = line.split('BYmax,BYmin:')[-1]
				vals = []
				for elem in part_last.split(' ') :
					if (len(elem) > 0) and not (elem == '\n') and not ( elem == ' ')  :
						vals.append(float(elem))
				summary.update( { 'bmax [T]' : vals[0] } )
				summary.update( { 'bmin [T]' : vals[1] } )
			if line.find('Power irradiated by the device [kWATT]:') >= 0 :
				part_last = line.split('Power irradiated by the device [kWATT]:')[-1]
				vals = []
				for elem in part_last.split(' ') :
					if (len(elem) > 0) and not (elem == '\n') and not ( elem == ' ')  :
						vals.append(float(elem))
				power = vals[0]
				summary.update( { 'power [kW]' : power } )
			if line.find('1. magnetic integral [T-m]:') >= 0 :
				part_last = line.split('1. magnetic integral [T-m]:')[-1]
				vals = []
				for elem in part_last.split(' ') :
					if (len(elem) > 0) and not (elem == '\n') and not ( elem == ' ')  :
						vals.append(float(elem))
				firstint = vals[1]
				summary.update( { 'first_int [Tm]' : firstint } )
			if line.find('2. mag. integral [T-m**2]:') >= 0 :
				part_last = line.split('2. mag. integral [T-m**2]:')[-1]
				vals = []
				for elem in part_last.split(' ') :
					if (len(elem) > 0) and not (elem == '\n') and not ( elem == ' ')  :
						vals.append(float(elem))
				scndint = vals[1]
				summary.update( { 'scnd_int [Tmm]' : scndint } )
		summary.update( { 'gamma' : gamma } )
		summary.update( { 'half opening angle [rad]' : open_angle } )
		if 'pinhole_x [m]' in summary :
			summary.update( { 'cone_radius_at_x [mm]' : math.tan( open_angle ) * summary['pinhole_x [m]'] * 1000 } )
		with open( folder+self._wave_api._wave_prog_paras.res_summary_file.get(), 'w') as o_f:
			for key, val in summary.items() :
				o_f.write( key + ' : ' + str(val) + '\n' )
		return summary
	
	def cleanup(self):
		"""
		Cleans up the WAVE run by removing the 'WAVE.mhb' file if it exists in the specified folder.

		Args:
			wave_folder (str): The folder containing the WAVE run files.
		"""
		wave_folder    = self._wave_api._wave_prog_paras.wave_prog_folder.get()
		wave_mbh_file = os.path.join(wave_folder, 'stage', 'WAVE.mhb')
		if os.path.exists(wave_mbh_file):
			os.remove(wave_mbh_file)
