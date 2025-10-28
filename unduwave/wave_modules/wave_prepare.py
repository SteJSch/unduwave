from unduwave.unduwave_incl import *
from unduwave.wave_modules.wave_parameters import *
import unduwave.analytical_module.ana_undu.bfield as bfield
import unduwave.helpers.file_folder_helpers as f_h

class wave_prepare():
	"""_summary_

		Args:
			wave_api (wave_api): Standard Parameters used for the simulation
	""" 
	def __init__(self, wave_api):
		"""
		Ini wave-prepare class for setting everything up for calcs
		"""
		self._wave_api = wave_api	
		self._all_paras = [ 
			self._wave_api._spectrometer_paras, self._wave_api._ebeam_paras, 
			self._wave_api._screen_paras, self._wave_api._bfield_paras, 
			self._wave_api._undu_paras, self._wave_api._prog_paras ]
		self._wave_api._prog_paras.update_values(
			bfield_paras=self._wave_api._bfield_paras
			)
		self._wave_api._ebeam_paras.update_values()
		self._wave_api._undu_paras.update_values()

	def create_wave_input(self):
		"""Creates all the files needed as input fro Wave.
		
		Loads the input file set in wave_paras, updates properties 
		based on other wave_paras properties,and copies the 
		resulting file to the WAVE program folder.
		"""

		wave_folder   = self._wave_api._prog_paras.wave_prog_folder.get()
		inp_folder    = self._wave_api._prog_paras.in_file_folder.get()
		b_type        = self._wave_api._prog_paras.b_type.get()
		configFile_in = self._wave_api._prog_paras.in_files.get()[b_type]
		# open the configuration file
		wave_in_file = []
		with open(inp_folder+configFile_in, 'r') as o_f:
			# read an store all lines into list
			wave_in_file = o_f.readlines()

		for ind, line in enumerate(wave_in_file) :

			if line.find('=') >= 0 :
				for para_list in self._all_paras :					
					for para in para_list.children():
						if not (para.get_in_name() is None) : 
							if line.find(f' {para.get_in_name()}=') >= 0 :
								stuff1 = line.split(f'{para.get_in_name()}=')
								stuff2 = stuff1[-1].split('!')
								wave_in_file[ind] = stuff1[0]+f'{para.get_in_name()}='+f'{para.get()*para.get_fac()}'+' !' + stuff2[-1]
		with open( wave_folder + 'stage/wave.in', 'w') as o_f:
			for ind, line in enumerate(wave_in_file) :
				o_f.write(line)

	def prepare_b_files_for_wave(self):
		"""Prepare the files for WAVE depending on the b type.
			
		Depending on which b_type, copies and 
		formats the b-field files needed
		"""
		wave_paras=self._wave_api._prog_paras
		wave_mode = wave_paras.wave_mode.get()
		wave_folder = wave_paras.wave_prog_folder.get()
		# field_files = self._wave_api._prog_paras.field_files.get()
		# field_folder = self._wave_api._prog_paras.field_folder.get()

		if wave_mode == 'four' :
			fourIn=wave_paras.four_file.get()
			fourInFile = os.path.basename(fourIn)
			fourInDirec=os.path.dirname(fourIn)+os.sep

			f_h.mv_cp_files(
				hints=[fourInFile], 
				exptns=[], 
				folder_in=fourInDirec, 
				folder_out=wave_folder + 'stage'+os.sep, 
				move=False, 
				add_string=''
				)
		elif wave_mode=='bfield' :
			b_type = self._wave_api._bfield_paras.field_type.get()
			bfield = self._wave_api._bfield_paras.bfield.get()
			if b_type == 'none' :
				return
			if bfield is None :
				return
			if b_type == 'By' :
				# unitsConv=None
				# if not ( bfield._unitsXB is None ) :
				# 	unitsConv = [ bfield._unitsXB[0]/0.001, bfield._unitsXB[1]/1.0 ]
				bfield.write_field_std(
					file=wave_folder + 'stage/btab.dat',
					unitsConv=bfield._unitsXB,
					whatStr='By',
					)
			elif b_type == 'Byz' :
				bfield.write_field_waveByz(
					filey=wave_folder+'stage/btab.dat',
					filez=wave_folder+'stage/bz.dat',
					unitsXB=None
					)
			elif b_type == 'Bxyz' :
				bfield.write_field_waveByz(
					filex=wave_folder+'stage/bx.dat',
					filey=wave_folder+'stage/btab.dat',
					filez=wave_folder+'stage/bz.dat',
					unitsXB=None
					)
			elif b_type == 'bmap' :
				bfield.write_field_map_wave(
					file=wave_folder + 'stage/bmap.ntup',
					unitsXB=None
					)
		if wave_paras.observationPnts.get() == 1 :
			obsrvtnFile=wave_paras.observationFile.get()
			shutil.copyfile(
				f_h.convert_path_to_win(obsrvtnFile), 
				f_h.convert_path_to_win(wave_folder + 'stage/observ.in')
				)

		# pdb.set_trace()	