from unduwave.unduwave_incl import *
from unduwave.wave_modules.wave_parameters import *
import unduwave.analytical_module.ana_undu.bfield as bfield

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
		self._wave_api._undu_paras.update_values(ebeam=self._wave_api._ebeam_paras)

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
		b_type = self._wave_api._bfield_paras.field_type.get()
		if b_type == 'none' :
			return

		wave_folder = self._wave_api._prog_paras.wave_prog_folder.get()
		bfield = self._wave_api._bfield_paras.bfield.get()
		# field_files = self._wave_api._prog_paras.field_files.get()
		# field_folder = self._wave_api._prog_paras.field_folder.get()

		if b_type == 'By' :
			unitsConv=None
			if not ( bfield._unitsXB is None ) :
				unitsConv = [ bfield._unitsXB[0]/0.001, bfield._unitsXB[1]/1.0 ]
			bfield.write_field_std(
				file=wave_folder + 'stage/btab.dat',
				unitsConv=unitsConv,
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

		pdb.set_trace()
		return

		# except for Bxyz, By is first field file
		y_pos = 0
		if (b_type == 'Byz') or (b_type == 'Bxyz') : 
			# last file is Bz
			field_file = field_files[-1]
			bfield.convert_x_mm_b_T_file_to_wave_std( folder_in = field_folder, file_in = field_file, out_path = wave_folder + 'stage/bz.dat' )
			if b_type == 'Bxyz' :
				# snd file is By
				y_pos = 1
				# first file is Bx
				field_file = field_files[0]
				bfield.convert_x_mm_b_T_file_to_wave_std( folder_in = field_folder, file_in = field_file, out_path = wave_folder + 'stage/bx.dat' )
			field_file = field_files[y_pos]
			bfield.convert_x_mm_b_T_file_to_wave_std( folder_in = field_folder, file_in = field_file, out_path = wave_folder + 'stage/btab.dat' )
		elif (b_type == 'By') :
			# first file is By
			if len(field_files) > 0 :
				field_file = field_files[0]
				field_file = field_file.replace("(", "\\(")
				field_file = field_file.replace(")", "\\)")    
				if os.name == 'nt' :
					shutil.copyfile(field_folder + field_file, wave_folder + 'stage/')
					shutil.move(wave_folder + 'stage/' + field_file, wave_folder + 'stage/btab.dat')
				else:
					os.system( 'cp ' + field_folder + field_file + ' ' + wave_folder + 'stage/' )
					os.system( 'mv ' + wave_folder + 'stage/' + field_file + ' ' + wave_folder + 'stage/btab.dat' )
		elif (b_type=='bmap') :
				field_file = field_files[-1]
				field_file = field_file.replace("(", "\\(")
				field_file = field_file.replace(")", "\\)")    
				if os.name == 'nt' :
					shutil.copyfile(field_folder + field_file, wave_folder + 'stage/')
					shutil.move(wave_folder + 'stage/' + field_file, wave_folder + 'stage/bmap.ntup')
				else:
					os.system( 'cp ' + field_folder + field_file + ' ' + wave_folder + 'stage/' )
					os.system( 'mv ' + wave_folder + 'stage/' + field_file + ' ' + wave_folder + 'stage/bmap.ntup' )

