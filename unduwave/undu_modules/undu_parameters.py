from unduwave.unduwave_incl import *
from unduwave.attribute_classes.attributes import _attribute
from unduwave.attribute_classes.attributes import _attribute_collection
import unduwave.helpers.file_folder_helpers as f_h
import unduwave.undu_modules.undu_blocks as undu_blocks

def load_old_para(folder=None,file='para.txt') : 
	para={}
	with open(folder+file) as para_file:
		line_num = 0
		for line in para_file :
			if line_num == 0 : 
				line_num = line_num + 1
				continue
			newline = line.rstrip("\n").split("\t")
			name = newline[0]				
			val = newline[1]
			try:
				para.update({name : int(val)})
			except:
				try : 
					para.update({name : float(val)})
				except : 
					a = str(val)
					if a == 'False' : 
						val = False
					elif a == 'True' : 
						val = True
					elif a == 'None' :
						val = None
					para.update({name : val})
	return para

def load_old_para(folder=None,file='para.txt') : 
	para={}
	with open(folder+file) as para_file:
		line_num = 0
		for line in para_file :
			if line_num == 0 : 
				line_num = line_num + 1
				continue
			newline = line.rstrip("\n").split("\t")
			name = newline[0]				
			val = newline[1]
			try:
				para.update({name : int(val)})
			except:
				try : 
					para.update({name : float(val)})
				except : 
					a = str(val)
					if a == 'False' : 
						val = False
					elif a == 'True' : 
						val = True
					elif a == 'None' :
						val = None
					para.update({name : val})
	return para

class undu_prog_parameters(_attribute_collection):
	"""
	Represents standard parameters for undumag simulations.
	"""

	def __init__(self) :
		self.undumag_prog_folder = _attribute('')
		self.undumag_stage_folder = _attribute('')
		self.undumag_curr_folder = _attribute('')
		self.in_file_folder = _attribute('')

		self.material_files_std_folders = _attribute([])
		self.magnetic_materials = _attribute([])

		self.in_file_nam = _attribute('')
		self.in_file_clc = _attribute('')
		self.in_file_clc_raw = _attribute('')
		self.in_file_clc_lines = _attribute('')
		self.copy_clc_folder = _attribute('')
		self.undu_mode = _attribute('')

		self.res_folder = _attribute('')
		self.undu_data_res_folder = _attribute('')
		self.pics_folder = _attribute('')
		self.res_summary_file = _attribute('')
		self.no_copy = _attribute([])
		self.undu_ending_extract = _attribute([])
		self.undu_ending_copy = _attribute([])
		self.undu_files_essentials = _attribute([])
		self.zipped = _attribute(False)
		self.undu_res_copy_behaviour = _attribute('')
		self.zip_res_folder = _attribute(0)

		self.periodLength =_attribute(0)
		self.nthreads = _attribute(1,in_name='nuthreads') # Number of threads to use
		self.writeGeometry = _attribute(1,in_name='iundugeo') # write geometry file, < 0 -> stop after writing
		self.plotGeometry = _attribute(1,in_name='iunduplot') # plot geometry, < 0 -> stop after plotting
		self.convergence_iron_residuals = _attribute(1,in_name='resiron') # convergence criterion, if iron-residual rms < convergence_iron_reiduals -> stop
		self.convergence_relative_b = _attribute(1,in_name='hconv') # convergence criterion, if rel. change of b-field less, stop
		self.beff_number_x = _attribute(1,in_name='nxbeff') # number of points at which to calc beff
		self.create_z_sym = _attribute(1,in_name='izsym') # mirror at x-y plane
		self.create_x_sym = _attribute(1,in_name='ixsym') # mirror at z-y plane
		self.create_y_sym = _attribute(1,in_name='iysym') # mirror at x-z plane
		self.bmap_dx = _attribute(1,in_name='dxmap') # step size for field map in x
		self.bmap_z_min = _attribute(1,in_name='zmapmin')
		self.bmap_z_max = _attribute(1,in_name='zmapmax')
		self.bmap_nz = _attribute(1,in_name='nzmap')
		self.bmap_x_min = _attribute(1,in_name='xmapmin')
		self.bmap_x_max = _attribute(1,in_name='xmapmax')
		self.bmap_nx = _attribute(1,in_name='nxmap')
		self.bmap_y_min = _attribute(1,in_name='ymapmin')
		self.bmap_y_max = _attribute(1,in_name='ymapmax')
		self.bmap_ny = _attribute(1,in_name='nymap')
		self.calc_force = _attribute(1,in_name='iforce')
		self.plot_force = _attribute(1,in_name='iplforce')
		self.force_box_center_x = _attribute(1,in_name='ubfcenx')
		self.force_box_center_y = _attribute(1,in_name='ubfceny')
		self.force_box_center_z = _attribute(1,in_name='ubfcenz')
		self.force_box_len_x = _attribute(1,in_name='ubflenx')
		self.force_box_len_y = _attribute(1,in_name='ubfleny')
		self.force_box_len_z = _attribute(1,in_name='ubflenz')
		self.force_segm_x = _attribute(1,in_name='mbforcex')
		self.force_segm_y = _attribute(1,in_name='mbforcey')
		self.force_segm_z = _attribute(1,in_name='mbforcez')
		self.center_magnet_struct = _attribute(1,in_name='kxcenter') # center magnetic struct
		self.bmap_not_inside_magn = _attribute(1,in_name='knomagmap')
		self.bmap_not_inside_pol = _attribute(1,in_name='knopolmap')
		self.kurad = _attribute(1,in_name='kurad')
		self.shuffle = _attribute(0,in_name='kshuffle')
		self.force_box_segm_y_distro = _attribute(1,in_name='ndivfboxy')
		self.kpreset = _attribute(1,in_name='kpreset')
		self.magnetCoating=_attribute(0.0,in_name='Mcoating')
		self.matrix = _attribute(1,in_name='matrix') # setting matrix mode (>1, calc and write interaction matrix,<0 read from file)
		self.undu_mode = _attribute('')
		super().__init__()

	def load_all_magnetic_material_files_from_folders(self,folders=None) :
		if folders is None :
			folders=self.material_files_std_folders()
		for folder in folders :
			directory = os.fsencode(folder)
				
			for file in os.listdir(directory):
				filename = os.fsdecode(file)
				if filename.endswith(".dat"): 
					fullFile=os.path.join(directory, file)
					mat=undu_blocks.magnetic_material(
						material_id=None, 
						base_material_type=None, # "magnet" or "pole"
						ksi_easy=None,
						ksi_perp=None,
						magnetization_data=None,
						magnetization_file=fullFile,
						)
					self.magnetic_materials.set(undu_blocks.magnetic_material.add_new_to_material_list(
							theList=self.magnetic_materials(),
							other=mat,
							)[0]
						)

		materialList=self.magnetic_materials()
		if materialList[0]._material_id().find('fm_vanadium_permendur')>=0:
			materialList.append(materialList[0])
			materialList.pop(0)

	def get_std_paras(self, undu_mode = 'from_clc_file' ): 
		"""
		"""

		self.undumag_curr_folder.set(Path(''))
		self.undumag_prog_folder.set(ROOT_DIR/"External-Software"/"UNDUMAG")
		self.undumag_stage_folder.set(ROOT_DIR/'External-Software'/'pristine_stages'/'undu')
		self.in_file_folder.set(ROOT_DIR/"UNDWAVE_IN_FILES"/"Undu-In-Files")

		self.material_files_std_folders.set([ROOT_DIR/"MATERIAL_FILES"])
		self.magnetic_materials = _attribute([])

		self.load_all_magnetic_material_files_from_folders(
			folders=self.material_files_std_folders()
			)

		self.in_file_nam.set('undumag.nam')
		self.in_file_clc.set('undumag.clc')
		self.in_file_clc_raw.set('undu_raw.clc')
		self.undu_mode.set(undu_mode)

		self.res_folder.set(Path(''))
		self.undu_data_res_folder.set(Path('undu_res'))
		self.pics_folder.set('')
		self.res_summary_file.set('res_summary.txt')
		self.no_copy.set(['WAVE_CODE.DAT', 'undumag_mu_77K.dat', 'undumag_mu_300K.dat',
								'iron_muinf_sat-2.34.dat', 'Vanadium_Permendur_Radia', 
								])
		self.undu_ending_extract.set([ 
				'dat', 
				'map', 
				'nb', 
				'py', 
				'txt', 
				'frc', 
				'beff',
				])
		self.undu_ending_copy.set([ 'clc', 'nam' ])
		self.undu_files_essentials.set([ 
			'undumag_on-axis.dat', 
			'undumag.beff', 
			'undumag_field_profile.dat', 
			'undumag_field_profile.eps',
			'urad_traxyz.dat',
			'undumag_trajectory.eps', 
			'undumag_y_z.eps', 
			'undumag_on-axis_by_bz.eps',
			'undumag.eps', 
			'undumag_on-axis_byint_bzint.eps', 
			'undumag_on-axis_byint2_bzint2.eps',
			'undumag.map',
			'undumag.frc',
			'undumag_mh_iron',
			'undumag_msh_radia.py',
			])
		self.zipped.set(False)
		self.undu_res_copy_behaviour.set('copy_essentials')
		self.zip_res_folder.set(0)
		self.in_file_clc_lines.set([])
		self.copy_clc_folder.set(Path(''))

		self.nthreads.set(6)
		self.periodLength.set(0.02)
		self.writeGeometry.set(0)
		self.plotGeometry.set(0)
		self.convergence_iron_residuals.set(1.0e-4)
		self.convergence_relative_b.set(1.0e-6)
		self.beff_number_x.set(101)
		self.create_z_sym.set(0)
		self.create_x_sym.set(0)
		self.create_y_sym.set(0)
		self.bmap_dx.set(1.0)
		self.bmap_z_min.set(-30)
		self.bmap_z_max.set(30)
		self.bmap_nz.set(1)
		self.bmap_x_min.set(9999.)
		self.bmap_x_max.set(9999.)
		self.bmap_nx.set(1)
		self.bmap_y_min.set(-30)
		self.bmap_y_max.set(30)
		self.bmap_ny.set(1)
		self.calc_force.set(0)
		self.plot_force.set(0)
		self.force_box_center_x.set(0.0)
		self.force_box_center_y.set(0.0)
		self.force_box_center_z.set(0.0)
		self.force_box_len_x.set(0.0)
		self.force_box_len_y.set(0.0)
		self.force_box_len_z.set(0.0)
		self.force_segm_x.set(20)
		self.force_segm_y.set(20)
		self.force_segm_z.set(20)
		self.center_magnet_struct.set(1)
		self.bmap_not_inside_magn.set(0)
		self.bmap_not_inside_pol.set(0)
		self.kurad.set(0)
		self.magnetCoating.set(0.014)
		self.shuffle.set(0)
		self.force_box_segm_y_distro.set(1)
		self.kpreset.set(0)
		self.matrix.set(2)

		# if self.undu_mode.get() == 'from_clc_file' :

		# elif self.undu_mode.get() == 'from_undu_magns' :

		return self

