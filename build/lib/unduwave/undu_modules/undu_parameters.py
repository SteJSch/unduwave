from unduwave.unduwave_incl import *
from unduwave.attribute_classes.attributes import _attribute
from unduwave.attribute_classes.attributes import _attribute_collection
import unduwave.helpers.file_folder_helpers as f_h

class undu_prog_parameters(_attribute_collection):
	"""
	Represents standard parameters for undumag simulations.
	"""

	undumag_prog_folder = _attribute('')
	in_file_folder = _attribute('')
	in_file_nam = _attribute('')
	in_file_clc = _attribute('')
	in_file_clc_raw = _attribute('')
	in_file_clc_lines = _attribute('')
	copy_clc_folder = _attribute('')
	undu_mode = _attribute('')

	res_folder = _attribute('')
	undu_data_res_folder = _attribute('')
	pics_folder = _attribute('')
	res_summary_file = _attribute('')
	no_copy = _attribute([])
	undu_ending_extract = _attribute([])
	undu_ending_copy = _attribute([])
	undu_files_essentials = _attribute([])
	zipped = _attribute(False)
	undu_res_copy_behaviour = _attribute('')
	zip_res_folder = _attribute(0)

	nthreads = _attribute(1,in_name='nuthreads') # Number of threads to use
	writeGeometry = _attribute(1,in_name='iundugeo') # write geometry file, < 0 -> stop after writing
	plotGeometry = _attribute(1,in_name='iunduplot') # plot geometry, < 0 -> stop after plotting
	convergence_iron_residuals = _attribute(1,in_name='resiron') # convergence criterion, if iron-residual rms < convergence_iron_reiduals -> stop
	convergence_relative_b = _attribute(1,in_name='hconv') # convergence criterion, if rel. change of b-field less, stop
	beff_number_x = _attribute(1,in_name='nxbeff') # number of points at which to calc beff
	create_z_sym = _attribute(1,in_name='izsym') # mirror at x-y plane
	create_x_sym = _attribute(1,in_name='ixsym') # mirror at z-y plane
	create_y_sym = _attribute(1,in_name='iysym') # mirror at x-z plane
	bmap_dx = _attribute(1,in_name='dxmap') # step size for field map in x
	bmap_z_min = _attribute(1,in_name='zmapmin')
	bmap_z_max = _attribute(1,in_name='zmapmax')
	bmap_nz = _attribute(1,in_name='nzmap')
	bmap_x_min = _attribute(1,in_name='xmapmin')
	bmap_x_max = _attribute(1,in_name='xmapmax')
	bmap_nx = _attribute(1,in_name='nxmap')
	bmap_y_min = _attribute(1,in_name='ymapmin')
	bmap_y_max = _attribute(1,in_name='ymapmax')
	bmap_ny = _attribute(1,in_name='nymap')
	calc_force = _attribute(1,in_name='iforce')
	plot_force = _attribute(1,in_name='iplforce')
	force_box_center_x = _attribute(1,in_name='ubfcenx')
	force_box_center_y = _attribute(1,in_name='ubfceny')
	force_box_center_z = _attribute(1,in_name='ubfcenz')
	force_box_len_x = _attribute(1,in_name='ubflenx')
	force_box_len_y = _attribute(1,in_name='ubfleny')
	force_box_len_z = _attribute(1,in_name='ubflenz')
	force_segm_x = _attribute(1,in_name='mbforcex')
	force_segm_y = _attribute(1,in_name='mbforcey')
	force_segm_z = _attribute(1,in_name='mbforcez')
	center_magnet_struct = _attribute(1,in_name='kxcenter') # center magnetic struct
	bmap_not_inside_magn = _attribute(1,in_name='knomagmap')
	bmap_not_inside_pol = _attribute(1,in_name='knopolmap')
	kurad = _attribute(1,in_name='kurad')
	force_box_segm_y_distro = _attribute(1,in_name='ndivfboxy')
	kpreset = _attribute(1,in_name='kpreset')
	matrix = _attribute(1,in_name='matrix') # setting matrix mode (>1, calc and write interaction matrix,<0 read from file)
	undu_mode = _attribute('')

	def get_std_paras(self, undu_mode = 'from_clc_file' ): 
		"""
		"""

		#std-bad!
		dir_path = os.path.dirname(os.path.realpath(__file__))
		self.undumag_prog_folder.set(dir_path+f_h.convert_path_to_win('/../../External-Software/Undumag/'))
		self.in_file_folder.set(dir_path+f_h.convert_path_to_win('/../UNDWAVE_IN_FILES/Undu-In-Files/'))
		# dir_path='/var/lib/unduwave'
		# self.wave_prog_folder.set(dir_path+f_h.convert_path_to_win('/External-Software/Undumag/'))
		# self.in_file_folder.set(dir_path+f_h.convert_path_to_win('/unduwave/UNDWAVE_IN_FILES/Undu-In-Files/'))

		self.in_file_nam.set('undumag.nam')
		self.in_file_clc.set('undumag.clc')
		self.in_file_clc_raw.set('undu_raw.clc')
		self.undu_mode.set(undu_mode)

		self.res_folder.set('')
		self.undu_data_res_folder.set('UNDU_DATA/')
		self.pics_folder.set('Pics/')
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
			'urad_traxyz.dat',
			'undumag_trajectory.eps', 
			'undumag_y_z.eps', 
			'undumag_on-axis_by_bz.eps',
		 	'undumag.eps', 
			'undumag_on-axis_byint_bzint.eps', 
			'undumag_on-axis_byint2_bzint2.eps',
			'undumag.map',
			'undumag.frc',
			])
		self.zipped.set(False)
		self.undu_res_copy_behaviour.set('copy_essentials')
		self.zip_res_folder.set(0)
		self.in_file_clc_lines.set([])
		self.copy_clc_folder.set('')

		self.nthreads.set(6)
		self.writeGeometry.set(0)
		self.plotGeometry.set(0)
		self.convergence_iron_residuals.set(1.0e-4)
		self.convergence_relative_b.set(1.0e-6)
		self.beff_number_x.set(101)
		self.create_z_sym.set(1)
		self.create_x_sym.set(0)
		self.create_y_sym.set(1)
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
		self.force_box_segm_y_distro.set(1)
		self.kpreset.set(0)
		self.matrix.set(2)

		# if self.undu_mode.get() == 'from_clc_file' :

		# elif self.undu_mode.get() == 'from_undu_magns' :

		return self

