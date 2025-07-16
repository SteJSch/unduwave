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
	writeGeo = _attribute(1,in_name='iundugeo')
	plotGeo = _attribute(1,in_name='iunduplot')
	resiron = _attribute(1,in_name='resiron')
	hconv = _attribute(1,in_name='hconv')
	nxbeff = _attribute(1,in_name='nxbeff')
	izsym = _attribute(1,in_name='izsym')
	ixsym = _attribute(1,in_name='ixsym')
	iysym = _attribute(1,in_name='iysym')
	dxmap = _attribute(1,in_name='dxmap')
	zmapmin = _attribute(1,in_name='zmapmin')
	zmapmax = _attribute(1,in_name='zmapmax')
	nzmap = _attribute(1,in_name='nzmap')
	xmapmin = _attribute(1,in_name='xmapmin')
	xmapmax = _attribute(1,in_name='xmapmax')
	nxmap = _attribute(1,in_name='nxmap')
	dxmap = _attribute(1,in_name='dxmap')
	ymapmin = _attribute(1,in_name='ymapmin')
	ymapmax = _attribute(1,in_name='ymapmax')
	nymap = _attribute(1,in_name='nymap')
	iforce = _attribute(1,in_name='iforce')
	iplforce = _attribute(1,in_name='iplforce')
	ubfcenx = _attribute(1,in_name='ubfcenx')
	ubfceny = _attribute(1,in_name='ubfceny')
	ubfcenz = _attribute(1,in_name='ubfcenz')
	ubflenx = _attribute(1,in_name='ubflenx')
	ubfleny = _attribute(1,in_name='ubfleny')
	ubflenz = _attribute(1,in_name='ubflenz')
	mbforcex = _attribute(1,in_name='mbforcex')
	mbforcey = _attribute(1,in_name='mbforcey')
	mbforcez = _attribute(1,in_name='mbforcez')
	kxcenter = _attribute(1,in_name='kxcenter')
	knomagmap = _attribute(1,in_name='knomagmap')
	knopolmap = _attribute(1,in_name='knopolmap')
	kurad = _attribute(1,in_name='kurad')
	ndivfboxy = _attribute(1,in_name='ndivfboxy')
	kpreset = _attribute(1,in_name='kpreset')
	matrix = _attribute(1,in_name='matrix')
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
		self.writeGeo.set(0)
		self.plotGeo.set(0)
		self.resiron.set(1.0e-4)
		self.hconv.set(1.0e-6)
		self.nxbeff.set(101)
		self.izsym.set(1)
		self.ixsym.set(0)
		self.iysym.set(1)
		self.dxmap.set(1.0)
		self.zmapmin.set(-30)
		self.zmapmax.set(30)
		self.nzmap.set(1)
		self.xmapmin.set(9999.)
		self.xmapmax.set(9999.)
		self.nxmap.set(1)
		self.dxmap.set(1)
		self.ymapmin.set(-30)
		self.ymapmax.set(30)
		self.nymap.set(1)
		self.iforce.set(0)
		self.iplforce.set(0)
		self.ubfcenx.set(0.0)
		self.ubfceny.set(0.0)
		self.ubfcenz.set(0.0)
		self.ubflenx.set(0.0)
		self.ubfleny.set(0.0)
		self.ubflenz.set(0.0)
		self.mbforcex.set(20)
		self.mbforcey.set(20)
		self.mbforcez.set(20)
		self.kxcenter.set(1)
		self.knomagmap.set(0)
		self.knopolmap.set(0)
		self.kurad.set(0)
		self.ndivfboxy.set(1)
		self.kpreset.set(0)
		self.matrix.set(2)

		# if self.undu_mode.get() == 'from_clc_file' :

		# elif self.undu_mode.get() == 'from_undu_magns' :

		return self

