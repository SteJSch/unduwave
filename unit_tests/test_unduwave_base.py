import unittest

import unduwave as uw
from unduwave.unduwave_incl import *
from unduwave import undu_blocks
import unduwave.helpers.file_folder_helpers as f_h

try :
	# works when calling script with python3 script_file
	dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
except:
	# works when calling script with exec from python console
	dir_path = Path(os.getcwd())

class unduwave_base_test(unittest.TestCase) :

	def test_undu_run(self) :

		self.res_folder='res'
		self.res_folder_dir=dir_path/f'{self.res_folder}/'

		undu = uw.undu(undu_mode='from_undu_magns')
		undu_prog_paras = undu._prog_paras
		undu_prog_paras.res_folder.set(self.res_folder_dir)
		undu_prog_paras.plotGeometry.set(1)
		undu_prog_paras.create_z_sym.set(0)

		undu_prog_paras.bmap_z_min.set(-10)
		undu_prog_paras.bmap_z_max.set(10)
		undu_prog_paras.bmap_y_min.set(0)
		undu_prog_paras.bmap_y_max.set(0)

		undu_prog_paras.bmap_ny.set(1)
		undu_prog_paras.bmap_nz.set(10)
		undu_prog_paras.bmap_nx.set(10)

		undu_data_res_fold=undu_prog_paras.undu_data_res_folder()
		res_folder_full=self.res_folder_dir/undu_data_res_fold

		pos=np.array([ 0.0, -15.0, 0.0 ])
		magn_paras=undu_blocks.magParameters(
			len_x_main=10, 
			len_y_main=20,
			len_z_main=30, 
			segm_x=2,
			segm_y=2,
			segm_z=2,
			frac_y=1,
			frac_z=1,
			chamf=0.3,
			material_id="pm_rec",
		)

		magnObject = undu_blocks.undumagBlockObject(
			center=pos,
			magnParas=magn_paras,
			name='MyMagnet',
			parentName='',
			)

		undu.set_magnet_objects(magn_objects=magnObject)

		undu.run()

		results = undu.get_results()

		trajx = results.get_result(which='trajx')
		by = results.get_result(which='by')
		bz = results.get_result(which='bz')
		intBy = results.get_result(which='intBy')
		intBz = results.get_result(which='intBz')
		intBy2 = results.get_result(which='intBy2')
		intBz2 = results.get_result(which='intBz2')

		zProfile = results.get_result(which='profz')
		ByProfile = results.get_result(which='profBy')
		BzProfile = results.get_result(which='profBz')

		bmap = results.get_result(which='bmap')

		self.assertEqual(os.path.exists(res_folder_full),True)
		self.assertEqual(os.path.exists(res_folder_full/"undumag.beff"),True)
		self.assertEqual(os.path.exists(res_folder_full/"undumag.map"),True)
		self.assertEqual(os.path.exists(res_folder_full/"undumag_msh_radia.py"),True)
		self.assertEqual(os.path.exists(res_folder_full/"undumag_on-axis.dat"),True)

		unduBfield=uw.bfield.bfield(
			unitsXB=[0.001,1.0] # setting the units
			)

		unduBfield.load_field_from_file(
					file=res_folder_full/"undumag_on-axis.dat", 
					fieldMap=False,
					unduFile = True, 
					radiaFile=False,
					header=None,
					skiprows=None,
				)
		"""
		Getting wave
		"""
		wave = uw.wave(wave_mode='bfield')

		bfield_paras = wave._bfield_paras # get bfield paras
		bfield_paras.bfield.set(unduBfield) # set the bfield
		bfield_paras.field_type.set('By') # tell wave which part of bfield to use for simu

		"""
		Setting Program Parameter
		"""

		wave_prog_paras = wave._prog_paras
		wave_prog_paras.res_folder.set(self.res_folder_dir)
		wave_prog_paras.calc_spectrum.set(True)
		wave_prog_paras.nthreads.set(6)

		"""
		Setting Spectrometer Parameter
		"""

		spectrometer_paras = wave._spectrometer_paras
		spectrometer_paras.spectrum_n_energies.set(10)
		spectrometer_paras.spectrum_min_energy.set(500)
		spectrometer_paras.spectrum_max_energy.set(550)

		"""
		Setting Screen Parameter
		"""
		screen_paras = wave._screen_paras
		screen_paras.screen_segm_hor.set(5) 
		screen_paras.screen_segm_vert.set(5)
		screen_paras.screen_extent_hor.set(40) # pinhole width mm
		screen_paras.screen_extent_vert.set(40) # pinhole height mm

		"""
		Setting Beam Parameter
		"""

		ebeam_paras = wave._ebeam_paras
		ebeam_paras.get_std_bessy_II_paras()

		"""
		Run
		"""
		wave_data_res_fold=wave_prog_paras.wave_data_res_folder()
		wave_folder_full=self.res_folder_dir/wave_data_res_fold

		wave.run()

		"""
		Get Results and Plot
		"""

		results = wave.get_results()
		nfig=5

		traj_x = results.get_result(which='traj_x')
		traj_y = results.get_result(which='traj_y')
		traj_z = results.get_result(which='traj_z')
		By = results.get_result(which='By')
		Bz = results.get_result(which='Bz')

		self.assertEqual(os.path.exists(wave_folder_full),True)
		self.assertEqual(os.path.exists(wave_folder_full/"stokes_dist_emittance_espread_10th_energy.wva"),True)
		self.assertEqual(os.path.exists(wave_folder_full/"wave.out"),True)
		self.assertEqual(os.path.exists(wave_folder_full/"photon_flux_(pinhole)_248001.wvh"),True)
		self.assertEqual(os.path.exists(wave_folder_full/"irradiated_power_dist.wva"),True)

		self.tearDown()

	def test_undu_mat(self) :

		self.res_folder='res'
		self.res_folder_dir=dir_path/f'{self.res_folder}/'

		material_folder=dir_path/'test_material/'

		material_combos=[ 
			["super_material","pm_rec","fm_vanadium_permendur"], 
			["super_material","super_material"], 
			["pm_rec","fm_vanadium_permendur"], 
			["pm_rec","pm_rec"], 
			["pm_rec","super_material"], 
			["fm_vanadium_permendur","super_material"], 
			]

		for combo in material_combos :

			undu = uw.undu(undu_mode='from_undu_magns')
			undu_prog_paras = undu._prog_paras
			undu_prog_paras.res_folder.set(self.res_folder_dir)
			undu_prog_paras.plotGeometry.set(1)
			undu_prog_paras.create_z_sym.set(0)

			undu_prog_paras.bmap_z_min.set(-20)
			undu_prog_paras.bmap_z_max.set(20)
			undu_prog_paras.bmap_y_min.set(0)
			undu_prog_paras.bmap_y_max.set(0)

			undu_prog_paras.bmap_ny.set(1)
			undu_prog_paras.bmap_nz.set(10)
			undu_prog_paras.bmap_nx.set(10)

			mat_folders=undu_prog_paras.material_files_std_folders()
			mat_folders.append(material_folder)
			undu_prog_paras.material_files_std_folders.set(mat_folders)

			pos0=np.array([ -15.0, 0.0, 0.0 ])
			magnet_blocks=[]
			for ind,el in enumerate(combo) :

				magn_paras=undu_blocks.magParameters(
					len_x_main=10, 
					len_y_main=20,
					len_z_main=30, 
					segm_x=2,
					segm_y=2,
					segm_z=2,
					frac_y=1,
					frac_z=1,
					chamf=0.3,
					material_id=el,
				)

				magnObject = undu_blocks.undumagBlockObject(
					center=pos0,
					magnParas=magn_paras,
					name=f'MyMagnet{ind}',
					parentName='',
					api=undu,
					)
				magnet_blocks.append(magnObject)
				pos0=pos0 + np.array([ 15.0, 0.0, 0.0 ])

			objs=undu_blocks.undumagObjectList(
				magnet_blocks=magnet_blocks,
				name='magnet_pole',
				parentName='',
				center=np.array([0.0,0.0,0.0]),
				)

			undu.set_magnet_objects(magn_objects=objs)

			undu.run()

		self.tearDown()

	# def test_undu_run(self) :

	# 	self.res_folder='res'
	# 	self.res_folder_dir=dir_path/f'{self.res_folder}/'

	# 	ROOT_DIR=uw.ROOT_DIR/'../scripts/Wave_Examples'
	# 	wave_script_folders=[
	# 		ROOT_DIR/'Spec_From_BField/spec_from_by',
	# 		ROOT_DIR/'Spec_From_BField/spec_from_byz',
	# 		]

	# 	for script in wave_script_folders :

	# 		os.chdir(script)
	# 		pyfiles=f_h.find_files_exptn(folder=script, hints=['.py'], exptns=[])
	# 		for pyfile in pyfiles:
	# 			subprocess.call(['python', pyfile])

	# 			pdb.set_trace()

	# 	self.tearDown()

	def setUp(self):
		self.resource = "Resource allocated"
		print("Setting up test resources...")

	def tearDown(self):
		print("Cleaning up resource...")

		ROOT_DIR_TEST = Path(__file__).resolve().parent
		del_dirs=[
			ROOT_DIR_TEST/'__pycache__',
			ROOT_DIR_TEST/self.res_folder_dir,
			]

		for del_dir in del_dirs :
			if os.path.exists(del_dir) :		
				shutil.rmtree(del_dir)

if __name__ == '__main__':
	unittest.main()
