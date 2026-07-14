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

		# pdb.set_trace()

		# traj_y.plot_over(x_quant=traj_x,nfig=nfig,nosave=True)
		# nfig=traj_z.plot_over(x_quant=traj_x,nfig=nfig)

		# By.plot_over(x_quant=traj_x,nfig=nfig,nosave=True)
		# nfig=Bz.plot_over(x_quant=traj_x,nfig=nfig)

		# power_z = results.get_result(which='power_z')
		# power_y = results.get_result(which='power_y')
		# power_distro = results.get_result(which='power_distribution')
		# nfig=power_distro.plot_over_3d(x_quant=power_y,y_quant=power_z,file_name=None,nosave=False,nfig=nfig)

		# en_flux = results.get_result(which='en_flux')
		# flux = results.get_result(which='flux')
		# nfig=flux.plot_over(x_quant=en_flux,file_name=None,nosave=False,nfig=nfig,loglog=False)

		# en_brill = results.get_result(which='en_brill')
		# brill0 = results.get_result(which='brill0')
		# brill0e = results.get_result(which='brill0e')
		# brill0f = results.get_result(which='brill0f')
		# brill0ef = results.get_result(which='brill0ef')
		# brill0.plot_over(x_quant=en_brill,nfig=nfig,nosave=True,loglog=True)
		# brill0e.plot_over(x_quant=en_brill,nfig=nfig,nosave=True,loglog=True)
		# brill0f.plot_over(x_quant=en_brill,nfig=nfig,nosave=True,loglog=True)
		# nfig=brill0ef.plot_over(x_quant=en_brill,nfig=nfig,loglog=True)

		# en_fd = results.get_result(which='en_fd')
		# flux_density_onaxis = results.get_result(which='flux_density')
		# nfig=flux_density_onaxis.plot_over(x_quant=en_fd,nfig=nfig,loglog=False)

		# flux_dens_distr_ens_loaded = results.find_load_flux_density_distribution(energies=[anaUndu.firstHarmEnergyStrong.get()])
		# fd_y = results.get_result(which='fd_y')
		# fd_z = results.get_result(which='fd_z')
		# for en in flux_dens_distr_ens_loaded :
		# 	fd = results.get_result(which=f'flux_density_distribution_{en:.2f}')
		# 	nfig=fd.plot_over_3d(x_quant=fd_y,y_quant=fd_z,file_name=None,nosave=False,nfig=nfig)

		# pdb.set_trace()

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
