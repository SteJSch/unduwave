import os
import pdb
import sys
import unittest

sys.path.insert(0, '../../../')

try :
	# works when calling script with python3 script_file
	dir_path = os.path.dirname(os.path.realpath(__file__))
except:
	# works when calling script with exec from python console
	dir_path = os.getcwd()

import unduwave as uw
from unduwave.unduwave_incl import *
from unduwave import undu_magnets
from unduwave import undu_representation

class undulatorObject(unittest.TestCase) :

	def test_undulatorObject(self) :

		res_folder='res'
		res_folder_full=dir_path+f'/{res_folder}/'

		undu = uw.undu(undu_mode='from_undu_magns')
		undu_prog_paras = undu._prog_paras
		undu_prog_paras.res_folder.set(res_folder_full)
		undu_prog_paras.plotGeometry.set(1)
		undu_prog_paras.create_z_sym.set(0)
		undu_prog_paras.create_y_sym.set(0)
		undu_prog_paras.bmap_nz.set(10)
		undu_prog_paras.bmap_nx.set(10)
		undu_prog_paras.nthreads.set(15)
		undu_prog_paras.center_magnet_struct.set(0)

		undu_center = undu_magnets.point_coords(x=0.0,y=0.0,z=0.0)

		def myMagnFunComplex(center,magnParas,name) :
			name_main = combine_strings(str_add='main',str_add_to=name)
			name_side = combine_strings(str_add='side',str_add_to=name)
			main_block = undu_magnets.undu_magnet_block_coords(
				p_center=center,
				name=name_main,
				mother=name,
				magnetParameters=magnParas,
				)
			# magnParasTmp=copy.deepcopy(magnParas)
			# p_side = copy.deepcopy(center)
			center._y=center._y - magnParas._len_y_main/2.0+magnParas._len_y_side/2.0
			center._z=center._z - magnParas._len_z_main/2.0-magnParas._len_z_side/2.0
			magnParas._len_y_main=magnParas._len_y_side
			magnParas._len_z_main=magnParas._len_z_side
			side_block = undu_magnets.undu_magnet_block_coords(
				p_center=center,
				name=name_side,
				mother=name,
				magnetParameters=magnParas,
				)
			return undu_magnets.undu_magnets(magnet_blocks=[main_block,side_block])

		magnParas=undu_magnets.magnetParameters(
			len_x_main=3, 
			len_y_main=20,
			len_z_main=25, 
			magnetization=1.22,
			magn_unit_vec='+x',
			segm_x=2,
			segm_y=3,
			segm_z=4,
			material="mag", #"mag","pol", 'NiCuFoil'
		)
		magnParas._len_y_side=5
		magnParas._len_z_side=2
		polParas=undu_magnets.magnetParameters(
			len_x_main=6,
			len_y_main=15, 
			len_z_main=35, 
			segm_x=2,
			segm_y=3,
			segm_z=3,
			frac_y=2,
			frac_z=3,
			material="pol", #"mag","pol", 'NiCuFoil'
		)
		polParas._len_y_side=6
		polParas._len_z_side=5

		magnetCenter=undu_magnets.point_coords()
		# magnetCenter._z=polParas._len_z_main/2.0-magnParas._len_z_main/2.0
		myMagnet=undu_representation.magnetRepresentation(
			center=magnetCenter,
			magnParas=magnParas,
			magnFun=myMagnFunComplex
			)

		polCenter=undu_magnets.point_coords()
		polCenter._y=magnParas._len_y_main/2.0-polParas._len_y_main/2.0
		# polCenter._z=magnParas._len_z_main/2.0-polParas._len_z_main/2.0
		myPol=undu_representation.magnetRepresentation(
			center=polCenter,
			magnParas=polParas,
			magnFun=myMagnFunComplex
			)

		undulatorMagMag=undu_representation.unduRepresentation(
			gap = 5.0, 
			shift=10,
			glueSlit = 1.0, 
			keeperSlit=2.0, 
			rowSlit=1.0, 
			nperiods = 3,
			unduType='ellipt', # "planar", 'ellipt' 
			magnetsPerKeeper=[myMagnet,myMagnet],
			endMagnets=[myMagnet],
			)
		undulatorMagPol=undu_representation.unduRepresentation(
			gap = 5.0, 
			shift=0.0,
			glueSlit = 1.0, 
			keeperSlit=2.0, 
			rowSlit=1.0, 
			nperiods = 5,
			unduType='planar', # "planar", 'ellipt' 
			magnetsPerKeeper=[myMagnet,myPol],
			endMagnets=[myMagnet,myPol],
			)
		undulatorRepresentation=undulatorMagMag# undulatorMagMag, undulatorMagPol
		undulator=undulatorRepresentation.createUndulator(
			center=undu_magnets.point_coords(), 
			magnetizations=[1.22,1.33],
			magnSeq=['-x','-y','+x','+y'], 
			nameCore='',
			onlyLL=False,
			)
		undulator.add_to_clc(api=undu)

		unduObj=uw.undulator.undulatorObj()
		unduObj.fromRepres(
			undulatorRepres=undulator,
			periodLength=None,
			nPeriods=None,
			unduAPI=None,
			setSym=True,
			simulate=True,
			beffs=None,
			resFolder=res_folder_full,
			) 

		pdb.set_trace()

		undu.run()

		wave = uw.wave(wave_mode='bfield')

		bfield_paras = wave._bfield_paras # get bfield paras
		bfield_paras.bfield.set(bfieldUE51) # set the bfield
		bfield_paras.field_type.set('By') # tell wave which part of bfield to use for simu

		"""
		Setting Program Parameter
		"""

		wave_prog_paras = wave._prog_paras
		wave_prog_paras.res_folder.set(res_folder)
		wave_prog_paras.calc_spectrum.set(True)
		wave_prog_paras.nthreads.set(6)
		wave_prog_paras.calc_emittance.set(1)

		"""
		Setting Spectrometer Parameter
		"""

		spectrometer_paras = wave._spectrometer_paras
		spectrometer_paras.spectrum_n_energies.set(100)
		spectrometer_paras.spectrum_min_energy.set(50)
		spectrometer_paras.spectrum_max_energy.set(200)
		spectrometer_paras.spectrum_undu_mode.set(1)

		"""
		Setting Screen Parameter
		"""
		screen_paras = wave._screen_paras
		screen_paras.screen_segm_hor.set(30) 
		screen_paras.screen_segm_vert.set(30)
		screen_paras.screen_extent_hor.set(40) # pinhole width mm
		screen_paras.screen_extent_vert.set(40) # pinhole height mm
		screen_paras.screen_pos_x.set(10)
		screen_paras.screen_pos_y.set(0.0)
		screen_paras.screen_pos_z.set(0.0)

		"""
		Setting Undulator Parameter
		"""

		undu_paras = wave._undu_paras # getting parameter object
		undu_paras.elliptUnduB0Y.set(1.18)
		undu_paras.elliptUnduB0Z.set(0.0)
		undu_paras.elliptUnduNumPeriods.set(10)
		undu_paras.elliptUnduPerLength.set(0.051)
		undu_paras.elliptUnduPerShift.set(0.5)

		"""
		Setting Beam Parameter
		"""

		ebeam_paras = wave._ebeam_paras
		ebeam_paras.beam_en.set(1.722) # [GeV]
		ebeam_paras.current.set(0.3) # [A]
		ebeam_paras.beamSizeHor.set(275e-6) # 
		ebeam_paras.beamDiveHor.set(28.1e-6) #
		ebeam_paras.beamSizeVer.set(22.5e-6) # 
		ebeam_paras.beamDiveVer.set(6.8e-6) # 
		ebeam_paras.espread.set(1e-3) # 
		ebeam_paras.emittanceHor.set(7.7e-9)
		ebeam_paras.emittanceVer.set(15.4e-11)
		ebeam_paras.update_values()

		"""
		Run
		"""

		wave.run()

		"""
		Get Results and Plot
		"""

		results = wave.get_results()
		nfig=0

		traj_x = results.get_result(which='traj_x')
		traj_y = results.get_result(which='traj_y')
		traj_z = results.get_result(which='traj_z')

		traj_y.plot_over(x_quant=traj_x,nfig=nfig,nosave=True,clear=True)
		nfig=traj_z.plot_over(x_quant=traj_x,nfig=nfig)

		By = results.get_result(which='By')
		By.plot_over(x_quant=traj_x,nfig=nfig,nosave=True,clear=True)
		Bz = results.get_result(which='Bz')
		nfig=Bz.plot_over(x_quant=traj_x,nfig=nfig)

		en_flux = results.get_result(which='en_flux')
		flux = results.get_result(which='flux')
		nfig=flux.plot_over(x_quant=en_flux,file_name=None,nosave=False,nfig=nfig,loglog=True,clear=True)

		en_fd = results.get_result(which='en_fd')
		flux_density_onaxis = results.get_result(which='flux_density')
		nfig=flux_density_onaxis.plot_over(x_quant=en_fd,nfig=nfig,loglog=True,clear=True)

		pdb.set_trace()
		power_z = results.get_result(which='power_z')
		power_y = results.get_result(which='power_y')
		power_distro = results.get_result(which='power_distribution')
		nfig=power_distro.plot_over_3d(x_quant=power_y,y_quant=power_z,file_name=None,nosave=False,nfig=nfig)

		en_brill = results.get_result(which='en_brill')
		brill0 = results.get_result(which='brill0')
		brill0e = results.get_result(which='brill0e')
		brill0f = results.get_result(which='brill0f')
		brill0ef = results.get_result(which='brill0ef')
		brill0.plot_over(x_quant=en_brill,nfig=nfig,nosave=True,loglog=True)
		brill0e.plot_over(x_quant=en_brill,nfig=nfig,nosave=True,loglog=True)
		brill0f.plot_over(x_quant=en_brill,nfig=nfig,nosave=True,loglog=True)
		nfig=brill0ef.plot_over(x_quant=en_brill,nfig=nfig,loglog=True)

		en_fd = results.get_result(which='en_fd')
		flux_density_onaxis = results.get_result(which='flux_density')
		nfig=flux_density_onaxis.plot_over(x_quant=en_fd,nfig=nfig,loglog=True)

		flux_dens_distr_ens_loaded = results.find_load_flux_density_distribution(energies=[414])
		fd_y = results.get_result(which='fd_y')
		fd_z = results.get_result(which='fd_z')
		for en in flux_dens_distr_ens_loaded :
			fd = results.get_result(which=f'flux_density_distribution_{en:.2f}')
			nfig=fd.plot_over_3d(x_quant=fd_y,y_quant=fd_z,file_name=None,nosave=False,nfig=nfig)

		pdb.set_trace()

if __name__ == '__main__':
	unittest.main()
