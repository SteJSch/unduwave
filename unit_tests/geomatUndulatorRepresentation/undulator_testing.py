import os
import pdb
import sys
import unittest
import numpy as np

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

class undulatorTests(unittest.TestCase) :

	def test_unduRepresentation(self) :
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
		undu_prog_paras.center_magnet_struct.set(0)

		undu_center = undu_magnets.point_coords(x=0.0,y=-10.0,z=-10)
		magn_unit_vec=undu_magnets.create_magnetization_unit_vec(magn_string='+x')

		magnParas=undu_magnets.magnetParameters(
			len_x_main=2.5,
			len_y_main=20, 
			len_z_main=30, 
			magnetization=0.6,
			magn_unit_vec=magn_unit_vec,
			segm_x=2,
			segm_y=3,
			segm_z=4,
			material="mag",
			)

		"""
		Creating 1 Magnet via the fundamental function
		"""

		# main_block = undu_magnets.undu_magnet_block_coords(
		# 	p_center=undu_center,
		# 	len_x=10,
		# 	len_y=20,
		# 	len_z=20,
		# 	magnetization=1.12,
		# 	magn_unit_vec=magn_unit_vec,
		# 	name="mag", # "mag", "pol"
		# 	mother="mom",
		# 	segm_x = 2, 
		# 	segm_y = 2, 
		# 	segm_z = 2, 
		# 	frac_y=1,
		# 	frac_z=1,
		# 	material = 'mag',
		# 	chamf = 0.3,
		# 	api=undu
		# 	)
		# main_block.add_to_clc(api=undu)

		# magnParas=undu_magnets.magnetParameters(
		# 	len_x_main=3, 
		# 	len_y_main=10, 
		# 	len_z_main=25, 
		# 	magnetization=1.22,
		# 	magn_unit_vec='+x',
		# 	segm_x=2,
		# 	segm_y=3,
		# 	segm_z=4,
		# 	material="mag", #"mag","pol", 'NiCuFoil'
		# )
		# myMagnetRepr=undu_representation.magnetRepresentation(
		# 	center=undu_center,
		# 	magnParas=magnParas,
		# 	magnFun=None
		# 	)
		# myMagnet=myMagnetRepr.createMagnet(name='mag')
		# myMagnet.add_to_clc(api=undu)

		# magnParasC1=undu_magnets.magnetParameters(
		# 	len_x_main=3, 
		# 	len_y_main=20, 
		# 	len_z_main=25, 
		# 	segm_x=2,
		# 	segm_y=3,
		# 	segm_z=4,
		# 	material="mag", #"mag","pol", 'NiCuFoil'
		# )
		# magnParasC1._len_y_side=10
		# magnParasC1._len_z_side=15

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

		# myMagnetRepr=undu_representation.magnetRepresentation(
		# 	center=undu_center,
		# 	magnParas=magnParasC1,
		# 	magnFun=myMagnFunComplex
		# 	)
		# myMagnet=myMagnetRepr.createMagnet(name='mag')
		# myMagnet.add_to_clc(api=undu)

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
			segm_y=2,
			segm_z=3,
			material="pol", #"mag","pol", 'NiCuFoil'
		)
		polParas._len_y_side=6
		polParas._len_z_side=5
		# def myMagnFun(center,magnParas,name) :
		# 	name_main = combine_strings(str_add='main',str_add_to=name)
		# 	p_main = copy.deepcopy(center)
		# 	main_block = undu_magnets.undu_magnet_block_coords(
		# 		p_center=p_main,
		# 		name=name_main,
		# 		mother=name,
		# 		magnetParameters=magnParas,
		# 		)
		# 	return undu_magnets.undu_magnets(magnet_blocks=[main_block])

		# def myPolFun(center,polParas,name) :
		# 	name_main = combine_strings(str_add='main',str_add_to=name)
		# 	p_main = copy.deepcopy(center)
		# 	main_block = undu_magnets.undu_magnet_block_coords(
		# 		p_center=p_main,
		# 		name=name_main,
		# 		mother=name,
		# 		magnetParameters=polParas,
		# 		)
		# 	return undu_magnets.undu_magnets(magnet_blocks=[main_block])

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

		undulatorMagPol=undu_representation.unduRepresentation(
			gap = 5.0, 
			shift=0.0,
			glueSlit = 1.0, 
			keeperSlit=2.0, 
			rowSlit=1.0, 
			nperiods = 3,
			unduType='planar', # "planar", 'ellipt' 
			magnetsPerKeeper=[myMagnet,myPol],
			endMagnets=[myMagnet,myPol],
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
		undulatorRepresentation=undulatorMagMag # undulatorMagMag, undulatorMagPol

		# keeperCenter=undu_magnets.point_coords(x=0.0,y=-magnParas._len_y_main/2.0,z=-magnParas._len_z_main/2.0)
		# keeper=undulatorMagPol.createKeeper(
		# 	center=keeperCenter, 
		# 	magnetizations=[1.22,1.33], 
		# 	magnSeq=['+y','+x','-y'],
		# 	startNum=0,
		# 	nameCore='keeper',
		# 	)
		# keeper.add_to_clc(api=undu)

		# periodCenter=undu_magnets.point_coords(x=0.0,y=-magnParas._len_y_main/2.0,z=-magnParas._len_z_main/2.0)
		# period=undulatorMagPol.createPeriod(
		# 		center=periodCenter, 
		# 		magnetizations=[1.22,1.33], 
		# 		magnSeq=['-x','-y','+x','+y'], 
		# 		nameCore='period_1',
		# 		)
		# period.add_to_clc(api=undu)

		# nPeriodCenter=undu_magnets.point_coords(x=0.0,y=-magnParas._len_y_main/2.0,z=-magnParas._len_z_main/2.0)
		# nperiods=undulatorRepresentation.createNPeriods(
		# 		center=nPeriodCenter, 
		# 		magnetizations=[1.22,1.33], 
		# 		magnSeq=['-x','-y','+x','+y'], 
		# 		nameCore='',
		# 		)
		# nperiods.move_it(vec=undu_magnets.point_coords(0.0,-undulatorRepresentation._gap/2.0,0.0))
		# nperiods.move_it(vec=undu_magnets.point_coords(0.0,0.0,-0.5))
		# nperiods.add_to_clc(api=undu)

		# llRow=undulatorRepresentation.createRow(
		# 		center=undu_center, 
		# 		magnetizations=[1.22,1.33], 
		# 		magnSeq=['-x','-y','+x','+y'], 
		# 		nameCore='ll',
		# 		)
		# llRow.add_to_clc(api=undu)

		undulator=undulatorRepresentation.createUndulator(
			center=undu_magnets.point_coords(), 
			magnetizations=[1.22,1.33],
			magnSeq=['-x','-y','+x','+y'], 
			nameCore='',
			onlyLL=False,
			)
		undulator.add_to_clc(api=undu)

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

		add=''

		nfig=0
		by.plot_over(
			x_quant=trajx,
			nfig=nfig,
			nosave=True,
			clear=True,
			dataFile=f'/{res_folder_full}by{add}.txt'
			)
		nfig=bz.plot_over(
			x_quant=trajx,
			nfig=nfig,
			file_name=f'byz{add}.png',
			plot=True,
			title=f'B$_y$ and B$_z$',
			dataFile=	f'/{res_folder_full}bz{add}.txt'
			)

		intBy.plot_over(
			x_quant=trajx,
			nfig=nfig,
			nosave=True,
			clear=True,
			dataFile=f'/{res_folder_full}intBy{add}.txt'
			)
		nfig=intBz.plot_over(
			x_quant=trajx,
			nfig=nfig,
			file_name=f'intByz{add}.png',
			plot=False,
			title=f'1. Integral of \n B$_y$ and B$_z$',
			dataFile=	f'/{res_folder_full}intBz{add}.txt'
			)

		intBy2.plot_over(
			x_quant=trajx,
			nfig=nfig,
			nosave=True,
			clear=True,
			dataFile=f'/{res_folder_full}int2By{add}.txt'
			)
		nfig=intBz2.plot_over(
			x_quant=trajx,
			nfig=nfig,
			file_name=f'int2Byz{add}.png',
			plot=False,
			title=f'2. Integral of \n B$_y$ and B$_z$',
			dataFile=	f'/{res_folder_full}int2Bz{add}.txt'
			)

		ByProfile.plot_over(
			x_quant=zProfile,
			nfig=nfig,
			nosave=True,
			clear=True,
			dataFile=f'/{res_folder_full}profileBy{add}.txt'
			)
		nfig=BzProfile.plot_over(
			x_quant=zProfile,
			nfig=nfig,
			file_name=f'profileBz{add}.png',
			plot=False,
			title=f'Profile of B$_y$ and B$_z$',
			dataFile=	f'/{res_folder_full}profileBz{add}.txt'
			)

		"""
		WAVE run
		"""
		pdb.set_trace()
		ebeam=uw.ebeam_parameters()
		ebeam.get_std_bessy_II_paras()

		unduBfield=uw.bfield.bfield(
			unitsXB=[0.001,1.0] # setting the units
			)

		unduBfield.load_field_from_file(
					file=dir_path+'/res/UNDU_DATA/undumag_on-axis.dat', 
					fieldMap=False,
					unduFile = True, 
					radiaFile=False,
					header=None,
					skiprows=None,
				)
		beff=unduBfield.calc_beff(prd_lngth=undulatorMagMag.getPeriodLength())
		
		anaUndu=uw.undulator.undulator(
			bEffY=beff,
			periodLength=undulatorMagMag.getPeriodLength()*0.001,
			ebeam=ebeam,
			)

		# aspec=uw.ana_eb.aspectrum(ebeam=ebeam,undulator=undu)
		# aspec.angular_power_distro_order(
		# 	r=6,
		# 	dX_scrn=0.04,
		# 	dY_scrn=0.04,
		# 	order=1,
		# 	num_integration_pnts=201,
		# 	omega=undu.properFrequencyStrong.get(),
		# 	sigmaOff=False,
		# 	piOff=False,
		# 	)
		# m0=undu.firstHarmEnergyStrong.get()
		# aspec.flux_factor_comparisson(order=1)
		# aspec.onAxisAngularSpectralFluxDensity(
		# 	ens=np.linspace(m0*0.5,m0*5.1,2000),
		# 	orders=[1,3,5],
		# 	freqs=None
		# 	)

		# order=1
		# X,Y,angularPowerDistro,avrg_total_power=aspec.angularSpectralFluxDistro_order(
		# 	r=6,
		# 	dX_scrn=0.04,
		# 	dY_scrn=0.04,
		# 	order=order,
		# 	num_integration_pnts=201,
		# 	omega=undu.properFrequencyStrong.get(),
		# 	sigmaOff=False,
		# 	piOff=False,
		# 	)

		# fig = plt.figure(figsize=(13*uw.cm_inch, 6.5*uw.cm_inch), dpi=150)
		# fig.suptitle(f"Power Distribution \n Avrg Pow: {avrg_total_power:2f} W", fontsize=14)
		# ax = fig.add_subplot(111, projection='3d')
		# ax.plot_surface(X, Y, angularPowerDistro, cmap=cm.coolwarm,
		#                        linewidth=0, antialiased=False)
		# ax.set_xlabel('x [m]', fontsize=12)
		# ax.set_ylabel('y [m]', fontsize=12)
		# ax.set_zlabel('Power [W]', fontsize=12)
		# plt.savefig(dir_path+f"/power_distro_order_{order}.png"   , bbox_inches='tight')
		# plt.show()

		# pdb.set_trace()

		# make the field known to wave

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
		wave_prog_paras.res_folder.set(res_folder_full)
		wave_prog_paras.calc_spectrum.set(True)
		wave_prog_paras.nthreads.set(6)

		"""
		Setting Spectrometer Parameter
		"""

		spectrometer_paras = wave._spectrometer_paras
		spectrometer_paras.spectrum_n_energies.set(100)
		spectrometer_paras.spectrum_min_energy.set(500)
		spectrometer_paras.spectrum_max_energy.set(1500)
		spectrometer_paras.spectrum_undu_mode.set(0)

		"""
		Setting Screen Parameter
		"""
		screen_paras = wave._screen_paras
		screen_paras.screen_segm_hor.set(30) 
		screen_paras.screen_segm_vert.set(30)
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

		# wave.run()

		"""
		Get Results and Plot
		"""

		results = wave.get_results()
		nfig=5

		traj_x = results.get_result(which='traj_x')
		traj_y = results.get_result(which='traj_y')
		traj_z = results.get_result(which='traj_z')

		traj_y.plot_over(x_quant=traj_x,nfig=nfig,nosave=True)
		nfig=traj_z.plot_over(x_quant=traj_x,nfig=nfig)

		By = results.get_result(which='By')
		By.plot_over(x_quant=traj_x,nfig=nfig,nosave=True)
		Bz = results.get_result(which='Bz')
		nfig=Bz.plot_over(x_quant=traj_x,nfig=nfig)
		pdb.set_trace()
		power_z = results.get_result(which='power_z')
		power_y = results.get_result(which='power_y')
		power_distro = results.get_result(which='power_distribution')
		nfig=power_distro.plot_over_3d(x_quant=power_y,y_quant=power_z,file_name=None,nosave=False,nfig=nfig)

		en_flux = results.get_result(which='en_flux')
		flux = results.get_result(which='flux')
		nfig=flux.plot_over(x_quant=en_flux,file_name=None,nosave=False,nfig=nfig,loglog=False)

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
		nfig=flux_density_onaxis.plot_over(x_quant=en_fd,nfig=nfig,loglog=False)

		flux_dens_distr_ens_loaded = results.find_load_flux_density_distribution(energies=[anaUndu.firstHarmEnergyStrong.get()])
		fd_y = results.get_result(which='fd_y')
		fd_z = results.get_result(which='fd_z')
		for en in flux_dens_distr_ens_loaded :
			fd = results.get_result(which=f'flux_density_distribution_{en:.2f}')
			nfig=fd.plot_over_3d(x_quant=fd_y,y_quant=fd_z,file_name=None,nosave=False,nfig=nfig)

		pdb.set_trace()

	def test_undulator_class_creation(self) :

		field_folder = f'/'
		res_folder = dir_path+'/res/'
		unduName='U49_2'

		ebeam=uw.ebeam_parameters()
		ebeam.get_std_bessy_II_paras()

		undu=uw.undulator.undulator(ebeam=ebeam)
		undu.set_paras_from_bessyII_undu_list(unduName=unduName)

		bfieldHarm=undu.getHarmonicBField(
			unitsXB=[1.0,1.0],
			wEnds=True
			)
		bfieldHarm.plot_fld(
			colx = 'x', 
			coly = 'By', 
			plt_extrm = False, 
			folder = dir_path+'/', 
			plot=False,
			save=True,
			filename=f'field_{unduName}.png'
			)

		# pdb.set_trace()
		# bfieldHarm.bvals

		# aspec=uw.ana_eb.aspectrum(ebeam=ebeam,undulator=undu)
		# # aspec.angular_power_distro_order(
		# # 	r=6,
		# # 	dX_scrn=0.04,
		# # 	dY_scrn=0.04,
		# # 	order=1,
		# # 	num_integration_pnts=201,
		# # 	omega=undu.properFrequencyStrong.get(),
		# # 	sigmaOff=False,
		# # 	piOff=False,
		# # 	)
		# # m0=undu.firstHarmEnergyStrong.get()
		# # aspec.flux_factor_comparisson(order=1)
		# # aspec.onAxisAngularSpectralFluxDensity(
		# # 	ens=np.linspace(m0*0.5,m0*5.1,2000),
		# # 	orders=[1,3,5],
		# # 	freqs=None
		# # 	)

		# order=1
		# X,Y,angularPowerDistro,avrg_total_power=aspec.angularSpectralFluxDistro_order(
		# 	r=6,
		# 	dX_scrn=0.04,
		# 	dY_scrn=0.04,
		# 	order=order,
		# 	num_integration_pnts=201,
		# 	omega=undu.properFrequencyStrong.get(),
		# 	sigmaOff=False,
		# 	piOff=False,
		# 	)

		# fig = plt.figure(figsize=(13*uw.cm_inch, 6.5*uw.cm_inch), dpi=150)
		# fig.suptitle(f"Power Distribution \n Avrg Pow: {avrg_total_power:2f} W", fontsize=14)
		# ax = fig.add_subplot(111, projection='3d')
		# ax.plot_surface(X, Y, angularPowerDistro, cmap=cm.coolwarm,
		#                        linewidth=0, antialiased=False)
		# ax.set_xlabel('x [m]', fontsize=12)
		# ax.set_ylabel('y [m]', fontsize=12)
		# ax.set_zlabel('Power [W]', fontsize=12)
		# plt.savefig(dir_path+f"/power_distro_order_{order}.png"   , bbox_inches='tight')
		# plt.show()

		# pdb.set_trace()

		# # make the field known to wave

		# """
		# Getting wave
		# """
		# wave = uw.wave(wave_mode='bfield')

		# bfield_paras = wave._bfield_paras # get bfield paras
		# bfield_paras.bfield.set(bfieldHarm) # set the bfield
		# bfield_paras.field_type.set('By') # tell wave which part of bfield to use for simu

		# """
		# Setting Program Parameter
		# """

		# wave_prog_paras = wave._prog_paras
		# wave_prog_paras.res_folder.set(res_folder)
		# wave_prog_paras.calc_spectrum.set(False)
		# wave_prog_paras.nthreads.set(6)

		# """
		# Setting Spectrometer Parameter
		# """

		# spectrometer_paras = wave._spectrometer_paras
		# spectrometer_paras.spectrum_n_energies.set(10)
		# spectrometer_paras.spectrum_min_energy.set(100)
		# spectrometer_paras.spectrum_max_energy.set(500)
		# spectrometer_paras.spectrum_undu_mode.set(0)

		# """
		# Setting Screen Parameter
		# """
		# screen_paras = wave._screen_paras
		# screen_paras.screen_segm_hor.set(30) 
		# screen_paras.screen_segm_vert.set(30)
		# screen_paras.screen_extent_hor.set(40) # pinhole width mm
		# screen_paras.screen_extent_vert.set(40) # pinhole height mm

		# """
		# Setting Undulator Parameter
		# """

		# undu_paras = wave._undu_paras # getting parameter object
		# undu_paras.elliptUnduB0Y.set(1.18)
		# undu_paras.elliptUnduB0Z.set(0.0)
		# undu_paras.elliptUnduNumPeriods.set(10)
		# undu_paras.elliptUnduPerLength.set(0.051)
		# undu_paras.elliptUnduPerShift.set(0.5)

		# """
		# Setting Beam Parameter
		# """

		# ebeam_paras = wave._ebeam_paras
		# ebeam_paras.beam_en.set(1.722) # [GeV]
		# ebeam_paras.current.set(0.3) # [A]
		# ebeam_paras.beamSizeHor.set(275e-6) # 
		# ebeam_paras.beamDiveHor.set(28.1e-6) #
		# ebeam_paras.beamSizeVer.set(22.5e-6) # 
		# ebeam_paras.beamDiveVer.set(6.8e-6) # 
		# ebeam_paras.espread.set(1e-3) # 
		# ebeam_paras.emittanceHor.set(7.7e-9)
		# ebeam_paras.emittanceVer.set(15.4e-11)
		# ebeam_paras.update_values()

		# """
		# Run
		# """

		# wave.run()

		# """
		# Get Results and Plot
		# """

		# results = wave.get_results()
		# nfig=0

		# traj_x = results.get_result(which='traj_x')
		# traj_y = results.get_result(which='traj_y')
		# traj_z = results.get_result(which='traj_z')

		# traj_y.plot_over(x_quant=traj_x,nfig=nfig,nosave=True)
		# nfig=traj_z.plot_over(x_quant=traj_x,nfig=nfig)

		# By = results.get_result(which='By')
		# By.plot_over(x_quant=traj_x,nfig=nfig,nosave=True)
		# Bz = results.get_result(which='Bz')
		# nfig=Bz.plot_over(x_quant=traj_x,nfig=nfig)
		# pdb.set_trace()
		# power_z = results.get_result(which='power_z')
		# power_y = results.get_result(which='power_y')
		# power_distro = results.get_result(which='power_distribution')
		# nfig=power_distro.plot_over_3d(x_quant=power_y,y_quant=power_z,file_name=None,nosave=False,nfig=nfig)

		# en_flux = results.get_result(which='en_flux')
		# flux = results.get_result(which='flux')
		# nfig=flux.plot_over(x_quant=en_flux,file_name=None,nosave=False,nfig=nfig,loglog=True)

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
		# nfig=flux_density_onaxis.plot_over(x_quant=en_fd,nfig=nfig,loglog=True)

		# flux_dens_distr_ens_loaded = results.find_load_flux_density_distribution(energies=[414])
		# fd_y = results.get_result(which='fd_y')
		# fd_z = results.get_result(which='fd_z')
		# for en in flux_dens_distr_ens_loaded :
		# 	fd = results.get_result(which=f'flux_density_distribution_{en:.2f}')
		# 	nfig=fd.plot_over_3d(x_quant=fd_y,y_quant=fd_z,file_name=None,nosave=False,nfig=nfig)

		# pdb.set_trace()

if __name__ == '__main__':
	unittest.main()
