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

class btests(unittest.TestCase) :

	def test_load_bfields(self) :

		field_folder = f'/'
		res_folder = dir_path+'/res/'

		"""
		Loading and setting a BField
		"""

		# Loading a file with x [mm] and By [T]

		mmT_bfield=uw.bfield.bfield(
			unitsXB=[0.001,1.0] # setting the units
			)

		mmT_bfield.load_field_from_file(
					file=dir_path+'/field_by_mm_T.dat', 
					fieldMap=False,
					cols=['x','By'],
					unduFile = False, 
					radiaFile=False,
					header=None,
					skiprows=None,
				)
		bfieldVals=mmT_bfield.bvals[454,:]
		bfieldVals=[round(val,6) for val in bfieldVals]
		self.assertEqual( bfieldVals,[0,round(-0.3842,6),0] )

		unduBfield=uw.bfield.bfield(
			unitsXB=[0.001,1.0] # setting the units
			)

		unduBfield.load_field_from_file(
					file=dir_path+'/undu_field.dat', 
					fieldMap=False,
					unduFile = True, 
					radiaFile=False,
					header=None,
					skiprows=None,
				)
		bfieldVals=unduBfield.bvals[405,:]
		bfieldVals=[round(val,6) for val in bfieldVals]
		self.assertEqual( bfieldVals,[0,round(6.0548586E-002,6),0] )

		# Loading a file with x [m] and By [kT]

		bfieldUE51=uw.bfield.bfield(
			unitsXB=[0.001,1.0] # setting the units
			)
		bfieldUE51.load_field_from_file(
					file=dir_path+'/UE51_field.dat', 
					fieldMap=False,
					cols=['x','By'],
					unduFile = True, 
					radiaFile=False,
					header=None,
					skiprows=None,
				)
		beff=bfieldUE51.calc_beff(prd_lngth=51.3, n_max = 10,colx = 'x')

		"""
		UE-51 Test
		"""

		mkT_bfield=uw.bfield.bfield(
			unitsXB=[1.0,1000.0] # setting the units
			)

		mkT_bfield.load_field_from_file(
					file=dir_path+'/field_by_m_kT.dat', 
					fieldMap=False,
					cols=['x','By'],
					unduFile = False, 
					radiaFile=False,
					header=None,
					skiprows=None,
				)
		bfieldVals=mkT_bfield.bvals[192,:]
		bfieldVals=[round(val,6) for val in bfieldVals]
		self.assertEqual( bfieldVals,[0,round(-0.0001196,6),0] )
		self.assertEqual( mkT_bfield.xvals[192],-0.738 )

		map_bfield=uw.bfield.bfield(
			unitsXB=[1.0,1.0] # setting the units
			)

		map_bfield.load_field_from_file(
					file=dir_path+'/field_map.map', 
					fieldMap=True,
					cols=['x','y','z','Bx','By','Bz'],
					unduFile = False, 
					radiaFile=False,
					header=None,
					skiprows=range(0,6),
				)
		bfieldVals=map_bfield.bvals[1,0,0,:]
		bfieldVals=[round(val,6) for val in bfieldVals]
		self.assertEqual( bfieldVals,[0,round(-3.647990759127e-03,6),0] )

		map_bfield.write_field_map_wave(file=dir_path+'/testMap.dat')
		"""
		add a read testMap.dat and then compare expectation to actual numbers
		"""

		"""
		Test the ana field creation 
		"""
		bfieldHarm=uw.bfield.bfield(
			unitsXB=[1.0,1.0] # setting the units
			)
		bfieldHarm=bfieldHarm.create_harm_field(
			periodLength=0.02, 
			amplitude=1.3, 
			nperiods=2, 
			phase_shift = 0, 
			num_pnts_per_period = 100, 
			colx = 'x', 
			coly = 'By',
			deltaX=None,
			)

		bfieldHarm2=uw.bfield.bfield(
			unitsXB=[1.0,1.0] # setting the units
			)
		bfieldHarm2=bfieldHarm2.create_harm_field(
			periodLength=0.02, 
			amplitude=1.3, 
			nperiods=2, 
			phase_shift = 0, 
			num_pnts_per_period = 100, 
			colx = 'x', 
			coly = 'By',
			deltaX=None,
			)

		addedFld=bfieldHarm.add_field(
			bfieldA=bfieldHarm2, 
			pos = 'front',#'front', 'back', None
			colx='x',
			pntsOut=None,
			)

		# addedFld.plot_fld( 
		# 	colx = 'x', 
		# 	coly = 'By', 
		# 	plt_extrm = False, 
		# 	folder = dir_path+'/', 
		# 	plot=True,
		# 	save=True,
		# 	filename='add.png'
		# 	) 

		endFld=uw.bfield.bfield(
			unitsXB=[1.0,1.0] # setting the units
			)
		endFld=endFld.create_harm_field_with_ends(
				periodLength=0.0494, 
				amplitude=0.788, 
				nperiods=83,
				num_pnts_per_period = 200, 
				colx = 'x', 
				coly = 'By',
				shift=0,
				)
		endFld.plot_fld(
			colx = 'x', 
			coly = 'By', 
			plt_extrm = False, 
			folder = dir_path+'/', 
			plot=False,
			save=True,
			filename='endFld.png'
			)
		beff=endFld.calc_beff(prd_lngth=0.0494, n_max = 10,colx = 'x')
		# pdb.set_trace()
		# make the field known to wave

		"""
		Getting wave
		"""
		# return
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
