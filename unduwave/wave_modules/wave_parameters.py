from unduwave.unduwave_incl import *
from unduwave.attribute_classes.attributes import _attribute
from unduwave.attribute_classes.attributes import _attribute_collection
import unduwave.helpers.file_folder_helpers as f_h
import unduwave.constants as uc
from unduwave.analytical_module import analytic_structures as anas

try :
	# works when calling script with python3 script_file
	dir_path = os.path.dirname(os.path.realpath(__file__))
except:
	# works when calling script with exec from python console
	dir_path = os.getcwd()

class ebeam_parameters(_attribute_collection):
	"""
	Defining basic electron-beam parameters
	beam_en - Beam energy in [GeV]
	current - current in [A]
	beamSizeHor - horizontal beam size [m]
	beamDiveHor - Horizontal beam divergence [rad]
	beamSizeVer - vertical beam size [m]
	beamDiveVer - vertical beam divergence [rad]
	espread - energy spread [%]
	emittanceHor/Ver - horizontal and vertical emittance [mrad]
	betaFunctionHor/Ver - horizontal and vertical beta functions [m]
	circumference - Ring circumference in [m]
	rdipol - Bending radius of dipoles [m]
	"""

	def __init__(self) :
		self.beam_en = _attribute(0,in_name='DMYENERGY')
		self.gammaFactor = _attribute(0)
		self.betaFactor = _attribute(0)
		self.current = _attribute(0,in_name='DMYCUR')
		self.beamSizeHor = _attribute(0,in_name='BSIGZ(1)') # horizontal beam size [m]
		self.beamDiveHor = _attribute(0,in_name='BSIGZP(1)') # hor beam divergence rad
		self.beamSizeVer = _attribute(0,in_name='BSIGY(1)') # Ver beam size m
		self.beamDiveVer = _attribute(0,in_name='BSIGYP(1)') # ver beam divergence rad
		self.espread = _attribute(0,in_name='ESPREAD') # ver beam divergence rad
		self.emittanceHor = _attribute(0,in_name='EPS0H')
		self.emittanceVer = _attribute(0,in_name='EPS0V')
		self.betaFunctionHor = _attribute(0,in_name='BETFUN')
		self.betaFunctionVer = _attribute(0,in_name='BETFUNV')
		self.circumference = _attribute(240,in_name='UMFANG')
		self.rdipol = _attribute(4.359,in_name='RDIPOL')
		super().__init__()

	def get_std_bessy_III_paras(self) :
		self.beam_en.set(2.5) # Beam energy in [GeV]
		self.current.set(0.3) # current in [A]
		self.beamSizeHor.set(275e-6) # horizontal beam size [m]
		self.beamDiveHor.set(28.1e-6) # Horizontal beam divergence [rad]
		self.beamSizeVer.set(22.5e-6) # vertical beam size [m]
		self.beamDiveVer.set(6.8e-6) # vertical beam divergence [rad]
		self.espread.set(0.963e-3) # energy spread [%]
		self.emittanceHor.set(9.804e-11) #  horizontal  emittance [mrad]
		self.emittanceVer.set(1.961e-12) #  vertical emittance [mrad]
		self.betaFunctionHor.set(0) #  horizontal beta functions [m]
		self.betaFunctionVer.set(0) #  vertical beta functions [m]
		self.circumference.set(350) # Ring circumference in [m]
		self.rdipol.set(2.78) # Bending radius of dipoles [m]
		self.update_values()
		return self

	def get_std_bessy_II_paras(self) :
		return self.get_std_paras()

	def get_std_paras(self): 
		self.beam_en.set(1.722) # Beam energy in [GeV]
		self.current.set(0.3) # current in [A]
		self.beamSizeHor.set(275e-6) # horizontal beam size [m]
		self.beamDiveHor.set(28.1e-6) # Horizontal beam divergence [rad]
		self.beamSizeVer.set(22.5e-6) # vertical beam size [m]
		self.beamDiveVer.set(6.8e-6) # vertical beam divergence [rad]
		self.espread.set(1e-3) # energy spread [%]
		self.emittanceHor.set(7.7e-9) #  horizontal  emittance [mrad]
		self.emittanceVer.set(15.4e-11) #  vertical emittance [mrad]
		self.betaFunctionHor.set(0) #  horizontal beta functions [m]
		self.betaFunctionVer.set(0) #  vertical beta functions [m]
		self.circumference.set(240) # Ring circumference in [m]
		self.rdipol.set(4.359) # Bending radius of dipoles [m]
		self.update_values()
		return self

	def update_values(self) : 
		self.gammaFactor.set( self.beam_en.get()*1e9*uc.q_el/(uc.m_el*uc.v_c**2) )
		self.betaFactor.set( math.sqrt( 1 - 1 / self.gammaFactor.get() ) )

class screen_parameters(_attribute_collection):
	"""
	Basic screen parameters
	screen_extent_hor/vert - width and height of pinhole [mm]
	screen_pos_x - distance of pinhole from center of undu [m]
	screen_segm_hor - number of points in z-direction
	screen_segm_vert - number of points in y-direction
	"""

	def __init__(self) :

		self.screen_extent_hor = _attribute(0,in_name='PINW',fac=1e-3)
		self.screen_extent_vert = _attribute(0,in_name='PINH',fac=1e-3)
		self.screen_pos_x = _attribute(0,in_name='PINCEN(1)')
		self.screen_pos_y = _attribute(9999.,in_name='PINCEN(2)')
		self.screen_pos_z = _attribute(9999.,in_name='PINCEN(3)')
		self.screen_segm_hor = _attribute(0,in_name='MPINZ')
		self.screen_segm_vert = _attribute(0,in_name='MPINY')
		super().__init__()

	def get_std_paras(self): 
		self.screen_extent_hor.set(3) # [mm]
		self.screen_extent_vert.set(3) # [mm]
		self.screen_pos_x.set(10) # [m]
		self.screen_segm_hor.set(10) # 
		self.screen_segm_vert.set(10) # 
		return self

class spectrometer_paras(_attribute_collection):
	"""
	Basic spectrometer parameters
	spectrum_min_energy, spectrum_max_energy - Energie at which to start/end spectrum calculation [eV]
	spectrum_n_energies - number of energies for which to calculate spectrum
	spectrum_undu_mode - undulator-mode (whole trajectory is source of radiation - coherent)
	spectrum_wigg_mode - wiggler-mode (only source-areas are considered and added incoherently)
	"""

	def __init__(self) :
		self.spectrum_min_energy = _attribute(0,in_name='FREQLOW')
		self.spectrum_max_energy = _attribute(0,in_name='FREQHIG')
		self.spectrum_n_energies = _attribute(0,in_name='NINTFREQ')
		self.spectrum_undu_mode = _attribute(0,in_name='IUNDULATOR')
		self.spectrum_wigg_mode = _attribute(0,in_name='IWIGGLER')
		self.spectrum_dipole_mode = _attribute(0,in_name='ISPECDIP')
		super().__init__()

	def get_std_paras(self): 
		self.spectrum_min_energy.set(300) # eV
		self.spectrum_max_energy.set(500) # eV
		self.spectrum_n_energies.set(5) # 
		self.spectrum_undu_mode.set(0) # undu-mode
		self.spectrum_wigg_mode.set(0) # wiggler-mode
		self.spectrum_dipole_mode.set(0) # dipole mode
		return self

class undu_paras(_attribute_collection):
	"""
	Parameters controlling the generation of the B-Field

	prog_parameters.undu_endp = 1
		planarUnduK - K-Parameter of Machine
		planarUnduB0 - B-Amplitude of Machine (either planarUnduK or this) [T]
		planarUnduPerLength- period length in x-direction [m]
		planarUnduNumPeriods - number of b-field peaks - each period contributes 2, the endpoles 3 -> odd number

	prog_parameters.undu_ellip = 1
		elliptUnduB0Y - B-Amplitude in y - [T]
		elliptUnduB0Z - B-Amplitude in z - [T]
		elliptUnduNumPeriods- numer of periods
		elliptUnduPerLength - period length - [m]
		elliptUnduPerShift - shift, % of period
	"""

	def __init__(self) :

		self.ebeam=_attribute(None)

		# Simple planar Undulator with endpoles
		self.planarUnduK = _attribute(0.0,in_name='PKHALBASY')
		self.planarUnduB0 = _attribute(1.0,in_name='B0HALBASY')
		self.planarUnduPerLength = _attribute(0.02,in_name='ZLHALBASY')
		self.planarUnduNumPeriods = _attribute(2,in_name='AHWPOL')
		self.undu_type = _attribute()

		# Simple elliptic undulator
		self.elliptUnduB0Y = _attribute(1.0,in_name='B0ELLIPV') # magn. field strength amplitude in vertical direction
		self.elliptUnduB0Z = _attribute(0.0,in_name='B0ELLIPH') # magn. field strength amplitude in horizontal direction
		self.elliptUnduNumPeriods = _attribute(2,in_name='PERELLIP') # number of periods
		self.elliptUnduPerLength = _attribute(0.02,in_name='XLELLIP') # period length
		self.elliptUnduPerShift = _attribute(0.0,in_name='ELLSHFT') # shift between the two magnetic arrays in fractions of one period

		self.bEffY = _attribute()
		self.bEffZ = _attribute()
		self.unduParameterKY = _attribute()
		self.unduParameterKZ = _attribute()
		self.shift = _attribute()
		self.periodLength = _attribute(None)
		self.numPeriods = _attribute(None)
		self.lengthEndPeriodsRelative = _attribute(None)
		super().__init__()

		# elliptic undulator with endpoles
		# b0ellana = _attribute(1.0,in_name='B0ELLANA') # Field Amplitude
		# nperella = _attribute(10,in_name='NPERELLA') # num of periods
		# xlellana = _attribute(0.02,in_name='XLELLANA') # period length in x in m
		# zlellana = _attribute(10,in_name='ZLELLANA') # period length in z ??? in m
		# x0ellana = _attribute(10,in_name='X0ELLANA') # x0 [m], distance of magnet center from device axis
		# gapell = _attribute(0.06,in_name='GAPELL') # gap [m]
		# refgapell = _attribute(0.06,in_name='REFGAPELL') # refrence gap [m] ?
		# shellana = _attribute(0.0,in_name='SHELLANA') # shift in units of zlellana
		# rowshella = _attribute(0.00,in_name='ROWSHELLA') # additional row shift of lower rows in units of zlellana
		# iells2s3 = _attribute(0,in_name='IELLS2S3') # >=0: S3-MODE - parallel; <0: S2-MODE - antiparallel
		# iellcoef = _attribute(0,in_name='IELLCOEF') # !>0: read IELLCOEF Fourier coefficients from file ellana.coef, =<0: First and second coefficients only with C0=0.5 and C1=1.

	def set_paras_from_bessyII_undu_list(self, unduName) :
		fullList=f_h.loadBessyIIundulatorList()
		indFnd=None
		for indEl, el in enumerate(fullList) : 
			if el['name'] == unduName :
				indFnd=indEl
				break
		if not (indFnd is None) :
			unduFnd=fullList[indFnd]
			if unduFnd['type'] == 'Apple II' : # elliptical undu
				self.undu_type.set("undu_ellip")
				self.bEffY.set(unduFnd['beffy [T]'])
				self.bEffZ.set(unduFnd['beffz [T]'])
				self.periodLength.set(unduFnd['period length [mm]']*1e-3)
				self.numPeriods.set(unduFnd['periods'])
			else : # planar (hybrid)
				self.undu_type.set("undu_endp")
				self.bEffY.set(unduFnd['beffy [T]'])
				self.periodLength.set(unduFnd['period length [mm]']*1e-3)
				self.numPeriods.set(unduFnd['periods'])
			self.update_values(thetaObservation=0.0)
		else:
			print(f"ana_undulator: set_paras_from_undu_list: Undulator {unduName} not found in list.")

	def get_std_paras(self,wave_mode,ebeam,thetaObservation=0.0): 
		"""
		getting standard undu parameters
		wave_mode - same as prog_parameters.wave_mode
		"""
		self.undu_type.set(wave_mode)

		self.bEffY.set(None)
		self.bEffZ.set(None)
		self.unduParameterKY.set(None)
		self.unduParameterKZ.set(None)
		self.shift.set(0.0)
		self.periodLength.set(0.02)
		self.numPeriods.set(78)
		self.lengthEndPeriodsRelative.set(1.5)
		self.ebeam.set(ebeam)
		self.update_values(thetaObservation=thetaObservation)
		return self

	def update_values(self,thetaObservation=0.0) :
		self._anasy=None
		self._anasz=None
		if not (self.unduParameterKY.get() is None) :
			self._anasy=anas.ana_undulator(
				bEff=None,
				unduK=self.unduParameterKY.get(),
				periodLength=self.periodLength.get(),
				numPeriods=self.numPeriods.get(),
				lengthEndPeriodsRelative=self.lengthEndPeriodsRelative.get(),
				ebeam=self.ebeam.get(),
				thetaObservation=thetaObservation
				)
		elif not (self.bEffY.get() is None) :
			self._anasy=anas.ana_undulator(
				bEff=self.bEffY.get(),
				periodLength=self.periodLength.get(),
				numPeriods=self.numPeriods.get(),
				lengthEndPeriodsRelative=self.lengthEndPeriodsRelative.get(),
				ebeam=self.ebeam.get(),
				thetaObservation=thetaObservation
				)
		if not (self.unduParameterKZ.get() is None) :
			self._anasz=anas.ana_undulator(
				bEff=None,
				unduK=self.unduParameterKZ.get(),
				periodLength=self.periodLength.get(),
				numPeriods=self.numPeriods.get(),
				lengthEndPeriodsRelative=self.lengthEndPeriodsRelative.get(),
				ebeam=self.ebeam.get(),
				thetaObservation=thetaObservation
				)
		elif not (self.bEffZ.get() is None) :
			self._anasz=anas.ana_undulator(
				bEff=self.bEffZ.get(),
				periodLength=self.periodLength.get(),
				numPeriods=self.numPeriods.get(),
				lengthEndPeriodsRelative=self.lengthEndPeriodsRelative.get(),
				ebeam=self.ebeam.get(),
				thetaObservation=thetaObservation
				)
		if self.undu_type.get() == 'undu_ellip' :
			shift=self.shift.get()
			bEffY=self.bEffY.get()
			if (bEffY is None) and (not (self._anasy is None)) :
				bEffY=self._anasy.bEff.get()
			bEffZ=self.bEffZ.get()
			if (bEffZ is None) and (not (self._anasz is None)) :
				bEffZ=self._anasz.bEff.get()
			if shift == 0.0 :
				self.elliptUnduB0Y.set(bEffY)
				self.elliptUnduB0Z.set(0.0)
			elif shift==0.25 :
				self.elliptUnduB0Y.set(bEffY)
				self.elliptUnduB0Z.set(bEffZ)
			elif shift==-0.5 : # antiparallel in this crazy world
				self.shift.set(0.0)
				self.elliptUnduB0Y.set(bEffY)
				self.elliptUnduB0Z.set(bEffY)
			elif shift==0.5 :
				self.elliptUnduB0Y.set(0.0)
				self.elliptUnduB0Z.set(bEffZ)
			elif shift==0.75:
				self.elliptUnduB0Y.set(bEffY)
				self.elliptUnduB0Z.set(bEffZ)
			self.elliptUnduPerShift.set(self.shift.get())
			self.elliptUnduNumPeriods.set(self.numPeriods.get()) 
			self.elliptUnduPerLength.set(self.periodLength.get())
		elif self.undu_type.get() == 'undu_endp' :
			if not (self.bEffY.get() is None) :
				self.planarUnduB0.set(self.bEffY.get())
			elif not (self.unduParameterKY.get() is None) :
				self.planarUnduK.set(self.unduParameterKY.get())
			self.planarUnduNumPeriods.set(2*self.numPeriods.get()+1) # we count the number of B-field peaks here - one extra for the end-fields (odd)
			self.planarUnduPerLength.set(self.periodLength.get())

class bfield_paras(_attribute_collection):

	def __init__(self) :
		self.field_folder = _attribute('/')
		self.bfield=_attribute()
		self.field_type=_attribute()
		"""
		Possible Types: 
		By / Byz / Bxyz / Bmap
		"""
		super().__init__()

	def get_std_paras(self,bfield=None): 	
		self.field_folder.set("/")
		self.bfield.set(bfield)
		self.field_type.set('By')

class wave_prog_parameters(_attribute_collection):
	"""
	Represents standard parameters for wave simulations.

	"""

	def __init__(self) :
		self.wave_prog_folder = _attribute('')
		self.in_file_folder = _attribute('')
		self.in_files = _attribute({})
		self.field_folder = _attribute('')
		self.field_files = _attribute([])

		self.four_file=_attribute('')
		self.res_folder = _attribute('')
		self.wave_data_res_folder = _attribute('')
		self.pics_folder = _attribute('')
		self.res_summary_file = _attribute('')
		self.no_copy = _attribute([])
		self.wave_ending_extract = _attribute([])
		self.wave_ending_copy = _attribute([])
		self.wave_files_essentials = _attribute([])
		self.wave_res_copy_behaviour = _attribute()
		self.zip_res_folder = _attribute(1)
		self.nthreads = _attribute(2,in_name='MTHREADS')
		self.zipped = _attribute(True)
		self.calc_spectrum = _attribute(False,in_name='ISPEC')
		self.calc_emittance = _attribute(1,in_name='IFOLD')
		self.calc_energy_fold = _attribute(1,in_name='IEFOLD')
		self.emittance_fold_with_sigmas = _attribute(1,in_name='ISIGUSR')
		self.source_point_calc = _attribute(9999.,in_name='WGWINFC')
		self.ihisascii = _attribute(111,in_name='IHISASCII')
		self.ntupgrid = _attribute(0,in_name='NTUPGRID')
		self.rayfile = _attribute(0,in_name='IWFILRAY')
		self.electron_intermediate_x = _attribute(-9999.,in_name='XINTER')
		self.electron_x0 = _attribute(9999.,in_name='XSTART')
		self.electron_y0 = _attribute(0.0,in_name='YSTART')
		self.electron_z0 = _attribute(0.0,in_name='ZSTART')
		self.electron_end_x = _attribute(9999.,in_name='XSTOP')
		self.electron_vx0 = _attribute(1.0,in_name='VXIN')
		self.electron_vy0 = _attribute(0.0,in_name='VYIN')
		self.electron_vz0 = _attribute(0.0,in_name='VZIN')

		"""
		Writing Field-Maps
		"""
		self.bxmapmin = _attribute(9999.,in_name='XMAPMN')
		self.bxmapmax = _attribute(9999.,in_name='XMAPMX')
		self.bxmapn = _attribute(-9999,in_name='NMAPX')
		self.bymapmin = _attribute(9999.,in_name='YMAPMN')
		self.bymapmax = _attribute(9999.,in_name='YMAPMX')
		self.bymapn = _attribute(-9999,in_name='NMAPY')
		self.bzmapmin = _attribute(9999.,in_name='ZMAPMN')
		self.bzmapmax = _attribute(9999.,in_name='ZMAPMX')
		self.bzmapn = _attribute(-9999,in_name='NMAPZ')

		# undu_type = _attribute('')

		self.b_type = _attribute('none')
		self.irbtab = _attribute(0,in_name='IRBTAB')
		self.irfileb0 = _attribute(0,in_name='IRFILB0')
		self.iwfilf = _attribute(0,in_name='IWFILF')
		self.irfilf = _attribute(0,in_name='IRFILF')
		self.nfourwls = _attribute(32,in_name='NFOURWLS')
		self.ifour0 = _attribute(0,in_name='ifour0')
		self.iwbmap = _attribute(0,in_name='IWBMAP')
		self.irbtabzy = _attribute(0,in_name='IRBTABZY')
		self.irbtabxyz = _attribute(0,in_name='IRBTABXYZ')
		self.undu_easy = _attribute(0,in_name='KHALBA') # magnetic structure without specific ends
		self.undu_endp = _attribute(0,in_name='KHALBASY') # magnetic structure with endpoles
		self.undu_gap = _attribute(0,in_name='KUNDUGAP') # undulator with analytic gap-variation
		self.undu_ellip = _attribute(0,in_name='KELLIP') # elliptic undulator
		self.undu_ellip_ana = _attribute(0,in_name='KELLANA') # elliptic undulator
		self.wave_mode = _attribute('undu_easy') # wave mode
		super().__init__()

	def get_std_paras(self, wave_mode = 'undu_easy' ): 
		"""
		Getting standard wave-parameters depending on mode
		- wave_mode = 'By' - takes by field data and runs with that
		- wave_mode = 'Byz'
		- wave_mode = 'Bxyz'
		- wave_mode = 'undu_ellip' - standard elliptical undulator
		- wave_mode = 'undu_easy'
		- wave_mode = 'undu_endp'
		- wave_mode = 'undu_gap'
		"""

		#std-bad!
		dir_path = os.path.dirname(os.path.realpath(__file__))
		self.wave_prog_folder.set(dir_path+f_h.convert_path_to_win('/../../External-Software/WAVE/'))
		self.in_file_folder.set(dir_path+f_h.convert_path_to_win('/../UNDWAVE_IN_FILES/WAVE-In-Files/'))
		#for containerization?
		# dir_path='/var/lib/unduwave'
		# self.wave_prog_folder.set(dir_path+f_h.convert_path_to_win('/External-Software/WAVE/'))
		# self.in_file_folder.set(dir_path+f_h.convert_path_to_win('/unduwave/UNDWAVE_IN_FILES/WAVE-In-Files/'))

		self.in_files.set({ 'By' : 'load_ext_on_axis_by_ALL_OUT.in', 
							'Byz' : 'load_ext_on_axis_byz_ALL_OUT.in', 
							'Bxyz' : 'load_ext_on_axis_bxyz_ALL_OUT.in', 
							'undu_ellip' : 'wave.in',
							'none' : 'wave.in',
							'bmap' : 'wave.in',
							})
		self.field_folder.set('')
		self.four_file.set('')
		self.field_files.set([])
		self.res_folder.set('')
		self.wave_data_res_folder.set('WAVE_DATA/')
		self.pics_folder.set('Pics/')
		self.res_summary_file.set('res_summary.txt')
		self.no_copy.set(['WAVE_CODE.DAT', 'undumag_mu_77K.dat', 'undumag_mu_300K.dat',
								'iron_muinf_sat-2.34.dat', 'Vanadium_Permendur_Radia', #'WAVE.mhb' 
								])
		self.wave_ending_extract.set([ 'dat', 'wva', 'wvh', 'out' ])
		self.wave_ending_copy.set([ 'in' ])
		self.wave_files_essentials.set([ 
				'stokes_dist_emittance_espread', 'trajectory',
				'irradiated_power_dist', 'brilliance_3702', 
				'photon_flux_(pinhole)_48000','wave_ray.dat', 'bmap.dat', 'btab.fou',
				'selected_s0_e_(folded)_x_1_e_6_180000', 'wave.out','WAVE.mhb' ])
		self.nthreads.set(2)
		self.iwfilf.set(0)
		self.irfilf.set(0)
		self.nfourwls.set(32)
		self.ifour0.set(0)
		self.zipped.set(False)
		self.wave_res_copy_behaviour.set('copy_essentials')
		self.zip_res_folder.set(0)
		self.calc_emittance.set(1)
		self.calc_energy_fold.set(1)
		self.emittance_fold_with_sigmas.set(1)
		self.calc_spectrum.set(1) # boolean
		self.wave_mode.set(wave_mode)

		self.b_type.set('none')
		self.irbtab.set(0)
		self.irbtabzy.set(0)
		self.irbtabxyz.set(0)
		self.undu_easy.set(0)
		self.undu_endp.set(0)
		self.undu_gap.set(0)
		self.undu_ellip.set(0)

	def update_values(self,bfield_paras=None) :

		wave_mode=self.wave_mode.get()
		if wave_mode == 'bfield' :
			if not bfield_paras is None :
				field_type=bfield_paras.field_type.get()
			else :
				return
			if field_type == 'By' :
				self.irbtab.set(-2)
			elif field_type == 'Byz' :
				self.irbtabzy.set(1)
			elif field_type == 'Bxyz' :
				self.irbtabxyz.set(1)
			elif field_type == 'bmap' :
				self.ntupgrid.set(-1)
				# self.irfileb0.set(-6) # loading field map
				self.irfileb0.set(6) # loading field map
		else :
			if wave_mode == 'undu_ellip' :
				self.undu_ellip.set(1)
			elif wave_mode == 'undu_easy' :
				self.undu_easy.set(1)
			elif wave_mode == 'undu_endp' :
				self.undu_endp.set(1)
			elif wave_mode == 'undu_gap' :
				self.undu_gap.set(1)
			elif wave_mode == 'undu_ellip_ana' :
				self.undu_ellip_ana.set(1)
			elif wave_mode == 'four' :
				self.iwfilf.set(0)
				self.irfilf.set(1)

		return self

