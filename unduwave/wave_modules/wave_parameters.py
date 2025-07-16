from unduwave.unduwave_incl import *
from unduwave.attribute_classes.attributes import _attribute
from unduwave.attribute_classes.attributes import _attribute_collection
import unduwave.helpers.file_folder_helpers as f_h

class ebeam_parameters(_attribute_collection):
	"""
	Defining basic electron-beam parameters
	beam_en - Beam energy in [GeV]
	current - current in [A]
	bsigz - horizontal beam size [m]
	bsigzp - Horizontal beam divergence [rad]
	bsigy - vertical beam size [m]
	bsigyp - vertical beam divergence [rad]
	espread - energy spread [%]
	emitt_h/v - horizontal and vertical emittance [mrad]
	betfunh/v - horizontal and vertical beta functions [m]
	circumference - Ring circumference in [m]
	rdipol - Bending radius of dipoles [m]
	"""
	beam_en = _attribute(0,in_name='DMYENERGY')
	current = _attribute(0,in_name='DMYCUR')
	bsigz = _attribute(0,in_name='BSIGZ(1)') # horizontal beam size [m]
	bsigzp = _attribute(0,in_name='BSIGZP(1)') # hor beam divergence rad
	bsigy = _attribute(0,in_name='BSIGY(1)') # Ver beam size m
	bsigyp = _attribute(0,in_name='BSIGYP(1)') # ver beam divergence rad
	espread = _attribute(0,in_name='ESPREAD') # ver beam divergence rad
	emitt_h = _attribute(0,in_name='EPS0H')
	emitt_v = _attribute(0,in_name='EPS0V')
	betfunh = _attribute(0,in_name='BETFUN')
	betfunv = _attribute(0,in_name='BETFUNV')
	circumference = _attribute(240,in_name='UMFANG')
	rdipol = _attribute(4.359,in_name='RDIPOL')

	def get_std_bessy_III_paras(self) :
		self.beam_en.set(2.5) # Beam energy in [GeV]
		self.current.set(0.3) # current in [A]
		self.bsigz.set(275e-6) # horizontal beam size [m]
		self.bsigzp.set(28.1e-6) # Horizontal beam divergence [rad]
		self.bsigy.set(22.5e-6) # vertical beam size [m]
		self.bsigyp.set(6.8e-6) # vertical beam divergence [rad]
		self.espread.set(0.963e-3) # energy spread [%]
		self.emitt_h.set(9.804e-11) #  horizontal  emittance [mrad]
		self.emitt_v.set(1.961e-12) #  vertical emittance [mrad]
		self.betfunh.set(0) #  horizontal beta functions [m]
		self.betfunv.set(0) #  vertical beta functions [m]
		self.circumference.set(350) # Ring circumference in [m]
		self.rdipol.set(2.78) # Bending radius of dipoles [m]
		return self

	def get_std_bessy_II_paras(self) :
		return self.get_std_paras()

	def get_std_paras(self): 
		self.beam_en.set(1.722) # Beam energy in [GeV]
		self.current.set(0.3) # current in [A]
		self.bsigz.set(275e-6) # horizontal beam size [m]
		self.bsigzp.set(28.1e-6) # Horizontal beam divergence [rad]
		self.bsigy.set(22.5e-6) # vertical beam size [m]
		self.bsigyp.set(6.8e-6) # vertical beam divergence [rad]
		self.espread.set(1e-3) # energy spread [%]
		self.emitt_h.set(7.7e-9) #  horizontal  emittance [mrad]
		self.emitt_v.set(15.4e-11) #  vertical emittance [mrad]
		self.betfunh.set(0) #  horizontal beta functions [m]
		self.betfunv.set(0) #  vertical beta functions [m]
		self.circumference.set(240) # Ring circumference in [m]
		self.rdipol.set(4.359) # Bending radius of dipoles [m]

		return self

class screen_parameters(_attribute_collection):
	"""
	Basic screen parameters
	pinh_w/h - width and height of pinhole [mm]
	pinh_x - distance of pinhole from center of undu [m]
	pinh_nz - number of points in z-direction
	pinh_ny - number of points in y-direction
	"""
	pinh_w = _attribute(0,in_name='PINW',fac=1e-3)
	pinh_h = _attribute(0,in_name='PINH',fac=1e-3)
	pinh_x = _attribute(0,in_name='PINCEN(1)')
	pinh_nz = _attribute(0,in_name='MPINZ')
	pinh_ny = _attribute(0,in_name='MPINY')

	def get_std_paras(self): 
		self.pinh_w.set(3) # [mm]
		self.pinh_h.set(3) # [mm]
		self.pinh_x.set(10) # [m]
		self.pinh_nz.set(10) # 
		self.pinh_ny.set(10) # 
		return self

class spectrometer_paras(_attribute_collection):
	"""
	Basic spectrometer parameters
	freq_low/high - Energie at which to start/end spectrum calculation [eV]
	freq_num - number of energies for which to calculate spectrum
	undu - undulator-mode (whole trajectory is source of radiation - coherent)
	wigg - wiggler-mode (only source-areas are considered and added incoherently)
	"""
	freq_low = _attribute(0,in_name='FREQLOW')
	freq_high = _attribute(0,in_name='FREQHIG')
	freq_num = _attribute(0,in_name='NINTFREQ')
	undu = _attribute(0,in_name='IUNDULATOR')
	wigg = _attribute(0,in_name='IWIGGLER')

	def get_std_paras(self): 
		self.freq_low.set(300) # eV
		self.freq_high.set(500) # eV
		self.freq_num.set(5) # 
		self.undu.set(0) # undu-mode
		self.wigg.set(0) # wiggler-mode
		return self

class undu_paras(_attribute_collection):
	"""
	Parameters controlling the generation of the B-Field

	prog_parameters.undu_endp = 1
		pkHalbasy - K-Parameter of Machine
		b0Halbasy - B-Amplitude of Machine (either pkHalbasy or this) [T]
		xlHalbasy- period length in x-direction [m]
		ahwpolHalbasy - number of main poles (odd number)

	prog_parameters.undu_ellip = 1
		b0y - B-Amplitude in y - [T]
		b0z - B-Amplitude in z - [T]
		nper- numer of periods
		perl_x - period length - [m]
		ell_shift - shift, % of period
	"""
	pkHalbasy = _attribute(0.0,in_name='PKHALBASY')
	b0Halbasy = _attribute(1.0,in_name='B0HALBASY')
	xlHalbasy = _attribute(0.018,in_name='ZLHALBASY')
	ahwpolHalbasy = _attribute(5,in_name='AHWPOL')
	undu_type = _attribute("undu_ellip")

	b0y = _attribute(0,in_name='B0ELLIPV') # magn. field strength amplitude in vertical direction
	b0z = _attribute(0,in_name='B0ELLIPH') # magn. field strength amplitude in horizontal direction
	nper = _attribute(0,in_name='PERELLIP') # number of periods
	perl_x = _attribute(0,in_name='XLELLIP') # period length
	ell_shift = _attribute(0.25,in_name='ELLSHFT') # shift between the two magnetic arrays in fractions of one period

	b0ellana = _attribute(1.0,in_name='B0ELLANA') # Field Amplitude
	nperella = _attribute(10,in_name='NPERELLA') # num of periods
	xlellana = _attribute(0.02,in_name='XLELLANA') # period length in x in m
	zlellana = _attribute(10,in_name='ZLELLANA') # period length in z ??? in m
	x0ellana = _attribute(10,in_name='X0ELLANA') # x0 [m], distance of magnet center from device axis
	gapell = _attribute(0.06,in_name='GAPELL') # gap [m]
	refgapell = _attribute(0.06,in_name='REFGAPELL') # refrence gap [m] ?
	shellana = _attribute(0.0,in_name='SHELLANA') # shift in units of zlellana
	rowshella = _attribute(0.00,in_name='ROWSHELLA') # additional row shift of lower rows in units of zlellana
	iells2s3 = _attribute(0,in_name='IELLS2S3') # >=0: S3-MODE - parallel; <0: S2-MODE - antiparallel
	iellcoef = _attribute(0,in_name='IELLCOEF') # !>0: read IELLCOEF Fourier coefficients from file ellana.coef, =<0: First and second coefficients only with C0=0.5 and C1=1.

	def get_std_paras(self,wave_mode): 
		"""
		getting standard undu parameters
		wave_mode - same as prog_parameters.wave_mode
		"""
		self.undu_type.set(wave_mode)
		if self.undu_type.get() == 'undu_ellip' :
			self.b0y.set(0.3)
			self.b0z.set(1.0)
			self.nper.set(5)
			self.perl_x.set(0.02)
			self.ell_shift.set(0.25)
		return self

class bfield_paras(_attribute_collection):
	field_folder = _attribute('/')

	def get_std_paras(self): 	
		self.field_folder.set("/")

class wave_prog_parameters(_attribute_collection):
	"""
	Represents standard parameters for wave simulations.

	"""

	wave_prog_folder = _attribute('')
	in_file_folder = _attribute('')
	in_files = _attribute({})
	field_folder = _attribute('')
	field_files = _attribute([])

	res_folder = _attribute('')
	wave_data_res_folder = _attribute('')
	pics_folder = _attribute('')
	res_summary_file = _attribute('')
	no_copy = _attribute([])
	wave_ending_extract = _attribute([])
	wave_ending_copy = _attribute([])
	wave_files_essentials = _attribute([])
	wave_res_copy_behaviour = _attribute()
	zip_res_folder = _attribute(1)
	nthreads = _attribute(2,in_name='MTHREADS')
	zipped = _attribute(True)
	spec_calc = _attribute(False,in_name='ISPEC')
	iemit = _attribute(0,in_name='IEMIT')
	iefold = _attribute(1,in_name='IEFOLD')
	isigusr = _attribute(1,in_name='ISIGUSR')
	ihisascii = _attribute(111,in_name='IHISASCII')
	ntupgrid = _attribute(0,in_name='NTUPGRID')
	rayfile = _attribute(0,in_name='IWFILRAY')
	xinter = _attribute(-9999.,in_name='XINTER')
	xstart = _attribute(9999.,in_name='XSTART')
	ystart = _attribute(0.0,in_name='YSTART')
	zstart = _attribute(0.0,in_name='ZSTART')
	xstop = _attribute(9999.,in_name='XSTOP')
	vxin = _attribute(1.0,in_name='VXIN')
	vyin = _attribute(0.0,in_name='VYIN')
	vzin = _attribute(0.0,in_name='VZIN')
	"""
	a comment
	"""

	bxmapmin = _attribute(9999.,in_name='XMAPMN')
	bxmapmax = _attribute(9999.,in_name='XMAPMX')
	bxmapn = _attribute(-9999,in_name='NMAPX')
	bymapmin = _attribute(9999.,in_name='YMAPMN')
	bymapmax = _attribute(9999.,in_name='YMAPMX')
	bymapn = _attribute(-9999,in_name='NMAPY')
	bzmapmin = _attribute(9999.,in_name='ZMAPMN')
	bzmapmax = _attribute(9999.,in_name='ZMAPMX')
	bzmapn = _attribute(-9999,in_name='NMAPZ')

	# undu_type = _attribute('')

	b_type = _attribute('none')
	irbtab = _attribute(0,in_name='IRBTAB')
	irfileb0 = _attribute(0,in_name='IRFILB0')
	iwbmap = _attribute(0,in_name='IWBMAP')	
	irbtabzy = _attribute(0,in_name='IRBTABZY')
	irbtabxyz = _attribute(0,in_name='IRBTABXYZ')
	undu_easy = _attribute(0,in_name='KHALBA') # magnetic structure without specific ends
	undu_endp = _attribute(0,in_name='KHALBASY') # magnetic structure with endpoles
	undu_gap = _attribute(0,in_name='KUNDUGAP') # undulator with analytic gap-variation
	undu_ellip = _attribute(0,in_name='KELLIP') # elliptic undulator
	undu_ellip_ana = _attribute(0,in_name='KELLANA') # elliptic undulator
	wave_mode = _attribute('undu_easy') # wave mode

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
				'photon_flux_(pinhole)_48000','wave_ray.dat', 'bmap.dat',
				'selected_s0_e_(folded)_x_1_e_6_180000', 'wave.out','WAVE.mhb' ])
		self.nthreads.set(2)
		self.zipped.set(False)
		self.wave_res_copy_behaviour.set('copy_essentials')
		self.zip_res_folder.set(0)
		self.iemit.set(0)
		self.iefold.set(1)
		self.isigusr.set(1)
		self.spec_calc.set(1) # boolean
		self.wave_mode.set(wave_mode)

		self.b_type.set('none')
		self.irbtab.set(0)
		self.irbtabzy.set(0)
		self.irbtabxyz.set(0)
		self.undu_easy.set(0)
		self.undu_endp.set(0)
		self.undu_gap.set(0)
		self.undu_ellip.set(0)

		if wave_mode == 'By' :
			self.b_type.set('By')
			self.irbtab.set(-2)
		elif wave_mode == 'Byz' :
			self.b_type.set('Byz')
			self.irbtabzy.set(1)
		elif wave_mode == 'Bxyz' :
			self.b_type.set('Bxyz')
			self.irbtabxyz.set(1)
		elif wave_mode == 'undu_ellip' :
			self.undu_ellip.set(1)
		elif wave_mode == 'undu_easy' :
			self.undu_easy.set(1)
		elif wave_mode == 'undu_endp' :
			self.undu_endp.set(1)
		elif wave_mode == 'undu_gap' :
			self.undu_gap.set(1)
		elif wave_mode == 'undu_ellip_ana' :
			self.undu_ellip_ana.set(1)
		elif wave_mode == 'bmap' :
			self.b_type.set('bmap')
			self.ntupgrid.set(-1)
			self.irfileb0.set(-6) # loading field map

		return self

