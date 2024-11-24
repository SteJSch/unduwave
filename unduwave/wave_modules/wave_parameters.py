from unduwave.unduwave_incl import *
from unduwave.attribute_classes.attributes import _attribute
from unduwave.attribute_classes.attributes import _attribute_collection

class ebeam_parameters(_attribute_collection):
	beam_en = _attribute(0,wave_in_name='DMYENERGY')
	current = _attribute(0,wave_in_name='DMYCUR')
	bsigz = _attribute(0,wave_in_name='BSIGZ(1)') # horizontal beam size [m]
	bsigzp = _attribute(0,wave_in_name='BSIGZP(1)') # hor beam divergence rad
	bsigy = _attribute(0,wave_in_name='BSIGY(1)') # Ver beam size m
	bsigyp = _attribute(0,wave_in_name='BSIGYP(1)') # ver beam divergence rad
	espread = _attribute(0,wave_in_name='ESPREAD') # ver beam divergence rad
	emitt_h = _attribute(0,wave_in_name='EPS0H')
	emitt_v = _attribute(0,wave_in_name='EPS0V')
	betfunh = _attribute(0,wave_in_name='BETFUN')
	betfunv = _attribute(0,wave_in_name='BETFUNV')
	circumference = _attribute(240,wave_in_name='UMFANG')
	rdipol = _attribute(4.359,wave_in_name='RDIPOL')

	def get_std_paras(self): 
		self.beam_en.set(1.722) # [GeV]
		self.current.set(0.3) # [A]

		self.bsigz.set(275e-6) # 
		self.bsigzp.set(28.1e-6) #
		self.bsigy.set(22.5e-6) # 
		self.bsigyp.set(6.8e-6) # 
		self.espread.set(1e-3) # 
		self.emitt_h.set(7.7e-9)
		self.emitt_v.set(15.4e-11)
		self.betfunh.set(0)
		self.betfunv.set(0)
		self.circumference.set(240)
		self.rdipol.set(4.359)
		return self

class screen_parameters(_attribute_collection):
	pinh_w = _attribute(0,wave_in_name='PINW',fac=1e-3)
	pinh_h = _attribute(0,wave_in_name='PINH',fac=1e-3)
	pinh_x = _attribute(0,wave_in_name='PINCEN(1)')
	pinh_nz = _attribute(0,wave_in_name='MPINZ')
	pinh_ny = _attribute(0,wave_in_name='MPINY')

	def get_std_paras(self): 
		self.pinh_w.set(3) # [mm]
		self.pinh_h.set(3) # [mm]
		self.pinh_x.set(10) # [m]
		self.pinh_nz.set(10) # 
		self.pinh_ny.set(10) # 
		return self

class spectrometer_paras(_attribute_collection):
	freq_low = _attribute(0,wave_in_name='FREQLOW')
	freq_high = _attribute(0,wave_in_name='FREQHIG')
	freq_num = _attribute(0,wave_in_name='NINTFREQ')
	undu = _attribute(0,wave_in_name='IUNDULATOR')
	wigg = _attribute(0,wave_in_name='IWIGGLER')

	def get_std_paras(self): 
		self.freq_low.set(300) # eV
		self.freq_high.set(500) # eV
		self.freq_num.set(20) # 
		self.undu.set(0)
		self.wigg.set(0)
		return self

class undu_paras(_attribute_collection):
	pkHalbasy = _attribute(0.0,wave_in_name='PKHALBASY')
	b0Halbasy = _attribute(1.0,wave_in_name='B0HALBASY')
	zlHalbasy = _attribute(0.018,wave_in_name='ZLHALBASY')
	ahwpolHalbasy = _attribute(5,wave_in_name='AHWPOL')
	undu_type = _attribute("undu_ellip")

	b0y = _attribute(0,wave_in_name='B0ELLIPV') # magn. field strength amplitude in vertical direction
	b0z = _attribute(0,wave_in_name='B0ELLIPH') # magn. field strength amplitude in horizontal direction
	nper = _attribute(0,wave_in_name='PERELLIP') # number of periods
	perl_x = _attribute(0,wave_in_name='XLELLIP') # period length
	ell_shift = _attribute(0.25,wave_in_name='ELLSHFT') # shift between the two magnetic arrays in fractions of one period

	def get_std_paras(self): 
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
	nthreads = _attribute(2,wave_in_name='MTHREADS')
	zipped = _attribute(True)
	spec_calc = _attribute(False,wave_in_name='ISPEC')
	iemit = _attribute(0,wave_in_name='IEMIT')
	iefold = _attribute(1,wave_in_name='IEFOLD')
	isigusr = _attribute(1,wave_in_name='ISIGUSR')
	ihisascii = _attribute(0,wave_in_name='IHISASCII')

	# undu_type = _attribute('')

	b_type = _attribute('none')
	irbtab = _attribute(0,wave_in_name='IRBTAB')
	irbtabzy = _attribute(0,wave_in_name='IRBTABZY')
	irbtabxyz = _attribute(0,wave_in_name='IRBTABXYZ')
	undu_easy = _attribute(0,wave_in_name='KHALBA') # magnetic structure without specific ends
	undu_endp = _attribute(0,wave_in_name='KHALBASY') # magnetic structure with endpoles
	undu_gap = _attribute(0,wave_in_name='KUNDUGAP') # undulator with analytic gap-variation
	undu_ellip = _attribute(0,wave_in_name='KELLIP') # elliptic undulator
	undu_mode = _attribute('undu_easy') # wave mode

	def get_std_paras(self, undu_mode = 'undu_easy' ): 
		"""
		Getting standard wave-parameters depending on mode
		- undu_mode = 'By' - takes by field data and runs with that
		- undu_mode = 'Byz'
		- undu_mode = 'Bxyz'
		- undu_mode = 'undu_ellip' - standard elliptical undulator
		- undu_mode = 'undu_easy'
		- undu_mode = 'undu_endp'
		- undu_mode = 'undu_gap'
		"""
		
		dir_path = os.path.dirname(os.path.realpath(__file__))
		self.wave_prog_folder.set(dir_path+'/../../External-Software/WAVE/')
		self.in_file_folder.set(dir_path+'/../UNDWAVE_IN_FILES/WAVE-In-Files/')
		self.in_files.set({ 'By' : 'load_ext_on_axis_by_ALL_OUT.in', 
							'Byz' : 'load_ext_on_axis_byz_ALL_OUT.in', 
							'Bxyz' : 'load_ext_on_axis_bxyz_ALL_OUT.in', 
							'undu_ellip' : 'wave.in',
							'none' : 'wave.in'
							})
		self.field_folder.set('')
		self.field_files.set([])
		self.res_folder.set('')
		self.wave_data_res_folder.set('WAVE_DATA/')
		self.pics_folder.set('Pics/')
		self.res_summary_file.set('res_summary.txt')
		self.no_copy.set(['WAVE_CODE.DAT', 'undumag_mu_77K.dat', 'undumag_mu_300K.dat',
								'iron_muinf_sat-2.34.dat', 'Vanadium_Permendur_Radia', 'WAVE.mhb' ])
		self.wave_ending_extract.set([ 'dat', 'wva', 'wvh', 'out' ])
		self.wave_ending_copy.set([ 'in' ])
		self.wave_files_essentials.set([ 'stokes_dist_emittance_espread', 'trajectory',
											'irradiated_power_dist', 'brilliance_3702', 
											'photon_flux_(pinhole)_48000',
											'selected_s0_e_(folded)_x_1_e_6_180000', 'wave.out'])
		self.nthreads.set(2)
		self.zipped.set(False)
		self.wave_res_copy_behaviour.set('copy_all')
		self.zip_res_folder.set(1)
		self.iemit.set(0)
		self.iefold.set(1)
		self.isigusr.set(1)
		self.spec_calc.set(1) # boolean
		self.undu_mode.set(undu_mode)

		self.b_type.set('none')
		self.irbtab.set(0)
		self.irbtabzy.set(0)
		self.irbtabxyz.set(0)
		self.undu_easy.set(0)
		self.undu_endp.set(0)
		self.undu_gap.set(0)
		self.undu_ellip.set(0)

		if undu_mode == 'By' :
			self.b_type.set('By')
			self.irbtab.set(-2)
		elif undu_mode == 'Byz' :
			self.b_type.set('Byz')
			self.irbtabzy.set(1)
		elif undu_mode == 'Bxyz' :
			self.b_type.set('Bxyz')
			self.irbtabxyz.set(1)
		elif undu_mode == 'undu_ellip' :
			self.undu_ellip.set(1)
		elif undu_mode == 'undu_easy' :
			self.undu_easy.set(1)
		elif undu_mode == 'undu_endp' :
			self.undu_endp.set(1)
		elif undu_mode == 'undu_gap' :
			self.undu_gap.set(1)

		return self

