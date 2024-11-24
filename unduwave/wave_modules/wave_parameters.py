from unduwave.unduwave_incl import *
from unduwave.attribute_classes.attributes import _attribute
from unduwave.attribute_classes.attributes import _attribute_collection

class wave_parameters(_attribute_collection):
    """
    Represents standard parameters for wave simulations.

    Args:
        b_type (str, optional): Type of the wave element.
    """
    b_type = _attribute('By')
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
    # WAVE-IN PARAMETERS
    zipped = _attribute(True)
    freq_low = _attribute(0,wave_in_name='FREQLOW')
    freq_high = _attribute(0,wave_in_name='FREQHIG')
    freq_num = _attribute(0,wave_in_name='NINTFREQ')
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

    pinh_w = _attribute(0,wave_in_name='PINW',fac=1e-3)
    pinh_h = _attribute(0,wave_in_name='PINH',fac=1e-3)
    spec_calc = _attribute(False,wave_in_name='ISPEC')
    pinh_x = _attribute(0,wave_in_name='PINCEN(1)')
    pinh_nz = _attribute(0,wave_in_name='MPINZ')
    pinh_ny = _attribute(0,wave_in_name='MPINY')

    # Simple Undulator Parameters
    undu_type = _attribute('')

    circumference = _attribute(240,wave_in_name='UMFANG')
    rdipol = _attribute(4.359,wave_in_name='RDIPOL')
    iemit = _attribute(0,wave_in_name='IEMIT')
    irbtab = _attribute(0,wave_in_name='IRBTAB')
    iefold = _attribute(1,wave_in_name='IEFOLD')
    isigusr = _attribute(1,wave_in_name='ISIGUSR')
    ihisascii = _attribute(0,wave_in_name='IHISASCII')
    undu = _attribute(0,wave_in_name='IUNDULATOR')
    wigg = _attribute(0,wave_in_name='IWIGGLER')
    undu_easy = _attribute(0,wave_in_name='KHALBA') # magnetic structure without specific ends
    undu_endp = _attribute(0,wave_in_name='KHALBASY') # magnetic structure with endpoles
    undu_gap = _attribute(0,wave_in_name='KUNDUGAP') # undulator with analytic gap-variation
    undu_ellip = _attribute(0,wave_in_name='KELLIP') # elliptic undulator

    pkHalbasy = _attribute(0.0,wave_in_name='PKHALBASY')
    b0Halbasy = _attribute(1.0,wave_in_name='B0HALBASY')
    zlHalbasy = _attribute(0.018,wave_in_name='ZLHALBASY')
    ahwpolHalbasy = _attribute(5,wave_in_name='AHWPOL')

    b0y = _attribute(0,wave_in_name='B0ELLIPV') # magn. field strength amplitude in vertical direction
    b0z = _attribute(0,wave_in_name='B0ELLIPH') # magn. field strength amplitude in horizontal direction
    nper = _attribute(0,wave_in_name='PERELLIP') # number of periods
    perl_x = _attribute(0,wave_in_name='XLELLIP') # period length
    ell_shift = _attribute(0.25,wave_in_name='ELLSHFT') # shift between the two magnetic arrays in fractions of one period

    def get_std_paras(self, b_type = 'By'): 
        """
		Returns a object containing WAVE standard parameters updated with standard values

		Args:
			b_type (str): Type of b-calculation. Options: 'By' (only By given), 'Byz', 'Bxyz' (each b-field in different file).
			wave_prog_folder (str): Main folder where wave is stored (ending on '/').
			in_file_folder (str): Folder in which the wave in-files are stored (ending on '/').
			in_files (dict): Dictionary of wave in-files for different b_type situations. Format: {b_type: wave_in_file}.
			field_folder (str): Folder where b-field files are stored.
			field_files (list): List of b-field files. Format: 2 cols - x[mm] and B[T], no separator, no headers.
			res_folder (str): Folder where results are stored (ending on '/').
			wave_data_res_folder (str): Subfolder of res_folder where wave data is stored.
			pics_folder (str): Subfolder of res_folder where pictures are stored.
			res_summary_file (str): Name of the summary file to be written.
			no_copy (list): List of file names not to be copied/moved after simulation from wave stage folder.
			wave_ending_extract (list): List of file endings to move from wave stage folder after simulation.
			wave_ending_copy (list): List of file endings of files to be copied (not moved).
			wave_res_copy_behaviour (str): Behavior for copying wave results: 'copy_all', 'copy_del_none', 'copy_essentials'.
			wave_files_essentials (list): List of essential files when wave_res_copy_behaviour is set to 'copy_essentials'.
			zip_res_folder (bool): Truth value, whether to zip results or not.
			freq_low (float): Lower frequency (energy) of spectrum to calculate [eV].
			freq_high (float): Upper frequency (energy) of spectrum to calculate [eV].
			freq_num (int): Number of frequencies to calculate.
			beam_en (float): Beam energy in [GeV].
			current (float): Current in [A].
			pinh_w (float): Pinhole width (horizontal-z) [mm].
			pinh_h (float): Pinhole height (vertical-y) [mm].
			spec_calc (bool): Truth value, whether to calculate spectrum or only write trajectory.
			pinh_x (float): Position of pinhole along optical axis [m].
			pinh_nz (int): Number of points in pinhole horizontally.
			pinh_ny (int): Number of points in pinhole vertically.
		"""
        
        dir_path = os.path.dirname(os.path.realpath(__file__))
        
        self.b_type.set(b_type)
        self.wave_prog_folder.set(dir_path+'/WAVE/')
        self.in_file_folder.set(dir_path+'/WAVE-In-Files/')
        self.in_files.set({ 'By' : 'load_ext_on_axis_by_ALL_OUT.in', 
			 					'Byz' : 'load_ext_on_axis_byz_ALL_OUT.in', 
								'Bxyz' : 'load_ext_on_axis_bxyz_ALL_OUT.in', 'undu_ellip' : 'wave.in' })
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
        self.zipped.set(True)
        self.wave_res_copy_behaviour.set('copy_all')
        self.zip_res_folder.set(1)
        self.freq_low.set(300) # eV
        self.freq_high.set(500) # eV
        self.freq_num.set(5) # 
        self.beam_en.set(1.722) # [GeV]
        self.current.set(0.3) # [A]
        self.irbtab.set(-2)

        self.bsigz.set(275e-6) # 
        self.bsigzp.set(28.1e-6) #
        self.bsigy.set(22.5e-6) # 
        self.bsigyp.set(6.8e-6) # 
        self.espread.set(1e-3) # 
        self.emitt_h.set(7.7e-9)
        self.emitt_v.set(15.4e-11)
        self.betfunh.set(0)
        self.betfunv.set(0)

        self.pinh_w.set(3) # [mm]
        self.pinh_h.set(3) # [mm]
        self.spec_calc.set(1) # boolean
        self.pinh_x.set(10) # [m]
        self.pinh_nz.set(21) # 
        self.pinh_ny.set(21) # 

        self.iemit.set(0)
        self.circumference.set(240)
        self.rdipol.set(4.359)
        self.iefold.set(1)
        self.isigusr.set(1)
        self.undu.set(0)
        self.wigg.set(0)
        self.undu_easy.set(0)
        self.undu_endp.set(0)
        self.undu_gap.set(0)

        if self.b_type.get() == 'undu_ellip' :
            self.undu_ellip.set(1)
            self.b0y.set(0.3)
            self.b0z.set(1.0)
            self.nper.set(5)
            self.perl_x.set(0.02)
            self.ell_shift.set(0.25)

        return self

