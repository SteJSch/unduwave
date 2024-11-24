
import pdb
import sys
sys.path.insert(0, '../')

try :
	# works when calling script with python3 script_file
	dir_path = os.path.dirname(os.path.realpath(__file__))
except:
	# works when calling script with exec from python console
	dir_path = os.getcwd()

import unduwave as uw

"""
A change 2
"""
wave = uw.wave()
wave_prog_paras = wave._wave_prog_paras
wave_prog_paras.res_folder.set(dir_path+'/res/')

wave.run()

# wave.set_ebam_settings()
# results = wave.run()

# flux_dens = results.get("flux_dens")
# flux_dens.plot_live()
# flux_dens_data = flux_dens.data

# #uw.hello()
# uw.api.hello()

pdb.set_trace()

# # Create wave instance for spectrum calculation from b-field data
# wave = WaveFromB('By') # get a wave instance that calculates from By and Bz

# # Get Wave-parameter
# paras = wave.get_paras()
# # Customize Parameters
# paras.field_folder.set('/home/jerostan/Dokumente/Very_Unimportant/Arbeit/HZB/Undus/Simus/Python_Code/wavepy_simone/examples/Fields/Interpolated/')# Field Folder
# paras.field_files.set( [ 'interp_b_field_gap_18.45_.dat'] )# The magnetic field files to be used in the simulation
# paras.res_folder.set('/home/jerostan/Dokumente/Very_Unimportant/Arbeit/HZB/Undus/Simus/Python_Code/wavepy_simone/examples/Spectrum_Results/') # where to store the wave results
# paras.spec_calc.set(1) # calculate spectrum? if not, only trajectory is calculated
# paras.freq_low.set(100)   # lower frequency for spectrum calc. [eV]
# paras.freq_high.set(105)  # upper frequency for spectrum calc. [eV]
# paras.freq_num.set(5) # number of frequencies for spectrum calc. [eV]
# paras.pinh_w.set(3)   # pinhole width (horizontal-z) [mm]
# paras.pinh_h.set(3)   # pinhole height (vertical-y) [mm]
# paras.pinh_x.set(10) 	  # pinhole distance in beam direction (x) [m]
# paras.pinh_nz.set(21)	  # number of horizontal points in pinhole for spec. calc.
# paras.pinh_ny.set(21)   # number of vetical points in pinhole for spec. calc.
# paras.wave_res_copy_behaviour.set('copy_essentials') #'copy_all'

# # Set Parameters
# wave.set_paras(paras)
# # Run WAVE
# wave.run_wave()