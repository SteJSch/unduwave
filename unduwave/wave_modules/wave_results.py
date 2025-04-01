from unduwave.unduwave_incl import *
from unduwave.wave_modules.wave_parameters import *
import unduwave.helpers.file_folder_helpers as f_h
from unduwave.quantities.quantities import *

class wave_results : 

	def __init__(self,wave_api) :
		"""
		Init wave-result class holding all the results
		"""
		self._wave_api = wave_api
		self._res_quantities = []
		self._res_folder      = self._wave_api._wave_prog_paras.res_folder.get()
		self._res_folder_wave = self._res_folder + self._wave_api._wave_prog_paras.wave_data_res_folder.get()
		self._res_folder_pics = self._res_folder + self._wave_api._wave_prog_paras.pics_folder.get()
		if not os.path.exists(self._res_folder_pics):
			os.makedirs(self._res_folder_pics)

	def load_from_res_folder(self) :
		"""
		goes through folder and loads all results it can find, flux, flux_d, brilliance
		"""
		self.find_load_trajectory_bfield_data()
		self.find_load_power_distribution()
		self.find_load_flux()
		self.find_load_brilliance()
		self.find_load_flux_density_on_axis()

	def find_load_flux_density_distribution(self, energies) : 
		"""
		Finds the distribution that is closed to each energy in energies, loads y/z data and flux distribution and adds them to the
		res_quantities list

		energies - list of energies in eV at which flux distribution is to be plotted
		"""

		file_name = 'stokes_dist_emittance_espread_'
		files = f_h.find_files_exptn(folder = self._res_folder_wave, hints = [file_name], exptns = [])
		en_files = []
		for file in files : 
			data = pd.read_csv( self._res_folder_wave + file, skiprows=3,header=None, dtype=object, delim_whitespace=True)
			data.columns = [ "z", "y", "S0", "S1", "S2", "S3" ]
			cols = data.columns
			for col in cols:
				data[col] = data[col].astype(float)            
			file_lines = []
			with open(self._res_folder_wave+file, 'r') as o_f:
				file_lines = o_f.readlines()
			energy_file = float(file_lines[2].split(':')[-1])
			en_files.append( { 'en' : energy_file, 'file' : file, 'data' : data } )
		en_files_pd = pd.DataFrame(en_files)
		en_files_pd = en_files_pd.sort_values(by=['en'])
		en_files = en_files_pd.to_dict('records')
		en_vals = en_files_pd['en'].to_list()    
		ens_loaded = []
		yz_loaded = False
		for energy in energies :
			dist_bst = None
			ind_bst = 0
			for ind, en_val in enumerate(en_vals) : 
				dist = abs(energy - en_val)
				if dist_bst is None : 
					dist_bst = dist
					ind_bst = ind
				elif dist_bst > dist :
					dist_bst = dist
					ind_bst = ind
			data = en_files[ind_bst]['data']
			name = f'E = {en_vals[ind_bst]:.2f} eV \n'
			save_n = f'_e_{en_vals[ind_bst]:.2f}_'

			ens_loaded.append(round(en_vals[ind_bst],2))
			if not yz_loaded :
				data['z'] = data['z'].apply( lambda x : round(x,2) )
				data['y'] = data['y'].apply( lambda x : round(x,2) )
				z_fluxd_dis = wave_quantity(
					wave_api=self._wave_api,
					data=data['z'].to_list(),
					name="fd_z",
					plot_name="z",
					description='z-Longitudinal Position',
					unit='mm',
					)
				y_fluxd_dis = wave_quantity(
					wave_api=self._wave_api,
					data=data['y'].to_list(),
					plot_name="y",
					name="fd_y",
					description='y-Vertical Position',
					unit='mm',
					)
				self._res_quantities = self._res_quantities + [z_fluxd_dis,y_fluxd_dis] 

				yz_loaded = True
			flux_dens = wave_quantity(
				wave_api=self._wave_api,
				data=data['S0'].to_list(),
				name=f"flux_density_distribution_{ens_loaded[-1]:.2f}",
				plot_name=f'Flux Density Distribution',
				description=f'Flux Density Distribution \n $E={ens_loaded[-1]:.2f}$ eV',
				unit='\n [$\dot{N}_\gamma/(0.1\%\cdot BW\cdot mm^2 )$]\t',
				)			
			self._res_quantities = self._res_quantities + [flux_dens] 

		return ens_loaded

	def find_load_flux_density_on_axis(self) :
		"""
		Loads on-axis flux-density
		"""

		for filename in os.listdir(self._res_folder_wave):
			if (filename.find('selected_s0_e_(folded)_x_1_e_6_180000') >= 0): 
				try:
					data = pd.read_csv( self._res_folder_wave + filename, skiprows=3,header=None, dtype=object, sep='\s+')
				except:
					return
				data.columns = [ 'en', 'flux_dens' ]
				cols = data.columns
				for col in cols:
					data[col] = data[col].astype(float)
				data = data[ data['flux_dens'] > 1e9 ]
				data['en'] = data['en'].apply( lambda x : x * 1e-3 )

				en_fd = wave_quantity(
					wave_api=self._wave_api,
					data=data['en'].to_list(),
					name="en_fd",
					plot_name="E",
					description='Photon Energy',
					unit='keV',
					)
				fd = wave_quantity(
					wave_api=self._wave_api,
					data=data['flux_dens'].to_list(),
					plot_name="Flux Density",
					name="flux_density",
					description='Flux Density - W. E-Spread,Emittance',
					unit='$\\dot{N}_\\gamma/(0.1\\%\\cdot BW \\cdot mm^2 )$',
					)
				self._res_quantities = self._res_quantities + [en_fd,fd] 

	def find_load_brilliance(self) :
		"""
		Loads brilliance data
		"""
		for filename in os.listdir(self._res_folder_wave):
			if (filename.find('brilliance_3702') >= 0): 
				try:
					data = pd.read_csv( self._res_folder_wave + filename, skiprows=6,header=None, dtype=object, delim_whitespace=True)
				except:
					return
				data.columns = [ "en", "b0", "b1", "b2", "b3", "b0f", "b1f", "b2f", "b3f", "b0e",\
							"b1e", "b2e", "b3e", "b0ef", "b1ef", "b2ef", "b3ef" ]
				cols = data.columns
				for col in cols:
					data[col] = data[col].astype(float)
				data['en'] = data['en'].apply( lambda x : x * 1e-3 )
				data['b0'] = data['b0'].apply( lambda x : x*1e-12 )
				data['b0e'] = data['b0e'].apply( lambda x : x*1e-12 )
				data['b0f'] = data['b0f'].apply( lambda x : x*1e-12 )
				data['b0ef'] = data['b0ef'].apply( lambda x : x*1e-12 )

				en_brill = wave_quantity(
					wave_api=self._wave_api,
					data=data['en'].to_list(),
					name="en_brill",
					plot_name="E",
					description='Photon Energy',
					unit='keV',
					)
				brill0 = wave_quantity(
					wave_api=self._wave_api,
					data=data['b0'].to_list(),
					plot_name="B0",
					name="brill0",
					description='Brilliance - No E-Spread,No Emittance',
					unit='$\dot{N}_\gamma/(0.1\% \cdot BW \cdot mm^2 \cdot mrad^2)$',
					)
				brill0e = wave_quantity(
					wave_api=self._wave_api,
					data=data['b0e'].to_list(),
					plot_name="B0e",
					name="brill0e",
					description='Brilliance - W. E-Spread,No Emittance',
					unit='$\dot{N}_\gamma/(0.1\% \cdot BW \cdot mm^2 \cdot mrad^2)$',
					)
				brill0ef = wave_quantity(
					wave_api=self._wave_api,
					data=data['b0ef'].to_list(),
					plot_name="B0ef",
					name="brill0ef",
					description='Brilliance - W. E-Spread and Emittance',
					unit='$\dot{N}_\gamma/(0.1\% \cdot BW \cdot mm^2 \cdot mrad^2)$',
					)
				brill0f = wave_quantity(
					wave_api=self._wave_api,
					data=data['b0f'].to_list(),
					plot_name="B0f",
					name="brill0f",
					description='Brilliance - No E-Spread, W. Emittance',
					unit='$\dot{N}_\gamma/(0.1\% \cdot BW \cdot mm^2 \cdot mrad^2)$',
					)
				self._res_quantities = self._res_quantities + [en_brill,brill0,brill0e,brill0ef,brill0f] 

	def find_load_flux(self) :
		"""
		Loads flux through pinhole
		"""
		for filename in os.listdir(self._res_folder_wave):
			if (filename.find('stokes_pinhole_emittance_espread') >= 0): 
				try:
					data = pd.read_csv( self._res_folder_wave + filename, skiprows=3,header=None, dtype=object, delim_whitespace=True)
				except:
					return
				data.columns = [ 'en', 'flux', 's1', 's2', 's3' ]
				cols = data.columns
				for col in cols:
					data[col] = data[col].astype(float)
				data['en'] = data['en'].apply( lambda x : x * 1e-3 )

				en_flux = wave_quantity(
					wave_api=self._wave_api,
					data=data['en'].to_list(),
					name="en_flux",
					plot_name="E",
					description='Photon Energy',
					unit='keV',
					)
				flux = wave_quantity(
					wave_api=self._wave_api,
					data=data['flux'].to_list(),
					plot_name="Flux",
					name="flux",
					description='Total Flux through Pinhole',
					unit='Nphot/(s 0.1% BW)',
					)
				self._res_quantities = self._res_quantities + [en_flux,flux] 

	def find_load_power_distribution(self) :
		"""
		Loads power-distribution data
		"""
		for filename in os.listdir(self._res_folder_wave):
			if (filename.find('irradiated_power_dist') >= 0): 
				try:
					data = pd.read_csv( self._res_folder_wave + filename, skiprows=2,header=None, dtype=object, delim_whitespace=True)
				except:
					return
				data.columns = [ "z", "y", "power" ]
				cols = data.columns
				for col in cols:
					data[col] = data[col].astype(float)

				data['z'] = data['z'].apply( lambda x : round(x,2) )
				data['y'] = data['y'].apply( lambda x : round(x,2) )
				z_power_d = wave_quantity(
					wave_api=self._wave_api,
					data=data['z'].to_list(),
					name="power_z",
					plot_name="z",
					description='z-Longitudinal Position',
					unit='mm',
					)
				y_power_d = wave_quantity(
					wave_api=self._wave_api,
					data=data['y'].to_list(),
					plot_name="y",
					name="power_y",
					description='y-Vertical Position',
					unit='mm',
					)
				power_d = wave_quantity(
					wave_api=self._wave_api,
					data=data['power'].to_list(),
					name="power_distribution",
					plot_name='Power Distribution',
					description='Total Power Distribution',
					unit='$W/m^2$',
					)
				self._res_quantities = self._res_quantities + [z_power_d,y_power_d,power_d] 

	def find_load_trajectory_bfield_data(self) :
		"""
		Loads trajectory-bfield data
		"""
		for filename in os.listdir(self._res_folder_wave):

			if (filename.find('trajectory.wva') >= 0) : 
				try:
					data = pd.read_csv( self._res_folder_wave + filename, skiprows=3,header=None, dtype=object, delim_whitespace=True)
				except:
					return
				data.columns = [ "x", "y", "z", "Bx", "By", "Bz" ]
				cols = data.columns
				for col in cols:
					data[col] = data[col].astype(float)
				# data = data.to_dict('records')

				x_quant = wave_quantity(
					wave_api=self._wave_api,
					data=data['x'].to_list(),
					plot_name="x",
					name='traj_x',
					description='x-Longitudinal Position',
					unit='m',
					)
				y_quant = wave_quantity(
					wave_api=self._wave_api,
					data=data['y'].to_list(),
					plot_name="y",
					name='traj_y',
					description='$y$-Vertical Position',
					unit='m',
					)
				z_quant = wave_quantity(
					wave_api=self._wave_api,
					data=data['z'].to_list(),
					plot_name="$z$",
					name='traj_z',
					description='$z$-Horizontal Position',
					unit='m',
					)
				Bx_quant = wave_quantity(
					wave_api=self._wave_api,
					data=data['Bx'].to_list(),
					name="Bx",
					plot_name='$B_x$',
					description='$B_x$-Longitudinal Magnetic Induction',
					unit='T',
					)
				By_quant = wave_quantity(
					wave_api=self._wave_api,
					data=data['By'].to_list(),
					name="By",
					plot_name='$B_y$',
					description='$B_y$-Vertical Magnetic Induction',
					unit='T',
					)
				Bz_quant = wave_quantity(
					wave_api=self._wave_api,
					data=data['Bz'].to_list(),
					name="Bz",
					plot_name='$B_z$',
					description='$B_z$-Horizontal Magnetic Induction',
					unit='T',
					)
				self._res_quantities = self._res_quantities + [x_quant,y_quant,z_quant,Bx_quant,By_quant,Bz_quant] 
				break

	def get_result(self,which) :
		"""
		Pass string "flux_dens",... to get object containing data and some functionality
		"""
		for quant in self._res_quantities :
			if quant._name == which :
				return quant
