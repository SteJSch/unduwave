from unduwave.unduwave_incl import *
from unduwave.wave_modules.wave_parameters import *
import unduwave.helpers.file_folder_helpers as f_h
from unduwave.quantities.quantities import *

class wave_results : 

	def __init__(self,wave_api) :
		"""
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

	def find_load_power_distribution(self) :
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
					data=data['z'].to_list(),
					name="power_z",
					plot_name="z",
					description='z-Longitudinal Position',
					unit='m',
					)
				y_power_d = wave_quantity(
					data=data['z'].to_list(),
					plot_name="y",
					name="power_y",
					description='y-Vertical Position',
					unit='m',
					)
				power_d = wave_quantity(
					data=data['power'].to_list(),
					name="power_distribution",
					plot_name='Power Distribution',
					description='Total Power Distribution',
					unit='$W/m^2$',
					)
				self._res_quantities = self._res_quantities + [z_power_d,y_power_d,power_d] 

	def find_load_trajectory_bfield_data(self) :
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
					data=data['x'].to_list(),
					plot_name="x",
					name='traj_x',
					description='x-Longitudinal Position',
					unit='m',
					)
				y_quant = wave_quantity(
					data=data['y'].to_list(),
					plot_name="y",
					name='traj_y',
					description='$y$-Vertical Position',
					unit='m',
					)
				z_quant = wave_quantity(
					data=data['z'].to_list(),
					plot_name="$z$",
					name='traj_z',
					description='$z$-Horizontal Position',
					unit='m',
					)
				Bx_quant = wave_quantity(
					data=data['Bx'].to_list(),
					name="Bx",
					plot_name='$B_x$',
					description='$B_x$-Longitudinal Magnetic Induction',
					unit='T',
					)
				By_quant = wave_quantity(
					data=data['By'].to_list(),
					name="By",
					plot_name='$B_y$',
					description='$B_y$-Vertical Magnetic Induction',
					unit='T',
					)
				Bz_quant = wave_quantity(
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

class ExtractAndPlot():

	def __init__(self, wave_paras):
		self.wave_paras = wave_paras

	
	def extract_wave_results(self, results: list, plot: bool = False, en_range: list = [], xlim: list = [], ylim: list = [], save_data_w_pics: bool = False, list_ex: bool = True, vert_lines : list = [], add_str = ''):
		"""
		Extracts the desired results and returns a DataFrame with them.

		Args:
			results (list): A list containing the desired result types. Valid options are:
				- 'traj_magn': Extracts Trajectory and magnetic field data.
				- 'power_distr': Power distribution [W/mm^2].
				- 'flux': Flux.
				- 'flux_dens': Flux Density.
				- 'flux_dens_distr': Flux Density Distribution.
				- 'power_dens_distr': Power density distribution.
				- 'brill': Brilliance.
			plot (bool): Determines if loaded data is plotted or not.
			en_range (list): For distribution, gives the energy range over which the data is integrated.
			xlim (list): Sets the limits of the distribution plots in the horizontal direction.
			ylim (list): Sets the limits of the distribution plots in the vertical direction. For 2D plots, only xlim is considered if given.
			save_data_w_pics (bool): If True, saves the data for each plot alongside the picture.
			
		Returns:
			list: A list of the loaded DataFrames.
		"""
		data_ret        = []
		res_folder      = self.wave_paras.res_folder.get()
		res_folder_wave = res_folder + self.wave_paras.wave_data_res_folder.get()
		res_folder_pics = res_folder + self.wave_paras.pics_folder.get()

		fontsize_max = 18
		plt.rcParams.update({'font.size': fontsize_max-2})

		if not os.path.exists(res_folder_pics):
			os.makedirs(res_folder_pics)

		if self.wave_paras.zipped.get() :
			zip_files = f_h.unzip_zip_in_folder(folder = res_folder_wave)
		for result in results :

			if (result == 'flux_dens_distr') : 
				ret = self.load_plot_stokes_distrib(folder = res_folder_wave, en_range = en_range, xlim = xlim, ylim = ylim, plot = plot, save_folder = res_folder_pics, 
								save_data_w_pics = save_data_w_pics, power_distr = False, list_ex = list_ex )
				data_ret.append( { 'result' : 'flux_dens_distr', 'data' : ret } )
				continue
			if (result == 'power_dens_distr') : 
				ret = self.load_plot_stokes_distrib(folder = res_folder_wave, en_range = en_range, xlim = xlim, ylim = ylim, plot = plot, save_folder = res_folder_pics, 
								save_data_w_pics = save_data_w_pics, power_distr = True )
				data_ret.append( { 'result' : 'power_dens_distr', 'data' : ret } )
				continue

			for filename in os.listdir(res_folder_wave):

				if (filename.find('trajectory.wva') >= 0) and (result == 'traj_magn') : 
					data = pd.read_csv( res_folder_wave + filename, skiprows=3,header=None, dtype=object, delim_whitespace=True)
					data.columns = [ "x", "y", "z", "Bx", "By", "Bz" ]
					data_mag = data[ ['x','By'] ]
					cols = data.columns
					for col in cols:
						data[col] = data[col].astype(float)
					data_ret.append( { 'result' : 'traj_magn', 'data' : data } )
					fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
					plt.plot(data['x'],data['Bx'], '-', label = '$B_x$')
					plt.plot(data['x'],data['By'], '-', label = '$B_y$')
					plt.plot(data['x'],data['Bz'], '-', label = '$B_z$')
					fig.suptitle("Magnetic Field", fontsize=14)
					plt.xlabel('x [m]', fontsize=12),
					plt.ylabel('Magnetic Field [T]', fontsize=12)
					ax = plt.gca()
					ax.legend(loc='best', bbox_to_anchor=(0.8, 0.5, 0.0, 0.0))	
					plt.savefig(res_folder_pics+f"magn-field{add_str}.png"   , bbox_inches='tight')

					fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
					data['z'] = data['z'].apply( lambda x : x * 1e6 )
					data['y'] = data['y'].apply( lambda x : x * 1e6 )
					plt.plot(data['x'],data['z'], '-', label = '$z(x)$')
					plt.plot(data['x'],data['y'], '-', label = '$y(x)$')
					fig.suptitle("Trajectory", fontsize=14)
					plt.xlabel('x [m]', fontsize=12),
					plt.ylabel('Deflection [$\mu$m]', fontsize=12)
					ax = plt.gca()
					ax.legend(loc='best', bbox_to_anchor=(0.8, 0.5, 0.0, 0.0))	
					plt.savefig(res_folder_pics+f"trajectory{add_str}.png"   , bbox_inches='tight')

					fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
					fig.suptitle("Trajectory-3D", fontsize=14)
					ax = fig.add_subplot(111, projection='3d')
					ax.plot3D(data['x'], data['y'], data['z'], '-')
					ax.set_xlabel('x [m]', fontsize=12)
					ax.set_ylabel('y\n[$\mu m$]', fontsize=12)
					ax.set_zlabel('z [$\mu m$]\t',linespacing=4, fontsize=12)
					plt.savefig(res_folder_pics+f"trajectory-3D{add_str}.png"   , bbox_inches='tight')

					if save_data_w_pics :
						data.columns = [ "x [m]", "y [mu m]", "z [mu m]", "Bx [T]", "By [T]", "Bz [T]" ]
						data.to_csv(res_folder_pics+f"traj_mag_{add_str}.csv",index=False, sep = ' ')
						d_list = data.to_dict('records')
						d_list = d_list[::2]
						new_list = []
						for el in d_list : 
							new_list.append( { 'x' : el['x [m]']*1e3, 'By':el['By [T]'] } )
							# el.update({'row1' : 3, 'row2': 4})
							# new_list.append( el )
						pd.DataFrame(new_list).to_csv(res_folder_pics+f"mag_{add_str}.csv",index=False, sep = ' ', header = False)
					break

				# if (filename.find('power_dist_emittance') >= 0) and (plot_hint == 'power_distr') : 
				if (filename.find('irradiated_power_dist') >= 0) and (result == 'power_distr') : 
					data = pd.read_csv( res_folder_wave + filename, skiprows=2,header=None, dtype=object, delim_whitespace=True)
					data.columns = [ "z", "y", "power" ]
					cols = data.columns
					for col in cols:
						data[col] = data[col].astype(float)
					data_ret.append( { 'result' : 'power_distr', 'data' : data } )
					fig = plt.figure(figsize=(13*cm_inch, 8.5*cm_inch), dpi=150)

					fig.suptitle("Power Density Distribution", fontsize=14)
					ax = fig.add_subplot(111, projection='3d')
					ax.plot_trisurf(data['y'], data['z'], data['power'], cmap=cm.jet, linewidth=0.2)

					if len(xlim) > 0 :
						ax.set_xlim( xlim )
					if len(ylim) > 0 :
						ax.set_ylim( ylim )
					ax.set_xlabel('y [mm]', fontsize=12)
					ax.set_ylabel('z [mm]', fontsize=12)
					# ax.set_zlabel('Power [$W/mm^2$]', fontsize=12)
					plt.tight_layout()
					plt.savefig(res_folder_pics+"power_distro.png"   , bbox_inches='tight')
					if save_data_w_pics :
						data_save=copy.deepcopy(data)
						data_save.columns = [ 'z [mm]', 'y [mm]', 'power [W/mm^2]' ]
						data_save.to_csv(res_folder_pics+"power.csv",index=False, sep = ' ')

					fig = plt.figure(figsize=(13*cm_inch, 8.5*cm_inch), dpi=150)
					fig.suptitle("Power Density Distribution", fontsize=14)
					ax = plt.gca()

					data = data[ (data['y'] >= -3) & (data['y'] <= 3) ]
					data['z'] = data['z'].apply( lambda x : round(x,2) )
					data['y'] = data['y'].apply( lambda x : round(x,2) )

					y_data=np.unique(data['y'].to_list())
					z_data=np.unique(data['z'].to_list())
					fun_data=np.array(data['power'].to_list())
					Y_data,Z_data = np.meshgrid(z_data,y_data)

					Funs=fun_data.reshape(len(y_data),len(z_data))

					plt.tight_layout()
					cmap = plt.colormaps["plasma"]
					cmap = cmap.with_extremes(bad=cmap(0))
					pcm = ax.pcolormesh(Y_data,Z_data,Funs, cmap=cmap)
					fig.colorbar(pcm, ax=ax, label='Power [$W/mm^2$]', pad=0)                    
					ax.set_xlabel('z [mm]', fontsize=12)
					ax.set_ylabel('y [mm]', fontsize=12)
					plt.savefig(res_folder_pics+"power_distro_heat.png"   , bbox_inches='tight')

					break

				if (filename.find('brilliance_3702') >= 0) and (result == 'brill') : 
					data = pd.read_csv( res_folder_wave + filename, skiprows=6,header=None, dtype=object, delim_whitespace=True)
					data.columns = [ "en", "b0", "b1", "b2", "b3", "b0f", "b1f", "b2f", "b3f", "b0e",\
								"b1e", "b2e", "b3e", "b0ef", "b1ef", "b2ef", "b3ef" ]
					cols = data.columns
					for col in cols:
						data[col] = data[col].astype(float)
					data['en'] = data['en'].apply( lambda x : x * 1e-3 )
					if len(xlim) > 0 :
						data = data[ (data['en'] >= xlim[0]) & (data['en'] <= xlim[1]) ]
					data_ret.append( { 'result' : 'brill', 'data' : data } )
					fig = plt.figure(figsize=(13*cm_inch, 10.5*cm_inch), dpi=150)
					data['b0'] = data['b0'].apply( lambda x : x*1e-12 )
					data['b0e'] = data['b0e'].apply( lambda x : x*1e-12 )
					data['b0f'] = data['b0f'].apply( lambda x : x*1e-12 )
					data['b0ef'] = data['b0ef'].apply( lambda x : x*1e-12 )
					# plt.plot(data['en'],data['b0'], 'x-', label = 'Flux')
					# plt.plot(data['en'],data['b0e'], 'x-', label = 'Flux')
					# plt.plot(data['en'],data['b0f'], 'x-', label = 'Flux')
					en_data = data['en'].to_list()
					plt_data=[]
					plt_en = []
					for ind_d, val in enumerate(data['b0ef'].to_list()) : 
						if val > 1 :
							plt_data.append(val)
							plt_en.append(en_data[ind_d])
					# plt.plot(data['en'],data['b0ef'], '-', label = 'Flux')
					plt.loglog(plt_en,plt_data, '-', label = 'Flux')
					if len(vert_lines) > 0 :
						for line in vert_lines :
							plt.axvline(x=line,color='red')
					fig.suptitle("Brilliance", fontsize=fontsize_max)
					if len(en_range) > 0 :
						lim_l = en_range[0] * 1e-3
						lim_u = en_range[1] * 1e-3
						ax = plt.gca()
						ax.set_xlim( [lim_l,lim_u] )
					plt.xlabel('Energy [keV]'),
					plt.ylabel('[$\dot{N}_\gamma/(0.1\% \cdot BW \cdot mm^2 \cdot mrad^2)$]')
					ax = plt.gca()
					if len(ylim) > 0 :
						ax.set_ylim( ylim )
					plt.tight_layout()
					plt.savefig(res_folder_pics+"brilliance_ef.png"   , bbox_inches='tight')
					data = pd.DataFrame( { 'en' : data['en'], 'brill' : data['b0ef'] } )
					if save_data_w_pics :
						data.columns = [ 'en [keV]', 'brill [Nphot/(s 0.1% BW mm^2 mrad^2) ]' ]
						data.to_csv(res_folder_pics+"brill.csv",index=False, sep = ' ')
					break

				# Flux with emittance and energy spread through pinnhole
				# if (filename.find('s0_e_(pinhole__folded)_80000') >= 0) and (result == 'flux') : 
				if (filename.find('stokes_pinhole_emittance_espread') >= 0) and (result == 'flux') : 
				# if (filename.find('stokes_pinhole.wva') >= 0) and (result == 'flux') : 
					data = pd.read_csv( res_folder_wave + filename, skiprows=3,header=None, dtype=object, delim_whitespace=True)
					data.columns = [ 'en', 'flux', 's1', 's2', 's3' ]
					cols = data.columns
					for col in cols:
						data[col] = data[col].astype(float)
					data['en'] = data['en'].apply( lambda x : x * 1e-3 )
					if len(xlim) > 0 :
						data = data[ (data['en'] >= xlim[0]) & (data['en'] <= xlim[1]) ]
					data_ret.append( { 'result' : 'flux', 'data' : data } )
					fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
					plt.plot(data['en'],data['flux'], '-', label = 'Flux')
					plt.loglog(data['en'],data['flux'], '-', label = 'Flux')
					if len(vert_lines) > 0 :
						for line in vert_lines :
							plt.axvline(x=line,color='red')
					fig.suptitle("Flux Through Pinhole", fontsize=14)
					if len(en_range) > 0 :
						lim_l = en_range[0] * 1e-3
						lim_u = en_range[1] * 1e-3
						ax = plt.gca()
						ax.set_xlim( [lim_l,lim_u] )
					plt.xlabel('Energy [keV]', fontsize=12),
					plt.ylabel('Flux [$\dot{N}_\gamma/(0.1\%\cdot BW )$]', fontsize=12)
					ax = plt.gca()
					if len(ylim) > 0 :
						ax.set_ylim( ylim )
					plt.tight_layout()
					plt.savefig(res_folder_pics+"flux_pinhole_ef.png"   , bbox_inches='tight')
					if save_data_w_pics :
						data.columns = [ 'en [keV]', 'flux [Nphot/(s 0.1% BW) ]', 's1', 's2', 's3' ]
						data.to_csv(res_folder_pics+"flux.csv",index=False, sep = ' ')
					break

				# Flux Density with emittance and energy spread through pinnhole
				if (filename.find('selected_s0_e_(folded)_x_1_e_6_180000') >= 0) and (result == 'flux_dens') : 
					data = pd.read_csv( res_folder_wave + filename, skiprows=3,header=None, dtype=object, delim_whitespace=True)
					data.columns = [ 'en', 'flux_dens' ]
					cols = data.columns
					for col in cols:
						data[col] = data[col].astype(float)
					data = data[ data['flux_dens'] > 1e9 ]
					data['en'] = data['en'].apply( lambda x : x * 1e-3 )
					data_y = data['flux_dens'].to_list()
					data_ret.append( { 'result' : 'flux_dens', 'data' : data } )
					fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
					if len(xlim) > 0 :
						data = data[ (data['en'] >= xlim[0]) & (data['en'] <= xlim[1]) ]
					plt.loglog(data['en'],data['flux_dens'], '-', label = 'Flux Density')
					# plt.plot(data['en'],data['flux_dens'], '-', label = 'Flux Density')
					# plt.xscale('log')
					# plt.yscale('log')
					if len(vert_lines) > 0 :
						for line in vert_lines :
							plt.axvline(x=line,color='red')
					ax = plt.gca()
					if len(en_range) > 0 :
						lim_l = en_range[0] * 1e-3
						lim_u = en_range[1] * 1e-3
						ax = plt.gca()
						ax.set_xlim( [lim_l,lim_u] )
					# fig.suptitle("Spectral Flux", fontsize=14)
					fig.suptitle("Flux Density", fontsize=14)
					plt.xlabel('Energy [keV]', fontsize=12),
					plt.ylabel('[$\\dot{N}_\\gamma/(0.1\\%\\cdot BW \\cdot mm^2 )$]', fontsize=12)
					# ax.set_ylim( [1e9,1e16] )
					if len(ylim) > 0 :
						ax.set_ylim( ylim )
					plt.savefig(res_folder_pics+"flux_dens_pinhole_ef.png"   , bbox_inches='tight')
					if save_data_w_pics :
						data.columns = [ 'en [keV]', 'flux_dens [Nphot/(s 0.1% BW mm^2) ]' ]
						data.to_csv(res_folder_pics+"flux_dens.csv",index=False, sep = ' ')
					break

		if plot : 
			# plt.pause(0.25)
			# plot.ion()
			plt.show()

		if self.wave_paras.zipped.get() :
			if len(zip_files) > 0 :
				f_h.del_all_files(exptns = [zip_files[-1]], folder = res_folder_wave)
		return data_ret

	def load_plot_stokes_distrib(self,folder, en_range = [], xlim = [], ylim = [], plot = True,\
			 save_folder = '',save_data_w_pics = False, power_distr = False, list_ex = False ) : 
		"""
		loads stokes (flux-density) data from a folder, integrates the flux-density distribution over the en_range 
		and cuts the plot horizontally and vertically to xlim and ylim, plots if plot = True and saves pic to save_folder
		returns the loaded and integrated data object
		If power_distr is True, then not the flux density but power density in the given energy
		range is calculated (multiplying all fluxes by the appropriate energies)
		"""
		fontsize_max = 18
		plt.rcParams.update({'font.size': fontsize_max-2})
		expn = []
		# Find all files that fit the description
		file_name = 'stokes_dist_emittance_espread_'
		files = f_h.find_files_exptn(folder = folder, hints = [file_name], exptns = [])
		en_files = []
		for file in files : 
			data = pd.read_csv( folder + file, skiprows=3,header=None, dtype=object, delim_whitespace=True)
			data.columns = [ "z", "y", "S0", "S1", "S2", "S3" ]
			cols = data.columns
			for col in cols:
				data[col] = data[col].astype(float)            
			file_lines = []
			with open(folder+file, 'r') as o_f:
				file_lines = o_f.readlines()
			energy = float(file_lines[2].split(':')[-1])
			en_files.append( { 'en' : energy, 'file' : file, 'data' : data } )
		en_files_pd = pd.DataFrame(en_files)
		# sort all files (flux data over pinhole for specific energy)
		en_files_pd = en_files_pd.sort_values(by=['en'])
		en_files = en_files_pd.to_dict('records')
		en_vals = en_files_pd['en'].to_list()    
		# en_spread = (en_vals[-1] - en_vals[0]) / ( len(en_vals) - 1 )
		# no energy range procided: then take whole measurement range
		if len(en_range) < 1 :
			en_range = [en_vals[0], en_vals[-1]]
		else :
			if en_range[0] < en_vals[0] : 
				en_range[0] = en_vals[0]
			if en_range[-1] > en_vals[-1] : 
				en_range[-1] = en_vals[-1]
		if list_ex :
			for en_val_plt in en_range :

				dist_bst = None
				ind_bst = 0
				for ind, en_val in enumerate(en_vals) : 
					dist = abs(en_val_plt - en_val)
					if dist_bst is None : 
						dist_bst = dist
						ind_bst = ind
					elif dist_bst > dist :
						dist_bst = dist
						ind_bst = ind
				spec = en_files[ind_bst]['data']
				name = f'E = {en_vals[ind_bst]:.2f} keV \n'
				save_n = f'_e_{en_vals[ind_bst]:.2f}_'

				fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
				# fig.suptitle("Flux Distribution " + name, fontsize=14)
				# fig.suptitle("Flux Distribution", fontsize=14)
				ax = fig.add_subplot(111, projection='3d')
				ax.plot_trisurf(spec['z'], spec['y'], spec['S0'], cmap=cm.jet, linewidth=0.2)
				ax.set_xlabel('z [mm]', fontsize=12)
				ax.set_ylabel('y [mm]', fontsize=12)
				ax.set_zlabel('\n [$\dot{N}_\gamma/(0.1\%\cdot BW\cdot mm^2 )$]\t',linespacing=4, fontsize=12)
				if len(xlim) > 0 :
					ax.set_xlim( xlim )
				if len(ylim) > 0 :
					ax.set_ylim( ylim )
				plt.tight_layout()
				plt.savefig(save_folder+"flux_density_distrib_ef" + save_n + ".png" , bbox_inches='tight')
				if save_data_w_pics :
					spec_save = copy.deepcopy(spec)
					spec_save.columns = [ 'z [mm]', 'y [mm]', 'S0', 'S1', 'S2', 'S3' ]
					spec_save.to_csv(save_folder+"flux_density_distrib_ef" + save_n + ".csv",index=False, sep = ' ')

				# fig = plt.figure(figsize=(15*cm_inch, 8.5*cm_inch), dpi=150)
				fig, ax = plt.subplots(figsize=(15*cm_inch, 8.5*cm_inch), dpi=150)  
				# fig.suptitle("Flux Distribution " + name, fontsize=fontsize_max)

				# spec = spec[ (spec['y'] >= -3) & (spec['y'] <= 3) ]
				# spec = spec[ (spec['z'] >= -3) & (spec['z'] <= 3) ]

				y_data=np.unique(spec['y'].to_list())
				z_data=np.unique(spec['z'].to_list())
				fun_data=np.array(spec['S0'].to_list())
				Y_data,Z_data = np.meshgrid(z_data,y_data)

				Funs=fun_data.reshape(len(y_data),len(z_data))
				plt.tight_layout()

				cmap = plt.colormaps["plasma"]
				cmap = cmap.with_extremes(bad=cmap(0))
				pcm = ax.pcolormesh(Y_data,Z_data,Funs, cmap=cmap)
				fig.colorbar(pcm, ax=ax, label='[$\dot{N}_\gamma/(0.1\%\cdot BW\cdot mm^2 )$]\t', pad=0)                    
				ax.set_xlabel('z [mm]', fontsize=12)
				ax.set_ylabel('y [mm]', fontsize=12)
				plt.tight_layout()
				plt.savefig(save_folder+"flux_density_distrib_ef_heat" + save_n + ".png"   , bbox_inches='tight')
				# pdb.set_trace()

			return spec
		else :
			# find all files that describe data within the given energy range
			en_vals_load = []
			for ind, en_val in enumerate(en_vals) : 
				if ( en_val >= en_range[0] ) and ( en_val <= en_range[-1] ) :
					en_vals_load.append( { 'en' : en_val, 'ind' : ind } )
				if en_val > en_range[-1] : 
					break
			# if no files found, then find the energy files that lie around the wanted interval (for interpol.)
			if len(en_vals_load) < 1 :
				en_middl = ( en_range[0] + en_range[-1] ) / 2
				for ind, en_val in enumerate(en_vals) : 
					en_vals_load.append( { 'en' : en_val, 'en_diff' : abs(en_val-en_middl), 'ind' : ind } )
				en_vals_load = pd.DataFrame(en_vals_load).sort_values(by=['en_diff'])
				en_vals_load = en_vals_load.to_dict('records')
				load_d = en_vals_load[0]
				en_vals_load = [ { 'en' : load_d['en'], 'ind' : load_d['ind'] } ]

			# add so many values at both edges of energy interval to energies for interpolation
			num_edge_ind = 10
			edge_ind_low = en_vals_load[0]['ind'] - num_edge_ind
			edge_ind_up = en_vals_load[-1]['ind'] + num_edge_ind
			# add edges
			if edge_ind_low < 0 :
				edge_ind_low = 0
			if edge_ind_up >= len( en_vals ) :
				edge_ind_up = len( en_vals ) - 1

			# after the real energy interval for which we integrate is found, 
			# find the files (for real) that are loaded
			en_vals_tmp = en_vals[edge_ind_low:edge_ind_up+1]
			en_vals_load = []
			for ind, en_val in enumerate(en_vals) : 
				if ( ind >= edge_ind_low ) and ( ind <= edge_ind_up ) :
					en_vals_load.append( { 'en' : en_val, 'ind' : ind } )
				if ind > edge_ind_up : 
					break
			# load all the spectrum data over pinhole, for each pinhole point
			# create list of flux values at the different energies present
			spec_tmp = []
			for ind, en_dict in enumerate(en_vals_load) : 
				en_load = en_files[ en_dict['ind'] ]['en']
				data_load = en_files[ en_dict['ind'] ]['data'].to_dict('records')
				if len(spec_tmp) < 1 :
					for line in data_load : 
						spec_tmp.append( { 'z' : line['z'], 'y' : line['y'], 'spec_d' : [] } )
				for ind_d, line in enumerate(data_load) : 
					if power_distr :
						spec_tmp[ind_d]['spec_d'].append( { 'en' : en_load, 'S0' : en_load*line['S0'] } )
					else :
						spec_tmp[ind_d]['spec_d'].append( { 'en' : en_load, 'S0' : line['S0'] } )
			# for all points of pinhole, integrate over the energy range
			spec = []
			for ind, line in enumerate(spec_tmp) : 
				if len(xlim) > 0 :
					if ( line['z'] < xlim[0] ) or ( line['z'] > xlim[1] ) :
						continue
				if len(ylim) > 0 :
					if ( line['y'] < ylim[0] ) or ( line['y'] > ylim[1] ) :
						continue
				data_spl = line['spec_d']
				data_spl = pd.DataFrame(data_spl)
				num_ens = len(data_spl['en'])
				if num_ens < 1000 : 
					num_ens = 1000
				cs_flux = CubicSpline(data_spl['en'] , data_spl['S0'])
				integr = quad( lambda x: cs_flux(x), en_range[0], en_range[-1], limit = num_ens )[0]
				spec.append( { 'z' : line['z'], 'y' : line['y'], 'S0' : integr } )
			spec = pd.DataFrame(spec)
			name = '-E$\in$ ['+str(round(en_range[0]/1e3,3))+','+str(round(en_range[1]/1e3,3))+'] keV \n'
			save_n = '_e_'+str(en_range[0])+'_'+str(en_range[1]) + '_'
		pdb.set_trace()
		fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		if power_distr :
			fig.suptitle("Power Distribution " + name, fontsize=fontsize_max)
		else :
			fig.suptitle("Flux Distribution " + name, fontsize=fontsize_max)
		ax = fig.add_subplot(111, projection='3d')
		ax.plot_trisurf(spec['z'], spec['y'], spec['S0'], cmap=cm.jet, linewidth=0.2)
		ax.set_xlabel('z [mm]')
		ax.set_ylabel('y [mm]')
		if power_distr :
			ax.set_zlabel('\n [$P/(0.1\%\cdot BW\cdot mm^2 )$]\t',linespacing=4)
		else :
			ax.set_zlabel('\n [$\dot{N}_\gamma/(0.1\%\cdot BW\cdot mm^2 )$]\t',linespacing=4)
		if len(xlim) > 0 :
			ax.set_xlim( xlim )
		if len(ylim) > 0 :
			ax.set_ylim( ylim )
		if power_distr :
			plt.savefig(save_folder+"power_density_distrib_ef" + save_n + ".eps" , bbox_inches='tight',format='eps')
		else :
			plt.savefig(save_folder+"flux_density_distrib_ef" + save_n + ".eps" , bbox_inches='tight',format='eps')
		if save_data_w_pics :
			if power_distr :
				spec.columns = [ 'z [mm]', 'y [mm]', 'power_dens_distr [En/(s 0.1% BW mm^2) ]' ]
				spec.to_csv(save_folder+"power_density_distrib_ef" + save_n + ".csv",index=False, sep = ' ')
			else :
				spec.columns = [ 'z [mm]', 'y [mm]', 'flux_dens_distr [Nphot/(s 0.1% BW mm^2) ]' ]
				spec.to_csv(save_folder+"flux_density_distrib_ef" + save_n + ".csv",index=False, sep = ' ')

		fig = plt.figure(figsize=(13*cm_inch, 8.5*cm_inch), dpi=150)
		fig.suptitle("Flux Distribution " + name, fontsize=fontsize_max)
		ax = plt.gca()

		y_data=np.unique(spec['y'])
		z_data=np.unique(spec['z'])
		fun_data=np.array(data['S0'].to_list())
		Y_data,Z_data = np.meshgrid(z_data,y_data)

		Funs=fun_data.reshape(len(y_data),len(z_data))

		plt.tight_layout()
		cmap = plt.colormaps["plasma"]
		cmap = cmap.with_extremes(bad=cmap(0))
		pcm = ax.pcolormesh(Y_data,Z_data,Funs, cmap=cmap)
		fig.colorbar(pcm, ax=ax, label='\n [$\dot{N}_\gamma/(0.1\%\cdot BW\cdot mm^2 )$]\t', pad=0)                    
		ax.set_xlabel('z [mm]', fontsize=12)
		ax.set_ylabel('y [mm]', fontsize=12)
		plt.savefig(save_folder+"flux_density_distrib_ef_heat" + save_n + ".eps"   , bbox_inches='tight')

		if plot : 
			# plt.pause(0.25)
			# plt.ion()
			plt.show()
		return spec		