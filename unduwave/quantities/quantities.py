"""
defines quantities which hold data and implements some basic plot-functions
"""

from unduwave.unduwave_incl import *

class wave_quantity :
	def __init__(self,wave_api=None,name=None,description=None,unit=None,data=None,plot_name=None) :
		"""
		wave_api - reference to the wave_api class
		name : name of the quantity
		description: some basic description
		unit : physical unit
		data: data-object
		plot_name: name shown in plot
		"""
		self._wave_api = wave_api
		self._name = name
		self._description = description
		self._unit = unit
		self._data = data
		self._plot_name = plot_name
		if self._plot_name is None :
			self._plot_name = self._name

	def plot_parametric_3d(self,x_quant,y_quant,title=None,file_name=None,nosave=False,nfig=None,plot=True) :
		"""
		Basic plot of data, 2d, 3d

		x_quant - the quantity plotted on the x-axis
		y_quant - quantity for y-axis plot
		title - title, if none is taken from description of this quantity

		Draws a parametetric curve (x_quant,y_quant,self(x_quant,y_quant))
		"""
		if nfig is None :
			fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		else :
			fig = plt.figure(nfig,figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		ax = fig.add_subplot(projection='3d')
		ax.plot(x_quant._data,y_quant._data,self._data, '-', label = self._name)
		if title is None :
			fig.suptitle(self._description, fontsize=14)
		else :
			fig.suptitle(title, fontsize=14)
		plt.xlabel(f'{x_quant._plot_name} [{x_quant._unit}]', fontsize=12),
		plt.ylabel(f'{y_quant._plot_name} [{y_quant._unit}]', fontsize=12),
		ax.set_zlabel(f'{self._plot_name} [{self._unit}]', fontsize=12)
		ax.legend(loc='best', bbox_to_anchor=(0.8, 0.5, 0.0, 0.0))  
		if file_name is None :
			if self._wave_api is None :
				file_name = f'{self._name}_over_{x_quant._name}_parametric.png'
			else:
				pics_folder = self._wave_api._wave_prog_paras.res_folder.get()+self._wave_api._wave_prog_paras.pics_folder.get()
				file_name = f'{pics_folder}{self._name}_over_{x_quant._name}_parametric.png'
		if not nosave :
			plt.savefig(file_name , bbox_inches='tight')
		if plot :
			plt.draw()
			plt.ion()
		nfig=nfig+1
		return nfig

	def plot_over(self,x_quant,title=None,file_name=None,nosave=False,nfig=None,loglog=False,plot=True) :
		"""
		Basic 2d plot of data

		x_quant - the quantity to be on the x-axis
		loglog - True if both axes to be logarithmic
		"""
		if nfig is None :
			fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		else :
			fig = plt.figure(nfig,figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		if loglog :
			plt.loglog(x_quant._data,self._data, '-', label = self._name)
		else :
			plt.plot(x_quant._data,self._data, '-', label = self._name)
		if title is None :
			fig.suptitle(self._description, fontsize=14)
		else :
			fig.suptitle(title, fontsize=14)
		plt.xlabel(f'{x_quant._plot_name} [{x_quant._unit}]', fontsize=12),
		plt.ylabel(f'{self._plot_name} [{self._unit}]', fontsize=12)
		ax = plt.gca()
		ax.legend(loc='best', bbox_to_anchor=(0.8, 0.5, 0.0, 0.0))  
		if file_name is None :
			if self._wave_api is None :
				file_name = f'{self._name}_over_{x_quant._name}.png'
			else:
				pics_folder = self._wave_api._wave_prog_paras.res_folder.get()+self._wave_api._wave_prog_paras.pics_folder.get()
				file_name = f'{pics_folder}{self._name}_over_{x_quant._name}.png'
		if not nosave :
			plt.savefig(file_name , bbox_inches='tight')
		if plot :
			plt.draw()
			plt.ion()
		nfig=nfig+1
		return nfig

	def plot_over_3d(self,x_quant,y_quant,title=None,file_name=None,nosave=False,nfig=None,plot=True) :
		"""
		Plots 3D plots and heat plots of data
		x_quant : quantity to be used as x-data
		y_quant : quantity to be used as y-data
		file_name : name under which to save plot
		nosave: True if you do not want to save plot

		Creates a 3D plot and a heat plot of the original data and of interpolated data.
		"""

		"""
		Original Data
		"""
		y_data=copy.deepcopy(np.unique(x_quant._data))
		y_data.sort()
		z_data=copy.deepcopy(np.unique(y_quant._data))
		z_data.sort()
		fun_data=copy.deepcopy(np.array(self._data))
		Y_data,Z_data = np.meshgrid(z_data,y_data)

		"""
		Creating interpolated data
		"""

		Y_data_intr = np.zeros( (len(y_data), len(z_data)) )
		Z_data_intr = np.zeros( (len(y_data), len(z_data)) )
		Fun_data_intr = np.zeros( (len(y_data), len(z_data)) )

		for ind, val in enumerate(fun_data.tolist()) :
			y_val = x_quant._data[ind]
			z_val = y_quant._data[ind]
			ind_g = 0
			ind_s = 0
			for ind_gt, gap_f in enumerate(y_data):
				if gap_f == y_val :
					ind_g = ind_gt
					break
			for ind_st, shift_f in enumerate(z_data):
				if shift_f == z_val :
					ind_s = ind_st
					break
			Y_data_intr[ind_g,ind_s] = y_val
			Z_data_intr[ind_g,ind_s] = z_val
			Fun_data_intr[ind_g,ind_s] = val

		interpol_f = RectBivariateSpline(y_data,z_data, Fun_data_intr )
		ynew = np.linspace(min(y_data), max(y_data), 150)
		znew = np.linspace(min(z_data), max(z_data), 150)

		Funs_intrpltd = interpol_f(ynew,znew)
		Y_data_intrpltd,Z_data_intrpltd = np.meshgrid(ynew, znew)

		"""
		Plotting Original 3D and Heat
		"""

		if nfig is None :
			fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		else :
			fig = plt.figure(nfig,figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		nfig=nfig+1
		if title is None :
			fig.suptitle(self._description, fontsize=14)
		else :
			fig.suptitle(title, fontsize=14)
		ax = fig.add_subplot(111, projection='3d')

		Funs=fun_data.reshape(len(y_data),len(z_data))

		ax.plot_trisurf(Y_data.flatten(), Z_data.flatten(), Funs.flatten(), cmap=cm.jet, linewidth=0.2)
		ax.set_xlabel(f'{x_quant._plot_name} [{x_quant._unit}]', fontsize=12)
		ax.set_ylabel(f'{y_quant._plot_name} [{y_quant._unit}]', fontsize=12)
		ax.set_zlabel(f'{self._plot_name} [{self._unit}]', fontsize=12)
		plt.tight_layout()

		if file_name is None :
			if self._wave_api is None :
				file_name = f'{self._name}_over_{x_quant._name}_{y_quant._name}_3d.png'
			else:
				pics_folder = self._wave_api._wave_prog_paras.res_folder.get()+self._wave_api._wave_prog_paras.pics_folder.get()
				file_name = f'{pics_folder}{self._name}_over_{x_quant._name}_{y_quant._name}_3d.png'
		if not nosave :
			plt.savefig(file_name   , bbox_inches='tight')

		if nfig is None :
			fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		else :
			fig = plt.figure(nfig,figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		nfig=nfig+1
		if title is None :
			fig.suptitle(self._description, fontsize=14)
		else :
			fig.suptitle(title, fontsize=14)
		ax = plt.gca()

		plt.tight_layout()
		cmap = plt.colormaps["plasma"]
		cmap = cmap.with_extremes(bad=cmap(0))
		pcm = ax.pcolormesh(Y_data,Z_data,Funs, cmap=cmap)
		fig.colorbar(pcm, ax=ax, label=f'{self._plot_name} [{self._unit}]', pad=0)                    
		ax.set_xlabel(f'{x_quant._plot_name} [{x_quant._unit}]', fontsize=12)
		ax.set_ylabel(f'{y_quant._plot_name} [{y_quant._unit}]', fontsize=12)
		file_name_h = file_name.split('.png')[0]+'_heat.png'
		if not nosave :
			plt.savefig(file_name_h , bbox_inches='tight')

		"""
		Plotting Interpolated 3D and Heat
		"""

		if nfig is None :
			fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		else :
			fig = plt.figure(nfig,figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		nfig=nfig+1

		if title is None :
			fig.suptitle(self._description, fontsize=14)
		else :
			fig.suptitle(title, fontsize=14)
		ax = fig.add_subplot(111, projection='3d')

		ax.plot_trisurf(Y_data_intrpltd.flatten(), Z_data_intrpltd.flatten(), Funs_intrpltd.flatten(), cmap=cm.jet, linewidth=0.2)
		ax.set_xlabel(f'{x_quant._plot_name} [{x_quant._unit}]', fontsize=12)
		ax.set_ylabel(f'{y_quant._plot_name} [{y_quant._unit}]', fontsize=12)
		ax.set_zlabel(f'{self._plot_name} [{self._unit}]', fontsize=12)
		plt.tight_layout()
		file_name_3d_intr = file_name.split('.png')[0]+'_interpolated.png'
		if not nosave :
			plt.savefig(file_name_3d_intr   , bbox_inches='tight')

		if nfig is None :
			fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		else :
			fig = plt.figure(nfig,figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		if title is None :
			fig.suptitle(self._description, fontsize=14)
		else :
			fig.suptitle(title, fontsize=14)
		ax = plt.gca()
		nfig=nfig+1

		plt.tight_layout()
		cmap = plt.colormaps["plasma"]
		cmap = cmap.with_extremes(bad=cmap(0))
		pcm = ax.pcolormesh(Y_data_intrpltd,Z_data_intrpltd,Funs_intrpltd, cmap=cmap)
		fig.colorbar(pcm, ax=ax, label=f'{self._plot_name} [{self._unit}]', pad=0)                    
		ax.set_xlabel(f'{x_quant._plot_name} [{x_quant._unit}]', fontsize=12)
		ax.set_ylabel(f'{y_quant._plot_name} [{y_quant._unit}]', fontsize=12)
		file_name_heat_intr = file_name.split('.png')[0]+'_heat_interpolated.png'
		if not nosave :
			plt.savefig(file_name_heat_intr , bbox_inches='tight')

		if plot :
			plt.draw()
			plt.ion()
		return nfig

	def save_plot(self) :
		"""
		"""
		pass