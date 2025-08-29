"""
defines quantities which hold data and implements some basic plot-functions
"""

from unduwave.unduwave_incl import *
from matplotlib.pylab import matplotlib
from matplotlib.ticker import ScalarFormatter

class quantity :
	def __init__(self,api=None,name=None,description=None,unit=None,data=None,plot_name=None) :
		"""
		api - reference to the api class (wave or undu)
		name : name of the quantity
		description: some basic description
		unit : physical unit
		data: data-object
		plot_name: name shown in plot
		"""
		self._api = api
		self._name = name
		self._description = description
		self._unit = unit
		self._data = data
		self._plot_name = plot_name
		if self._plot_name is None :
			self._plot_name = self._name

	def save_data(self,file) :
		dfd=pd.DataFrame({self._name:self._data})
		header=[self._name]
		dfd.to_csv( file , index = False, sep = ' ', header=header)

	def load_data(self,file,x_quant,y_quant=None) :
		if y_quant is None :
			dfd=pd.DataFrame({self._name:self._data,x_quant._name:x_quant._data})
			header=[self._name,x_quant._name]
		else :
			dfd=pd.DataFrame({self._name:self._data,x_quant._name:x_quant._data,y_quant._name:y_quant._data})
			header=[self._name,x_quant._name,y_quant._name]
		dfd.to_csv( file , index = False, sep = ' ', header=header )

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
		ax.set_box_aspect(aspect=None, zoom=0.8)
		if title is None :
			fig.suptitle(self._description, fontsize=14)
		else :
			fig.suptitle(title, fontsize=14)
		plt.xlabel(f'{x_quant._plot_name} [{x_quant._unit}]', fontsize=12),
		plt.ylabel(f'{y_quant._plot_name} [{y_quant._unit}]', fontsize=12),
		ax.set_zlabel(f'{self._plot_name} {self._unit}', fontsize=12)
		ax.legend(loc='best', bbox_to_anchor=(0.8, 0.5, 0.0, 0.0))  
		pics_folder='pics/'
		if file_name is None :
			if self._api is None :
				file_name = f'{self._name}_over_{x_quant._name}_parametric.png'
			else:
				pics_folder = self._api._prog_paras.res_folder.get()+self._api._prog_paras.pics_folder.get()
				file_name = f'{pics_folder}{self._name}_over_{x_quant._name}_parametric.png'
		else:
			file_name=pics_folder+file_name
		if not nosave :
			plt.savefig(file_name , bbox_inches='tight')
		if plot :
			plt.draw()
			plt.ion()
		nfig=nfig+1
		return nfig

	def plot_over(self,
			x_quant,
			title=None,
			file_name=None,
			nosave=False,
			nfig=None,
			loglog=False,
			plot=True,
			leg=True,
			clear=False,
			dataFile=None,
			xlim=None,
			ylim=None,
			addPicsFolder=True,
			) :
		"""
		Basic 2d plot of data

		x_quant - the quantity to be on the x-axis
		loglog - True if both axes to be logarithmic
		"""
		if nfig is None :
			fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		else :
			fig = plt.figure(nfig,figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		if clear:
			plt.clf()
		if loglog :
			plt.loglog(x_quant._data,self._data, '-', label = self._name)
		else :
			plt.plot(x_quant._data,self._data, '-', label = self._name)
		if title is None :
			fig.suptitle(self._description, fontsize=12, y=1.2)
		else :
			fig.suptitle(title, fontsize=12, y=1.2)
		plt.xlabel(f'{x_quant._plot_name} [{x_quant._unit}]', fontsize=10),
		plt.ylabel(f'{self._plot_name} [{self._unit}]', fontsize=10)
		ax = plt.gca()

		if not (xlim is None) :
			ax.set_xlim(xlim)
		if not (ylim is None) :
			ax.set_ylim(ylim)		

		ax.tick_params(axis='both', which='major', labelsize=8)
		ax.tick_params(axis='both', which='minor', labelsize=8)
		plt.xticks(fontsize=8)
		if leg:
			ax.legend(loc='best', bbox_to_anchor=(0.8, 0.5, 0.0, 0.0))  
		pics_folder = self._api._prog_paras.res_folder.get()+self._api._prog_paras.pics_folder.get()
		if file_name is None :
			if self._api is None :
				file_name = f'{self._name}_over_{x_quant._name}.png'
			else:
				file_name = f'{pics_folder}{self._name}_over_{x_quant._name}.png'
		else:
			if addPicsFolder:
				file_name=pics_folder+file_name
		if not nosave :
			try:
				plt.savefig(file_name , bbox_inches='tight')
			except:
				pdb.set_trace()
		if dataFile is None:
			dataFile=f'{pics_folder}{self._name}_over_{x_quant._name}.dat'
		pd.DataFrame({'x':x_quant._data,'y':self._data}).to_csv(dataFile, sep = ' ',header=['x','y'], index=False)
		if plot :
			plt.draw()
			plt.ion()
		nfig=nfig+1
		return nfig

	def plot_over_3d(self,
			x_quant,
			y_quant,
			title=None,
			file_name=None,
			nosave=False,
			nfig=None,
			plot=True,
			zLabelAdd='',
			clear=False,
			dataFile=None,
			addPicsFolder=True,
			) :
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
		leny=len(y_data)
		lenz=len(z_data)
		ynew = np.linspace(min(y_data), max(y_data), leny*4)
		znew = np.linspace(min(z_data), max(z_data), lenz*4)

		Funs_intrpltd = interpol_f(ynew,znew)
		Y_data_intrpltd,Z_data_intrpltd = np.meshgrid(ynew, znew)

		"""
		Plotting Original 3D and Heat
		"""
		yTitlePos=0.8
		if nfig is None :
			fig = plt.figure(figsize=(13*cm_inch, 8.5*cm_inch), dpi=150)
		else :
			fig = plt.figure(nfig,figsize=(13*cm_inch, 8.5*cm_inch), dpi=150)
		if clear:
			plt.clf()
		nfig=nfig+1
		if title is None :
			fig.suptitle(self._description, fontsize=12, y=yTitlePos)
		else :
			fig.suptitle(title, fontsize=12, y=yTitlePos)
		ax = fig.add_subplot(111, projection='3d')

		Funs=fun_data.reshape(len(y_data),len(z_data))
		Funs=Funs.T
		ax.plot_trisurf(Y_data.flatten(), Z_data.flatten(), Funs.flatten(), cmap=cm.jet, linewidth=0.2)
		ax.set_ylabel(f'{x_quant._plot_name} [{x_quant._unit}]', fontsize=8)
		ax.set_xlabel(f'{y_quant._plot_name} [{y_quant._unit}]', fontsize=8)
		ax.set_zlabel(zLabelAdd+f'{self._plot_name} [{self._unit}]', fontsize=8)
		ax.tick_params(axis='both', which='major', labelsize=8)
		ax.set_box_aspect(aspect=None, zoom=0.7)
		ax.zaxis.get_offset_text().set_fontsize(8)
		plt.tight_layout()

		pics_folder='pics/'
		if file_name is None :
			if self._api is None :
				file_name = f'{self._name}_over_{x_quant._name}_{y_quant._name}_3d.png'
			else:
				pics_folder = self._api._prog_paras.res_folder.get()+self._api._prog_paras.pics_folder.get()
				file_name = f'{pics_folder}{self._name}_over_{x_quant._name}_{y_quant._name}_3d.png'
		else:
			if addPicsFolder :
				file_name=pics_folder+file_name
		if not nosave :
			plt.savefig(file_name   , bbox_inches='tight')
			# if dataFile is None:
			# 	dataFile=f'{pics_folder}{self._name}_over_{x_quant._name}_{y_quant._name}_3d.dat'
			# pd.DataFrame({'x':x_quant._data,'y':y_quant._data,'z':self._data}).to_csv(dataFile, sep = ' ',header=['x','y','z'], index=False)

		if nfig is None :
			fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		else :
			fig = plt.figure(nfig,figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		if clear:
			plt.clf()
		nfig=nfig+1
		if title is None :
			fig.suptitle(self._description, fontsize=12)
		else :
			fig.suptitle(title, fontsize=12)
		ax = plt.gca()

		plt.tight_layout()
		cmap = plt.colormaps["plasma"]
		cmap = cmap.with_extremes(bad=cmap(0))
		pcm = ax.pcolormesh(Y_data,Z_data,Funs, cmap=cmap)
		# pcm = ax.pcolormesh(Z_data,Y_data,Funs, cmap=cmap)
		cb=fig.colorbar(pcm, ax=ax, label=f'{self._plot_name}\n [{self._unit}]')                    

		axCb = cb.ax
		text = axCb.yaxis.label
		axCb.yaxis.get_offset_text().set_fontsize(8)
		axCb.tick_params(axis='both', which='major', labelsize=8)
		font = matplotlib.font_manager.FontProperties(size=8)
		text.set_font_properties(font)

		ax.set_xlabel(f'{x_quant._plot_name} [{x_quant._unit}]', fontsize=8)
		ax.set_ylabel(f'{y_quant._plot_name} [{y_quant._unit}]', fontsize=8)
		ax.tick_params(axis='both', which='major', labelsize=8)
		file_name_h = file_name.split('.png')[0]+'_heat.png'
		if not nosave :
			plt.savefig(file_name_h , bbox_inches='tight')
			# if dataFile is None:
			# 	dataFile=f'{pics_folder}{self._name}_over_{x_quant._name}_{y_quant._name}_heat.dat'
			# else:
			# 	dataFile=dataFile+'_heat.dat'
			# pd.DataFrame({'x':x_quant._data,'y':y_data._data,'z':self._data}).to_csv(dataFile, sep = ' ',header=['x','y','z'], index=False)

		"""
		Plotting Interpolated 3D and Heat
		"""

		if nfig is None :
			fig = plt.figure(figsize=(1.25*13*cm_inch, 1.25*8.5*cm_inch), dpi=150)
		else :
			fig = plt.figure(nfig,figsize=(1.25*13*cm_inch, 1.25*8.5*cm_inch), dpi=150)
		if clear:
			plt.clf()
		nfig=nfig+1

		if title is None :
			fig.suptitle(self._description, fontsize=12, y=yTitlePos)
		else :
			fig.suptitle(title, fontsize=12, y=yTitlePos)
		ax = fig.add_subplot(111, projection='3d')

		Funs_intrpltd=Funs_intrpltd.T
		ax.plot_trisurf(Y_data_intrpltd.flatten(), Z_data_intrpltd.flatten(), Funs_intrpltd.flatten(), cmap=cm.jet, linewidth=0.2)
		ax.set_xlabel(f'{x_quant._plot_name} [{x_quant._unit}]', fontsize=8)
		ax.set_ylabel(f'{y_quant._plot_name} [{y_quant._unit}]', fontsize=8)
		ax.set_zlabel(zLabelAdd+f'{self._plot_name} [{self._unit}]', fontsize=8)
		ax.zaxis.get_offset_text().set_fontsize(8)
		ax.tick_params(axis='both', which='major', labelsize=8)
		plt.tight_layout()
		ax.set_box_aspect(aspect=None, zoom=0.7)
		file_name_3d_intr = file_name.split('.png')[0]+'_interpolated.png'
		if not nosave :
			plt.savefig(file_name_3d_intr   , bbox_inches='tight')
		if dataFile is None:
			dataFile3d=f'{pics_folder}{self._name}_over_{x_quant._name}_{y_quant._name}_3D_interpolated.dat'
		else:
			dataFile3d=dataFile+'_3d_interpolated.dat'
		pd.DataFrame({'x':Y_data_intrpltd.flatten(),'y':Z_data_intrpltd.flatten(),'z':Funs_intrpltd.flatten()}).to_csv(dataFile3d, sep = ' ',header=['x','y','z'], index=False)

		if nfig is None :
			fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		else :
			fig = plt.figure(nfig,figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		if clear:
			plt.clf()
		if title is None :
			fig.suptitle(self._description, fontsize=12)
		else :
			fig.suptitle(title, fontsize=12)
		ax = plt.gca()
		nfig=nfig+1

		plt.tight_layout()
		cmap = plt.colormaps["plasma"]
		cmap = cmap.with_extremes(bad=cmap(0))
		pcm = ax.pcolormesh(Y_data_intrpltd,Z_data_intrpltd,Funs_intrpltd, cmap=cmap)
		cb=fig.colorbar(pcm, ax=ax, label=f'{self._plot_name}\n [{self._unit}]', pad=0)                    
		axCb = cb.ax
		text = axCb.yaxis.label
		axCb.yaxis.get_offset_text().set_fontsize(8)
		axCb.tick_params(axis='both', which='major', labelsize=8)
		font = matplotlib.font_manager.FontProperties(size=8)
		text.set_font_properties(font)

		ax.set_xlabel(f'{x_quant._plot_name} [{x_quant._unit}]', fontsize=8)
		ax.set_ylabel(f'{y_quant._plot_name} [{y_quant._unit}]', fontsize=8)
		ax.tick_params(axis='both', which='major', labelsize=8)
		file_name_heat_intr = file_name.split('.png')[0]+'_heat_interpolated.png'
		if not nosave :
			plt.savefig(file_name_heat_intr , bbox_inches='tight')
		if dataFile is None:
			dataFile=f'{pics_folder}{self._name}_over_{x_quant._name}_{y_quant._name}_heat_interpolated.dat'
		else:
			dataFile=dataFile+'_heat_interpolated.dat'
		pd.DataFrame({'x':Y_data_intrpltd.flatten(),'y':Z_data_intrpltd.flatten(),'z':Funs_intrpltd.flatten()}).to_csv(dataFile, sep = ' ',header=['x','y','z'], index=False)

		if plot :
			plt.draw()
			plt.ion()
		return nfig

	def save_plot(self) :
		"""
		"""
		pass