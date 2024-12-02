from unduwave.unduwave_incl import *

class wave_quantity :
	def __init__(self,name=None,description=None,unit=None,data=None,plot_name=None) :
		self._name = name
		self._description = description
		self._unit = unit
		self._data = data
		self._plot_name = plot_name
		if self._plot_name is None :
			self._plot_name = self._name

	def plot_over(self,x_quant,file_name=None,nosave=False) :
		"""
		Basic plot of data, 2d, 3d
		"""
		fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		plt.plot(x_quant._data,self._data, '-', label = self._name)
		fig.suptitle(self._description, fontsize=14)
		plt.xlabel(f'{x_quant._name} [{x_quant._unit}]', fontsize=12),
		plt.ylabel(f'{self._plot_name} [{self._unit}]', fontsize=12)
		ax = plt.gca()
		ax.legend(loc='best', bbox_to_anchor=(0.8, 0.5, 0.0, 0.0))  
		if file_name is None :
			file_name = f'{self._name}_over_{x_quant._name}.png'
		if not nosave :
			plt.savefig(file_name , bbox_inches='tight')
		plt.show()

	def plot_over_3d(self,x_quant,y_quant,file_name=None,nosave=False) :
		fig = plt.figure(figsize=(13*cm_inch, 8.5*cm_inch), dpi=150)

		if file_name is None :
			file_name = f'{self._name}_over_{x_quant._name}_{y_quant._name}_3d.png'

		fig.suptitle(self._description, fontsize=14)
		ax = fig.add_subplot(111, projection='3d')

		y_data=np.unique(x_quant._data)
		z_data=np.unique(y_quant._data)
		fun_data=np.array(self._data)
		Y_data,Z_data = np.meshgrid(z_data,y_data)

		Funs=fun_data.reshape(len(y_data),len(z_data))

		ax.plot_trisurf(Y_data.flatten(), Z_data.flatten(), Funs.flatten(), cmap=cm.jet, linewidth=0.2)

		# if len(xlim) > 0 :
		# 	ax.set_xlim( xlim )
		# if len(ylim) > 0 :
		# 	ax.set_ylim( ylim )
		ax.set_xlabel(f'{x_quant._name} [{x_quant._unit}]', fontsize=12)
		ax.set_ylabel(f'{y_quant._name} [{y_quant._unit}]', fontsize=12)
		ax.set_zlabel(f'{self._plot_name} [{self._unit}]', fontsize=12)
		plt.tight_layout()

		plt.savefig(file_name   , bbox_inches='tight')
		# if save_data_w_pics :
		# 	data_save=copy.deepcopy(data)
		# 	data_save.columns = [ 'z [mm]', 'y [mm]', 'power [W/mm^2]' ]
		# 	data_save.to_csv(res_folder_pics+"power.csv",index=False, sep = ' ')

		fig = plt.figure(figsize=(13*cm_inch, 8.5*cm_inch), dpi=150)
		fig.suptitle(self._description, fontsize=14)
		ax = plt.gca()

		plt.tight_layout()
		cmap = plt.colormaps["plasma"]
		cmap = cmap.with_extremes(bad=cmap(0))
		pcm = ax.pcolormesh(Y_data,Z_data,Funs, cmap=cmap)
		fig.colorbar(pcm, ax=ax, label=f'{self._plot_name} [{self._unit}]', pad=0)                    
		ax.set_xlabel(f'{x_quant._name} [{x_quant._unit}]', fontsize=12)
		ax.set_ylabel(f'{y_quant._name} [{y_quant._unit}]', fontsize=12)
		file_name = file_name.split('.png')[0]+'_heat.png'
		plt.savefig(file_name , bbox_inches='tight')
		plt.show()

	def save_plot(self) :
		"""
		"""
		pass