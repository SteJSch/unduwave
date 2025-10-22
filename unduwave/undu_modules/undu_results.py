"""
This file is nice
"""

from unduwave.unduwave_incl import *
from unduwave.undu_modules.undu_parameters import *
import unduwave.helpers.file_folder_helpers as f_h
from unduwave.quantities.quantities import *

class undu_results : 

	def __init__(self,undu_api: int) :
		"""
		This class is for this and that

		:param parameter1: my favorite int
		:param parameter2: my favorite class
		"""
		self._undu_api = undu_api
		self._res_quantities = []
		self._summary=None
		self._res_folder      = self._undu_api._prog_paras.res_folder.get()
		self._res_folder_undu = self._res_folder + self._undu_api._prog_paras.undu_data_res_folder.get()
		self._res_folder_pics = self._res_folder + self._undu_api._prog_paras.pics_folder.get()
		if not os.path.exists(self._res_folder_pics):
			os.makedirs(self._res_folder_pics)

	def load_from_res_folder(self) :
		"""
		This fun is for this and that

		:param parameter1: my favorite int
		:param parameter2: my favorite class
		:return: returns some nothing
		"""

		self.load_on_axis_undumag_file()
		self.load_field_profile()
		self.load_field_map()
		self.extract_summary()

	def load_field_map(self) : 

		file_name = 'undumag.map'
		files = f_h.find_files_exptn(folder = self._res_folder_undu, hints = [file_name], exptns = [])

		data = pd.read_csv( self._res_folder_undu+files[0], skiprows=range(0, 3), dtype=object, delimiter=r"\s+",header=None)
		data.columns = [ 'imoth', 'imag', 'mat', 'ityp', 'matmod', 'x', 'y', 'z', 'Bx', 'By', 'Bz', 'B', 'Hx', 'Hy', 'Hz', 'H', 'Mx', 'My', 'Mz', 'M', 'BxDip', 'ByDip', 'BzDip', 'ifail', 'kfail', 'cmag', 'cmoth' ]
		cols_float = ['x', 'y', 'z', 'Bx', 'By', 'Bz']
		for col in cols_float:
			data[col] = data[col].astype(float)

		bmap = quantity(
			api=self._undu_api,
			data=data.to_dict('records'),
			plot_name="bmap",
			name='bmap',
			description='magnetic field map',
			unit='',
			)
		self._res_quantities = self._res_quantities + [bmap]

	def load_field_profile(self) :
		file_name = 'undumag_field_profile'
		files = f_h.find_files_exptn(folder = self._res_folder_undu, hints = [file_name], exptns = [])
		if len(files) < 1 :
			return
		try :
			data = pd.read_csv( self._res_folder_undu+files[0], skiprows=3, dtype=object, delimiter=r"\s+")
		except: 
			return
		data.columns = [ 'x', 'y', 'z', 'Bx', 'By', 'Bz' ]
		cols_float = data.columns
		for col in cols_float:
			data[col] = data[col].astype(float)

		z_profile = quantity(
			api=self._undu_api,
			data=data['z'].to_list(),
			plot_name="z",
			name='profz',
			description='z-vertical position',
			unit='mm',
			)

		by_profile = quantity(
			api=self._undu_api,
			data=data['By'].to_list(),
			plot_name="B$_y$",
			name='profBy',
			description='By-along z',
			unit='T',
			)

		bz_profile = quantity(
			api=self._undu_api,
			data=data['Bz'].to_list(),
			plot_name="B$_z$",
			name='profBz',
			description='Bz-along z',
			unit='T',
			)
		self._res_quantities = self._res_quantities + [z_profile,by_profile,bz_profile]

	def load_on_axis_undumag_file(self) : 
		file_name = 'undumag_on-axis.dat'
		files = f_h.find_files_exptn(folder = self._res_folder_undu, hints = [file_name], exptns = [])
		if len(files) < 1 :
			return
		data = pd.read_csv( self._res_folder_undu+files[0], dtype=object, delimiter=r"\s+")
		data.columns = ['x','By','Bz','intBy','intBz','int2By','int2Bz','quark']
		cols = data.columns
		for col in cols:
			data[col] = data[col].astype(float)

		x_quant = quantity(
			api=self._undu_api,
			data=data['x'].to_list(),
			plot_name="x",
			name='trajx',
			description='x-Longitudinal Position',
			unit='mm',
			)

		by_quant = quantity(
			api=self._undu_api,
			data=data['By'].to_list(),
			plot_name="B$_y$",
			name='by',
			description='By-vertical induction',
			unit='T',
			)

		bz_quant = quantity(
			api=self._undu_api,
			data=data['Bz'].to_list(),
			plot_name="B$_z$",
			name='bz',
			description='Bz-horizontal induction',
			unit='T',
			)

		byInt_quant = quantity(
			api=self._undu_api,
			data=data['intBy'].to_list(),
			plot_name="1. Integral B$_y$",
			name='intBy',
			description='first integral By',
			unit='Tmm',
			)

		bzInt_quant = quantity(
			api=self._undu_api,
			data=data['intBz'].to_list(),
			plot_name="1. Integral B$_z$",
			name='intBz',
			description='first integral Bz',
			unit='Tmm',
			)

		byInt2_quant = quantity(
			api=self._undu_api,
			data=data['int2By'].to_list(),
			plot_name="2. Integral B$_y$",
			name='intBy2',
			description='second integral By',
			unit='Tmm$^2$',
			)

		bzInt2_quant = quantity(
			api=self._undu_api,
			data=data['int2Bz'].to_list(),
			plot_name="2. Integral B$_z$",
			name='intBz2',
			description='second integral Bz',
			unit='Tmm$^2$',
			)
		self._res_quantities = self._res_quantities + [x_quant,by_quant,bz_quant,byInt_quant,bzInt_quant,byInt2_quant,bzInt2_quant] 

	def load_force_undumag_file(self) : 
		file_name = 'undumag.frc'
		files = f_h.find_files_exptn(folder = self._res_folder_undu, hints = [file_name], exptns = [])
		if len(files) < 1 :
			return
		with open(self._res_folder_undu+files[-1], 'r') as o_f:
			# read an store all lines into list
			out_file_lines_tmp = o_f.readlines()
		for ind, line in enumerate(out_file_lines_tmp) : 
			if line.find(' * Fx, Fy, Fz [N]') >= 0 :
				text_forces = out_file_lines_tmp[ind+1]
			if line.find(' * Tx, Ty, Tz [Nmm]') >= 0 :
				text_torques = out_file_lines_tmp[ind+1]
		force_components = []
		torque_components = []
		text_forces = text_forces.split('.') # gives 4 elements, each number has 3 after comma
		num1 = float(text_forces[0].strip() + '.' + text_forces[1].strip()[0:3])
		num2 = float(text_forces[1].strip()[3:] + '.' + text_forces[2].strip()[0:3])
		num3 = float(text_forces[2].strip()[3:] + '.' + text_forces[3].strip()[0:3])
		force_components = [num1, num2, num3]
	
		# text_torques = text_torques.split('.') # gives 4 elements, each number has 3 after comma
		# num1 = float(text_torques[0].strip() + text_torques[1].strip()[0:3])
		# num2 = float(text_torques[1].strip()[3:] + text_torques[2].strip()[0:3])
		# num3 = float(text_torques[2].strip()[3:] + text_torques[3].strip()[0:3])
		# torque_components = [num1, num2, num3]

		return force_components, torque_components

	def extract_summary(self) : 
		file_name = 'undumag.beff'
		files = f_h.find_files_exptn(folder = self._res_folder_undu, hints = [file_name], exptns = [])
		if len(files) < 1 :
			return
		# loop all resulting files we want to save
		out_file_lines_tmp = []
		with open(self._res_folder_undu+files[-1], 'r') as o_f:
			# read an store all lines into list
			out_file_lines_tmp = o_f.readlines()
		ind_f = -1
		ind_f_byint1 = -1
		ind_f_by = -1
		ind_f_bz = -1
		# search for the right line in file
		for ind_l, line in enumerate(out_file_lines_tmp) :
			if line.find('* Beff = Sqrt( ByEff**2 + BzEff**2 ), Keff, 1. Harm. [eV]:') >= 0 :
				ind_f = ind_l + 1
			if line.find('* ByInt1, BzInt1, ByInt2, BzInt2:') >= 0 :
				ind_f_byint1 = ind_l + 1
			if line.find('* ByMin, ByMax, (ByMax-ByMin)/2, ByEff:') >= 0 :
				ind_f_by = ind_l + 1
			if line.find('* BzMin, BzMax, (BzMax-BzMin)/2, BzEff:') >= 0 :
				ind_f_bz = ind_l + 1
			if ( ind_f >= 0 ) and ( ind_f_byint1 >= 0 ) and ( ind_f_by >= 0 ) and ( ind_f_bz >= 0 ) :
				break
		res = {}
		# line found?
		if ind_f >= 0 :
			# getting rid of many spaces
			line_vals = out_file_lines_tmp[ind_f]
			parts = line_vals.split(' ')
			vals = []
			for elem in parts :
				if (len(elem) > 0) and not (elem == '\n')  :
					vals.append(float(elem))
			res.update( { 'Beff' : vals[0], 'Keff' : vals[1], 'firstHarm' : vals[2] } )
		if ind_f_byint1 >= 0 :
			# getting rid of many spaces
			line_vals = out_file_lines_tmp[ind_f_byint1]
			parts = line_vals.split(' ')
			vals = []
			for elem in parts :
				if (len(elem) > 0) and not (elem == '\n')  :
					vals.append(float(elem))
			res.update( { 'ByInt1' : vals[0], 'ByInt2' : vals[2],'BzInt1' : vals[1], 'BzInt2' : vals[3] } )
		if ind_f_by >= 0 :
			# getting rid of many spaces
			line_vals = out_file_lines_tmp[ind_f_by]
			parts = line_vals.split(' ')
			vals = []
			for elem in parts :
				if (len(elem) > 0) and not (elem == '\n')  :
					vals.append(float(elem))
			res.update( { 'ByMin' : vals[0], 'ByMax' : vals[1] , 'ByEff' : vals[3] } )
		if ind_f_bz >= 0 :
			# getting rid of many spaces
			line_vals = out_file_lines_tmp[ind_f_bz]
			parts = line_vals.split(' ')
			vals = []
			for elem in parts :
				if (len(elem) > 0) and not (elem == '\n')  :
					vals.append(float(elem))
			res.update( { 'BzMin' : vals[0], 'BzMax' : vals[1] , 'BzEff' : vals[3] } )

		if not ( self._undu_api._prog_paras.calc_force.get() == 0 ) :
			force_components, torque_components = self.load_force_undumag_file()
			res.update( { 'force' : force_components, 'torque' : torque_components } )

		self._summary=res

		dataFile=f'{self._undu_api._prog_paras.res_summary_file.get()}'
		with open( self._undu_api._prog_paras.res_folder.get()+dataFile, 'w') as o_f:
			for key, val in res.items() :
				o_f.write( key + ' : ' + str(val) + '\n' )

	def get_result(self,which) :
		"""
		Pass string "flux_dens",... to get object containing data and some functionality
		"""
		for quant in self._res_quantities :
			if quant._name == which :
				return quant
