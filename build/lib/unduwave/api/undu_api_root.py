"""
The basic api
"""
from unduwave.unduwave_incl import *
from unduwave.wave_modules.wave_parameters import *
from unduwave.wave_modules.wave_prepare import *
from unduwave.wave_modules.wave_control import *
from unduwave.wave_modules.wave_postprocess import *
from unduwave.wave_modules.wave_results import *

def create_magnetization_unit_vec(magn_string) : 
	magn_unit_vec = point_coords(x=1.0,y=0.0,z=0.0)
	if magn_string == '+x' : 
		magn_unit_vec = point_coords(x=1.0,y=0.0,z=0.0)
	elif magn_string == '-x' : 
		magn_unit_vec = point_coords(x=-1.0,y=0.0,z=0.0)
	if magn_string == '+y' : 
		magn_unit_vec = point_coords(x=0.0,y=1.0,z=0.0)
	elif magn_string == '-y' : 
		magn_unit_vec = point_coords(x=0.0,y=-1.0,z=0.0)
	if magn_string == '+z' : 
		magn_unit_vec = point_coords(x=0.0,y=0.0,z=1.0)
	elif magn_string == '-z' : 
		magn_unit_vec = point_coords(x=0.0,y=0.0,z=-1.0)
	return magn_unit_vec

class point_coords : 

	def __init__(self,x=0.0,y=0.0,z=0.0) : 
		self._x = x
		self._y = y
		self._z = z

	def __sub__(self,pnt) :
		point = point_coords()
		point._x = self._x-pnt._x
		point._y = self._y-pnt._y
		point._z = self._z-pnt._z
		return point

	def __add__(self,pnt) :
		point = point_coords()
		point._x = self._x+pnt._x
		point._y = self._y+pnt._y
		point._z = self._z+pnt._z
		return point

class undu_api :
	def __init__(self,undu_mode='undu_easy',res_folder='') :
		"""
		"""

		# self._undu_prog_paras = undu_prog_parameters()
		# self._wave_prog_paras.get_std_paras(undu_mode=undu_mode)

		# self._undu_paras = undu_paras()
		# self._undu_paras.get_std_paras(undu_mode=undu_mode)

		self.nam_d = {}
		self.clc_d = {}
		self.para = self.get_std_para()
		self.para['res_folder'] = res_folder
		self._undu_prog_paras = self.para
		self._elements=[]

	def add_element(self,element) : 
		self._elements.append(element)

	def get_para(self) : 
		return self.para

	def set_para(self, para) : 
		self.para = dict(para)
		# check if something has str. val, if there is another para with that name
		# if yes, take its value
		for key, val in self.para.items() : 
			if isinstance(val,str) : 
				for key2, val2 in self.para.items() : 
					if val == key2 : 
						self.para[key] = self.para[key2]
						break
		if 'para_fun' in self.para.keys() : 
			self.para = self.para['para_fun'](self.para)

	# change the input file
	def create_fresh_nam( self,file_lines, undu_paras = None ) :
		if undu_paras is None : 
			undu_paras = self.para
		perL = undu_paras['periodL']
		numP = undu_paras['periodN']
		halfL = round(0.5*(numP+20)*perL*1.1,2) # +2 poles for endstructure and enlarge 10%
		nam_items = copy.deepcopy(self.nam_d)
		for ind, line in enumerate(file_lines) :
			if undu_paras['xmapmin'] is None :
				if line.find('xmapmin=') >= 0 :
					file_lines[ind] = ' xmapmin=' + str(-halfL) + ' ! xmin for field map\n'
					continue
				if line.find('xmapmax=') >= 0 :
					file_lines[ind] = ' xmapmax=' + str(halfL) + ' ! xmax for field map\n'
					continue
			for key_nam, val_nam in nam_items.items() : 
				if not (undu_paras[val_nam] is None) :
					if line.find(f'{key_nam}=') >= 0 :
						line_comm = ''
						if line.find('!') >= 0 :
							line_comm = ' ! ' + line.split('!')[-1]
						file_lines[ind] = f' {key_nam}=' + str(undu_paras[val_nam]) + line_comm+'\n'
						del nam_items[key_nam]
						break

	# change the input file
	def create_fresh_clc(self, file_lines, undu_paras = None ) :
		if undu_paras is None : 
			undu_paras = self.para

		no_lines = False

		for ind, line in enumerate(file_lines) :

			for key_nam, val_nam in self.clc_d.items() : 
				if line.find(f'${key_nam} =') >= 0 :
					line_comm = ''
					if line.find('!') >= 0 :
						line_comm = ' ! ' + line.split('!')[-1]
					file_lines[ind] = f'${key_nam} = ' + str(undu_paras[val_nam]) + line_comm+'\n'
					break

	# magnet length is not(!) part of parameters but is the dependent quantity
	def get_std_para(self) : 
		dir_path = os.path.dirname(os.path.realpath(__file__))
		paras = {}
		paras.update( { 'undumag_prog_fold' : dir_path+'/../../External-Software/Undumag/' } )
		paras.update( { 'input_folder' : dir_path+'/../../unduwave/UNDWAVE_IN_FILES/Undu-In-Files/' } )
		paras.update( { 'clc' : 'undu_tmp.clc', 'nam' : 'undumag.nam' } )
		paras.update( { 'res_folder' : 'res' } )
		paras.update( { 'res_file' : 'res.txt' } )
		# paras.update( { 'undu_res_files_save' : [ 'undumag.beff' ] } 
		paras.update( { 'undu_res_files_save' : \
				[ 'undumag_on-axis.dat', 'undumag.beff', 'undumag_field_profile.dat', 'undumag.clc',\
				'urad_traxyz.dat','undumag_trajectory.eps', 'undumag_y_z.eps', 'undumag_on-axis_by_bz.eps',\
				 'undumag.eps', 'undumag_on-axis_byint_bzint.eps', 'undumag_on-axis_byint2_bzint2.eps' ] } )
		paras.update( { 'periodN' : 1 } )
		paras.update( { 'nuthreads' : 6 } )
		paras.update( { 'writeGeo' : 0 } )
		paras.update( { 'plotGeo' : 1 } )
		paras.update( { 'nxbeff' : 101 } )
		paras.update( { 'ixsym' : 0 } )
		paras.update( { 'iysym' : 0 } )
		paras.update( { 'izsym' : 0 } )

		paras.update( { 'dxmap' : 1.0 } )
		paras.update( { 'periodL' : 20 } )
		paras.update( { 'resiron' : 1.0e-4, 'hconv' : 1.0e-6 } )

		paras.update( { 'zmapmin' : -30 } )
		paras.update( { 'zmapmax' : 30 } )
		paras.update( { 'nzmap' : 1 } )

		paras.update( { 'xmapmin' : None } )
		paras.update( { 'xmapmax' : None } )
		paras.update( { 'nxmap' : 1 } )
		paras.update( { 'dxmap' : 1 } )

		paras.update( { 'ymapmin' : -30 } )
		paras.update( { 'ymapmax' : 30 } )
		paras.update( { 'nymap' : 1 } )

		paras.update( { 'kxcenter' : 1 } )

		paras.update( { 'iforce' : 0 } )
		paras.update( { 'iplforce' : 0 } )
		paras.update( { 'ubfcenx' : 0.0 } )
		paras.update( { 'ubfceny' : 0.0 } )
		paras.update( { 'ubfcenz' : 0.0 } )
		paras.update( { 'ubflenx' : 100.0 } )
		paras.update( { 'ubfleny' : 100.0 } )
		paras.update( { 'ubflenz' : 100.0 } )
		paras.update( { 'perlen' : '9999.' } )

		paras.update( { 'kurad' : 0 } )
		paras.update( { 'ndivfboxy' : 1 } )
		paras.update( { 'kpreset' : 0 } )
		paras.update( { 'matrix' : 2 } )

		paras.update( { 'mbforcex' : 10 } )
		paras.update( { 'mbforcey' : 10 } )
		paras.update( { 'mbforcez' : 10 } )

		paras.update( { 'knomagmap' : 0 } )
		paras.update( { 'knopolmap' : 0 } )

		self.nam_d.update({ 'iundugeo' : 'writeGeo', 'iunduplot' : 'plotGeo', 'resiron' : 'resiron', 'hconv' : 'hconv',\
			'nxbeff' : 'nxbeff', 'izsym' : 'izsym', 'ixsym' : 'ixsym', 'iysym' : 'iysym', 'dxmap' : 'dxmap', 'nuthreads' : 'nuthreads',\
			'zmapmin' : 'zmapmin', 'zmapmax' : 'zmapmax', 'nzmap' : 'nzmap',
			'xmapmin' : 'xmapmin', 'xmapmax' : 'xmapmax', 'nxmap' : 'nxmap', 'dxmap' : 'dxmap',
			'ymapmin' : 'ymapmin', 'ymapmax' : 'ymapmax', 'nymap' : 'nymap',
			'iforce':'iforce',
			'iplforce':'iplforce', 'ubfcenx' :'ubfcenx', 'ubfceny':'ubfceny','ubfcenz':'ubfcenz',
			'ubflenx':'ubflenx','ubfleny':'ubfleny','ubflenz':'ubflenz', 'mbforcex' : 'mbforcex', 'mbforcey' : 'mbforcey',
			'mbforcez' : 'mbforcez', 'kxcenter' : 'kxcenter', 'knomagmap' : 'knomagmap', 'knopolmap' : 'knopolmap', 
			'perlen' : 'perlen', 'kurad' : 'kurad', 'ndivfboxy' : 'ndivfboxy', 
			'kpreset' : 'kpreset', 'matrix' : 'matrix' })
		self.clc_d.update({ 'nPeriods' : 'periodN', 'PerLen' : 'periodL' })

		return paras

	# if copy >= 0: result files from undu_res_files_save are saved with added the val copy added to name for distinction
	def run( self, undu_paras = None, copy = '', freate_fresh_clc = True ) : 

		load_clc_file = self.para['input_folder']+'undu_raw.clc'

		with open(load_clc_file, 'r') as o_f:
			load_clc = o_f.readlines()

		for el in self._elements :
			el.add_to_clc(
				clc_txt=load_clc,
				magns_ignore=None
				)

		with open( self.para['input_folder'] + 'undu_tmp.clc', 'w') as o_f:
			for ind, line in enumerate(load_clc) :
				o_f.write(line)

		if undu_paras is None : 
			undu_paras = self.para
		undu_folder = undu_paras['undumag_prog_fold']
		inp_folder = undu_paras['input_folder']
		configFile_clc = undu_paras['clc']
		configFile_nam = undu_paras['nam']
		res_folder = undu_paras['res_folder']
		undu_res_files_save = undu_paras['undu_res_files_save']
		# open the configuration file
		out_file_lines_clc = []

		with open(inp_folder+configFile_clc, 'r') as o_f:
			# read an store all lines into list
			out_file_lines_clc = o_f.readlines()
		out_file_lines_nam = []
		with open(inp_folder+configFile_nam, 'r') as o_f:
			# read an store all lines into list
			out_file_lines_nam = o_f.readlines()
		if freate_fresh_clc :
			self.create_fresh_clc( file_lines = out_file_lines_clc, undu_paras = undu_paras )
		self.create_fresh_nam( file_lines = out_file_lines_nam, undu_paras = undu_paras )
		# write new input file
		with open( undu_folder + 'stage/undumag.clc', 'w') as o_f:
			for ind, line in enumerate(out_file_lines_clc) :
				o_f.write(line)
		with open( undu_folder + 'stage/undumag.nam', 'w') as o_f:
			for ind, line in enumerate(out_file_lines_nam) :
				o_f.write(line)
		#copy to undu folder and run
		script_folder = os.getcwd()
		os.chdir( undu_folder + 'stage/' )
		os.system("." + "/../bin/undumag.exe")
		os.chdir( script_folder )
		if not ( copy is None ) : 
			if not os.path.exists(res_folder):
				os.makedirs(res_folder)
			# self.save_para(file = 'para_'+copy+'.txt')
			for resFile in undu_res_files_save :
				if len(copy) > 0 :
					resFile_n = resFile.split('.')[0]
					resFile_e = resFile.split('.')[-1]
					resFile_tmp = resFile_n + '_' + copy + '.' + resFile_e
				else :
					resFile_tmp = resFile
				# copy res files to folder
				os.system( 'cp ' + undu_folder + 'stage/' + resFile + ' ' + res_folder + resFile_tmp )

		beffFile = 'undumag.beff'
		res_beff = self.load_beff_undumag_file(file = undu_folder + 'stage/'+ beffFile)
		res_onax = self.load_on_axis_undumag_file(file = undu_folder + 'stage/undumag_on-axis.dat')
		if not ( self.para['iforce'] == 0 ) :
			force_components, torque_components = self.load_force_undumag_file(file = undu_folder + 'stage/undumag.frc')
			return { 'beff_data' : res_beff, 'on_axis' : res_onax, 'force' : force_components, 'torque' : torque_components }
		return { 'beff_data' : res_beff, 'on_axis' : res_onax }

	def load_on_axis_undumag_file(self, file) : 
		data = pd.read_csv( file, dtype=object, delim_whitespace=True)
		data.columns = ['x','By','Bz','intBy','intBz','int2By','int2Bz','quark']
		cols = data.columns
		for col in cols:
			data[col] = data[col].astype(float)
		return data

	def load_force_undumag_file(self, file) : 
		with open(file, 'r') as o_f:
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

	def load_beff_undumag_file(self, file) : 
		# loop all resulting files we want to save
		out_file_lines_tmp = []
		with open(file, 'r') as o_f:
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
			res.update( { 'ByInt1' : vals[0], 'ByInt2' : vals[2] } )
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
		return res


class point_coords : 

	def __init__(self,x=0.0,y=0.0,z=0.0) : 
		self._x = x
		self._y = y
		self._z = z

	def __sub__(self,pnt) :
		point = point_coords()
		point._x = self._x-pnt._x
		point._y = self._y-pnt._y
		point._z = self._z-pnt._z
		return point

	def __add__(self,pnt) :
		point = point_coords()
		point._x = self._x+pnt._x
		point._y = self._y+pnt._y
		point._z = self._z+pnt._z
		return point

def rotate(pnts,degrees,axis=point_coords(0,0,0),plane='yz'):
	pnts_tmp = copy.deepcopy(pnts)
	degrees_rad = 2*math.pi*degrees/360
	new_pnts = []
	for pnt in pnts_tmp :
		yres = pnt._y*math.cos(degrees_rad)-pnt._z*math.sin(degrees_rad)
		zres = pnt._y*math.sin(degrees_rad)+pnt._z*math.cos(degrees_rad)
		pnt._y = yres
		pnt._z = zres
		if not ( (axis._y == 0) and (axis._z == 0)) :
			yres = axis._y*(1+math.cos(degrees_rad))-axis._z*math.sin(degrees_rad)
			zres = axis._y*math.sin(degrees_rad)+axis._z*(1+math.cos(degrees_rad))
			pnt._y = yres - pnt._y
			pnt._z = zres - pnt._z

	return pnts_tmp

class undu_magnets : 
	"""
	Collects a group of magnets into a list and implements some functionality on that list, like move,mirror,clc-creation,
	extent determination
	magnet_blocks - list of undu_magnet_block_coords objects
	"""
	def __init__(self,magnet_blocks) : 
		self._magnet_blocks = magnet_blocks

	# def create_edge_points(self) : 
	# 	for mag in self._magnet_blocks :
	# 		# mag.create_edge_points()

	def move_it(self,vec) : 
		for mag in self._magnet_blocks :
			mag.move_it(vec=vec)

	def find_magn_name(self,names,fnd_list=None):
		if fnd_list is None :
			fnd_list = []
		for mag in self._magnet_blocks :
			mag.find_magn_name(names=names,fnd_list=fnd_list)
		return fnd_list

	def find_all_names(self,names=None) :
		if names is None :
			names = []
		for mag in self._magnet_blocks :
			mag.find_all_names(names=names)
		return names

	def find_all_mag_blocks(self, mag_blocks = []) : 
		for mag in self._magnet_blocks :
			mag_blocks = mag.find_all_mag_blocks(mag_blocks=mag_blocks)
		return mag_blocks

	def rotate(self,degrees,axis,plane='yz'):
		for mag in self._magnet_blocks :
			mag.rotate(degrees=degrees,axis=axis,plane=plane)

	def mirror(self,coord='x') :
		for mag in self._magnet_blocks :
			mag.mirror(coord=coord)

	def change_segm(self,segm_x,segm_y,segm_z,frac_y=None,frac_z=None) :
		for mag in self._magnet_blocks :
			mag.change_segm(segm_x,segm_y,segm_z,frac_y,frac_z)

	def set_inactive(self) :
		for mag in self._magnet_blocks :
			mag.set_inactive()

	def add_to_clc(self,clc_txt, magns_ignore = None) : 
		for mag in self._magnet_blocks :
			clc_txt = mag.add_to_clc(clc_txt=clc_txt, magns_ignore = magns_ignore)
		return clc_txt

	def get_max_extent(self,maxs=None,mins=None) : 
		for mag in self._magnet_blocks : 
			# print(f"Maxs {maxs} and Mins {mins}")
			p_center, maxs, mins = mag.get_max_extent(maxs=maxs,mins=mins)
		return p_center, maxs, mins

	def set_magnetization(self,magnetization,magn_unit_vec) : 
		for mag in self._magnet_blocks : 
			mag.set_magnetization(magnetization=magnetization,magn_unit_vec=magn_unit_vec)

	def get_center(self) : 
		centers = []
		for mag in self._magnet_blocks :
			if isinstance(mag,undu_magnets) :
				centers.append(mag.get_center())
			else: 
				centers.append(mag._p_center)
		center_res = centers[0]
		for center in centers[1:]:
			center_res._x = 0.5*(center_res._x+center._x)
			center_res._y = 0.5*(center_res._y+center._y)
			center_res._z = 0.5*(center_res._z+center._z)
		return center_res

class undu_magnet_block_coords :
	"""
	Implements basic undumag magnet block
	can be moved, mirrored, incorporated into undumag-clc and the extent can be calculated
	"""

	def __init__(
				self,p_center,pnts=None,len_x=None,len_y=None,len_z=None, magnetization=None, 
				magn_unit_vec=None,	name='name',mother='mother',
				segm_x = 1, segm_y = 1, segm_z = 1, frac_y=1,frac_z=1,
				material = "mag", chamf = None 
				) : 
		"""
		material - "magnet" or "pole"
		chamf - if some float - chamfer is added
		"""
		self._p_center = p_center
		self._len_x = len_x
		self._len_y = len_y
		self._len_z = len_z
		self._magnetization = magnetization
		self._magn_unit_vec = magn_unit_vec
		self._segm_x = segm_x
		self._segm_y = segm_y
		self._segm_z = segm_z
		self._material = material
		self._chamf = chamf
		self._name = name
		self._mother = mother
		self._frac_y = frac_y
		self._frac_z = frac_z
		self._pnts = pnts
		self._inactive = False
		if self._pnts is None :
			self.create_edge_points()

	def find_all_names(self,names=None) :
		if names is None: 
			names = []
		names.append(self._name)
		return names

	def find_all_mag_blocks(self, mag_blocks = []) : 
		mag_blocks.append(self)
		return mag_blocks

	def set_inactive(self) :
		self._inactive = True

	def change_segm(self,segm_x,segm_y,segm_z,frac_y=None,frac_z=None) :
		self._segm_x=segm_x
		self._segm_y=segm_y
		self._segm_z=segm_z
		if not (frac_y is None) :
			if self._frac_y < 1 :
				if frac_y > 1 :
					self._frac_y = 1/frac_y
				else :
					self._frac_y = frac_y
			else :
				if frac_y < 1 :
					self._frac_y = 1/frac_y
				else :
					self._frac_y = frac_y
		if not (frac_z is None) :
			if self._frac_z < 1 :
				if frac_z > 1 :
					self._frac_z = 1/frac_z
				else :
					self._frac_z = frac_z
			else :
				if frac_z < 1 :
					self._frac_z = 1/frac_z
				else :
					self._frac_z = frac_z


	def find_magn_name(self,names,fnd_list=None) :
		if fnd_list is None :
			fnd_list=[]
		not_fnd = False
		for name in names :
			if len(name) > 0 :
				if not (self._name.find(name) >= 0) :
					not_fnd=True
					break
		if not not_fnd :
			fnd_list.append(self)

	def set_magnetization(self,magnetization,magn_unit_vec) : 
		self._magnetization = magnetization
		self._magn_unit_vec = magn_unit_vec

	def create_edge_points(self) : 
		p1 = point_coords(
			x=-self._len_x/2.0,
			y=-self._len_y/2.0,
			z=-self._len_z/2.0,
			)

		p2 = point_coords(
			x=-self._len_x/2.0,
			y=+self._len_y/2.0,
			z=-self._len_z/2.0,
			)

		p3 = point_coords(
			x=-self._len_x/2.0,
			y=-self._len_y/2.0,
			z=+self._len_z/2.0,
			)

		p4 = point_coords(
			x=-self._len_x/2.0,
			y=+self._len_y/2.0,
			z=+self._len_z/2.0,
			)

		p5 = point_coords(
			x=self._len_x/2.0,
			y=-self._len_y/2.0,
			z=-self._len_z/2.0,
			)

		p6 = point_coords(
			x=self._len_x/2.0,
			y=+self._len_y/2.0,
			z=-self._len_z/2.0,
			)

		p7 = point_coords(
			x=self._len_x/2.0,
			y=-self._len_y/2.0,
			z=+self._len_z/2.0,
			)

		p8 = point_coords(
			x=self._len_x/2.0,
			y=+self._len_y/2.0,
			z=+self._len_z/2.0,
			)

		self._pnts = [p1,p2,p3,p4,p5,p6,p7,p8]

	def move_it(self,vec) : 
		p_center = self._p_center
		p_center._x = p_center._x + vec._x
		p_center._y = p_center._y + vec._y
		p_center._z = p_center._z + vec._z
		self._p_center = p_center

	def rotate(self,degrees,axis,plane='yz'):
		self._pnts = rotate(pnts=self._pnts,degrees=degrees)
		new_center = rotate(pnts=[self._p_center],degrees=degrees,axis=axis,plane=plane)[-1]
		self._p_center = new_center
		# need to update the lengths here in order to gete right dimensions later on
		degrees_rad = 2*math.pi*degrees/360
		if ( not (self._len_y is None)) and ( not (self._len_z is None)) :
			new_ly_z = abs(math.sin(degrees_rad)*self._len_y)
			new_lz_z = abs(math.cos(degrees_rad)*self._len_z)
			new_lz = max(new_ly_z,new_lz_z)

			new_ly_y = abs(math.cos(degrees_rad)*self._len_y)
			new_lz_y = abs(math.sin(degrees_rad)*self._len_z)
			new_ly = max(new_ly_y,new_lz_y)
			self._len_y = new_ly
			self._len_z = new_lz

	def mirror(self,coord='x') :
		if coord == 'x' :  
			self._p_center._x = - self._p_center._x
			for ind, pnt in enumerate( self._pnts) : 
				self._pnts[ind]._x = -pnt._x
		elif coord == 'y' :
			self._frac_y = 1/self._frac_y
			self._p_center._y = - self._p_center._y
			for ind, pnt in enumerate( self._pnts) : 
				self._pnts[ind]._y = -pnt._y
		elif coord == 'z' :
			self._frac_z = 1/self._frac_z
			self._p_center._z = - self._p_center._z
			for ind, pnt in enumerate( self._pnts) : 
				self._pnts[ind]._z = -pnt._z

	def add_to_clc(self,clc_txt, magns_ignore = None) : 
		if self._inactive :
			return clc_txt
		if not (magns_ignore is None):
			for ignore in magns_ignore:
				if self._mother==ignore :
					return clc_txt
		my_txt = self.create_clc_txt()
		ind_end = -1
		for ind, line in enumerate(clc_txt) :
			if line.find('*PER END') >= 0 :
				ind_end = ind-1
				break
		if ind_end < 0 : 
			return
		clc_txt[ind_end:ind_end] = my_txt
		return clc_txt

	def create_clc_txt(self) : 
		chamfer = False
		block_type = 'BlockChamf'
		if self._chamf is None :
			corners = True
			block_type = 'Corners'
		else :
			corners = False
		if self._material == 'mag' :
			material = 'Magnet'
			color = 'ColorMag_Hybrid'
		elif self._material == 'pol' :
			material = 'Pole'
			color = 'ColorPol_Hybrid'
		txt = []

		txt.append(f'& {material}\n')
		txt.append(f'{block_type} {self._name} {self._mother} ${color}  !key, name, mother, color\n')
		# txt.append(f'0.0 0.0 0.0\n')
		txt.append(f'{self._p_center._x} {self._p_center._y} {self._p_center._z}\n')
		if material == 'Magnet' :
			txt.append(f'{self._magnetization} {self._magn_unit_vec._x} {self._magn_unit_vec._y} {self._magn_unit_vec._z} $RECIndex_Hybrid  !length bc and comp. of magnetization, material index\n')
		elif material == 'Pole' :
			txt.append(f'$IronIndex_Hybrid  !material index\n')
		if not corners :
			txt.append(f'{self._len_x} {self._len_y} {self._len_z} {self._chamf:.2f}  !dimensions \n')
		txt.append(f'{self._segm_x} {self._segm_y} {self._segm_z} {self._frac_y} {self._frac_z}  !segmentation\n')
		if corners :
			txt.append(f'{len(self._pnts)}\n')
			for pnt in self._pnts : 
				txt.append(f'{pnt._x} {pnt._y} {pnt._z}\n')
		txt.append(' \n')
		return txt

	def get_max_extent(self,maxs=None,mins=None) : 
		"""
		determines the extent of this magnet block and compares to maxs and mins vals given, returns
		the max and min vals
		"""
		max_x = self._pnts[0]._x
		min_x = self._pnts[0]._x
		max_y = self._pnts[0]._y
		min_y = self._pnts[0]._y
		max_z = self._pnts[0]._z
		min_z = self._pnts[0]._z
		for pnt in self._pnts[1:] :
			if pnt._x > max_x :
				max_x = pnt._x
			elif pnt._x < min_x :
				min_x = pnt._x
			if pnt._y > max_y :
				max_y = pnt._y
			elif pnt._y < min_y :
				min_y = pnt._y
			if pnt._z > max_z :
				max_z = pnt._z
			elif pnt._z < min_z :
				min_z = pnt._z

		max_x = self._p_center._x + max_x
		min_x = self._p_center._x + min_x
		max_y = self._p_center._y + max_y
		min_y = self._p_center._y + min_y
		max_z = self._p_center._z + max_z
		min_z = self._p_center._z + min_z
		# max_x = self._p_center._x + self._len_x/2.0
		# min_x = self._p_center._x - self._len_x/2.0
		# max_y = self._p_center._y + self._len_y/2.0
		# min_y = self._p_center._y - self._len_y/2.0
		# max_z = self._p_center._z + self._len_z/2.0
		# min_z = self._p_center._z - self._len_z/2.0
		if maxs is None : 
			maxs = [max_x, max_y, max_z]
		else:
			if max_x > maxs[0] : 
				maxs[0] = max_x
			if max_y > maxs[1] : 
				maxs[1] = max_y
			if max_z > maxs[2] : 
				maxs[2] = max_z
		if mins is None :
			mins = [min_x, min_y, min_z]
		else:
			if min_x < mins[0] : 
				mins[0] = min_x
			if min_y < mins[1] : 
				mins[1] = min_y
			if min_z < mins[2] : 
				mins[2] = min_z
		p_center = point_coords(x=0.5*(maxs[0]+mins[0]),
			y=0.5*(maxs[1]+mins[1]),
			z=0.5*(maxs[2]+mins[2])
			)
		# print(f"Mag fun Maxs {maxs} and Mins {mins}")
		return p_center, maxs, mins

