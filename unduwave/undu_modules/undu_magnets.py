"""
Undu_magnet definitions allowing to construct a 3D magnetic system comprissed of magnet blocks (w/o chamfers)
"""

from unduwave.unduwave_incl import *
from unduwave.undu_modules.undu_parameters import *
import unduwave.helpers.file_folder_helpers as f_h
from unduwave.quantities.quantities import *

class magnetParameters : 
	def __init__( self,
				len_x_main, 
				len_y_main, 
				len_z_main, 
				magnetization=1.0,
				magn_unit_vec='+x',
				segm_x=1,
				segm_y=1,
				segm_z=1,
				material="mag", #"mag","pol", 'NiCuFoil'
				chamf=None,
				frac_y=1,
				frac_z=1,
				) :

		self._len_x_main = len_x_main
		self._len_y_main = len_y_main
		self._len_z_main = len_z_main
		self._magnetization=magnetization
		if isinstance(magn_unit_vec,str) :
			magn_unit_vec=create_magnetization_unit_vec(magn_string=magn_unit_vec)
		self._magn_unit_vec=magn_unit_vec
		self._segm_x = segm_x
		self._segm_y = segm_y
		self._segm_z = segm_z
		self._material=material
		self._chamf=chamf
		self._frac_y=frac_y
		self._frac_z=frac_z

	def __eq__(self,other) :
		if self._len_x_main==other._len_x_main :
			if self._len_y_main==other._len_y_main :
				if self._len_z_main==other._len_z_main :
					if self._magnetization==other._magnetization :
						if self._magn_unit_vec==other._magn_unit_vec :
							if self._material==other._material :
								return True
		return False

class point_coords : 

	def __init__(self,x=0.0,y=0.0,z=0.0) : 
		"""
		Initializes simple 3d vector class

		:param float x: x-component
		:param float y: y-component
		:param float z: z-component
		"""
		self._x = x
		self._y = y
		self._z = z

	def __sub__(self,pnt) :
		"""
		Substraction of vector pnt

		:param point_coords pnt: vector to substract
		"""
		point = point_coords()
		point._x = self._x-pnt._x
		point._y = self._y-pnt._y
		point._z = self._z-pnt._z
		return point

	def __add__(self,pnt) :
		"""
		Addition of vector pnt

		:param point_coords pnt: vector to add
		"""
		point = point_coords()
		point._x = self._x+pnt._x
		point._y = self._y+pnt._y
		point._z = self._z+pnt._z
		return point

	def __eq__(self,other) :
		if self._x == other._x :
			if self._y == other._y :
				if self._z == other._z :
					return True
		return False

def rotate(pnts,degrees,axis=point_coords(0,0,0),plane='yz'):
	"""
	Takes a list of vectors and rotates them (their component in the yz-plane).
	The x-component of the vectors remains unchanged. Returns the rotated vectors

	:param list pnts: List of point_coords to be rotated
	:param float degrees: Degrees by which to rotate [°]
	:param 
	"""
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
	def __init__(self,magnet_blocks,api=None) : 
		self._magnet_blocks = magnet_blocks
		self._api=api
		self.set_api()

	def set_api(self,api=None) :
		if api is None :
			api=self._api 
		for mag in self._magnet_blocks :
			mag.set_api(api=api)

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

	def add_to_clc(self, magns_ignore = None, api=None) : 
		if not api is None :
			self.set_api(api=api)
		for mag in self._magnet_blocks :
			mag.add_to_clc(magns_ignore = magns_ignore)

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

	def get_period_length(self) : 
		"""
		n_magn / n_magn_per_per = n_per
		"""
		names=self.find_all_names()
		posMagn=[]
		for name in names :
			magn=self.find_magn_name(names=[name])[-1]
			posMagn.append({'magn':magn,'pos': magn._p_center._x })
		posMagn=pd.DataFrame(posMagn).sort_values(by=['pos']).to_dict('records')
		indStrt=int(len(posMagn)/2)
		while posMagn[indStrt]['magn']._magnParas._material=='pol' :
			indStrt=indStrt+1
		magnParas0=posMagn[indStrt]['magn']._magnParas
		pos0=posMagn[indStrt]['pos']
		fnd=False
		indEnd=indStrt
		while True :
			indEnd=indEnd+1
			if indStrt >= len(posMagn) :
				break
			magnParas=posMagn[indEnd]['magn']._magnParas
			if (posMagn[indEnd]['pos'] == pos0) :
				continue
			if not (magnParas == magnParas0) :
				continue
			fnd=True
			break
		period=-1
		if fnd :
			period=posMagn[indEnd]['magn']._p_center._x-posMagn[indStrt]['magn']._p_center._x
		pdb.set_trace()
		return period

class undu_magnet_block_coords :
	"""
	Implements basic undumag magnet block
	can be moved, mirrored, incorporated into undumag-clc and the extent can be calculated
	"""

	def __init__(self,
				p_center,
				magnetParameters,
				pnts=None,
				name='name',
				mother='mother',
				api=None
				) : 
		"""
		material - "magnet" or "pole"
		chamf - if some float - chamfer is added
		"""
		self._p_center = copy.deepcopy(p_center)
		self._magnParas=copy.deepcopy(magnetParameters)
		self._name = name
		self._mother = mother
		self._pnts = pnts
		self.set_api(api=api)
		self._inactive = False
		if self._pnts is None :
			self.create_edge_points()

	def set_api(self,api) :
		if not (api is None):
			self._api=api

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
		self._magnParas._segm_x=segm_x
		self._magnParas._segm_y=segm_y
		self._magnParas._segm_z=segm_z
		if not (frac_y is None) :
			if self._magnParas._frac_y < 1 :
				if frac_y > 1 :
					self._magnParas._frac_y = 1/frac_y
				else :
					self._magnParas._frac_y = frac_y
			else :
				if frac_y < 1 :
					self._magnParas._frac_y = 1/frac_y
				else :
					self._magnParas._frac_y = frac_y
		if not (frac_z is None) :
			if self._magnParas._frac_z < 1 :
				if frac_z > 1 :
					self._magnParas._frac_z = 1/frac_z
				else :
					self._magnParas._frac_z = frac_z
			else :
				if frac_z < 1 :
					self._magnParas._frac_z = 1/frac_z
				else :
					self._magnParas._frac_z = frac_z


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
		self._magnParas._magnetization = magnetization
		self._magnParas._magn_unit_vec = magn_unit_vec

	def create_edge_points(self) : 
		p1 = point_coords(
			x=-self._magnParas._len_x_main/2.0,
			y=-self._magnParas._len_y_main/2.0,
			z=-self._magnParas._len_z_main/2.0,
			)

		p2 = point_coords(
			x=-self._magnParas._len_x_main/2.0,
			y=+self._magnParas._len_y_main/2.0,
			z=-self._magnParas._len_z_main/2.0,
			)

		p3 = point_coords(
			x=-self._magnParas._len_x_main/2.0,
			y=-self._magnParas._len_y_main/2.0,
			z=+self._magnParas._len_z_main/2.0,
			)

		p4 = point_coords(
			x=-self._magnParas._len_x_main/2.0,
			y=+self._magnParas._len_y_main/2.0,
			z=+self._magnParas._len_z_main/2.0,
			)

		p5 = point_coords(
			x=self._magnParas._len_x_main/2.0,
			y=-self._magnParas._len_y_main/2.0,
			z=-self._magnParas._len_z_main/2.0,
			)

		p6 = point_coords(
			x=self._magnParas._len_x_main/2.0,
			y=+self._magnParas._len_y_main/2.0,
			z=-self._magnParas._len_z_main/2.0,
			)

		p7 = point_coords(
			x=self._magnParas._len_x_main/2.0,
			y=-self._magnParas._len_y_main/2.0,
			z=+self._magnParas._len_z_main/2.0,
			)

		p8 = point_coords(
			x=self._magnParas._len_x_main/2.0,
			y=+self._magnParas._len_y_main/2.0,
			z=+self._magnParas._len_z_main/2.0,
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
		if ( not (self._magnParas._len_y_main is None)) and ( not (self._magnParas._len_z_main is None)) :
			new_ly_z = abs(math.sin(degrees_rad)*self._magnParas._len_y_main)
			new_lz_z = abs(math.cos(degrees_rad)*self._magnParas._len_z_main)
			new_lz = max(new_ly_z,new_lz_z)

			new_ly_y = abs(math.cos(degrees_rad)*self._magnParas._len_y_main)
			new_lz_y = abs(math.sin(degrees_rad)*self._magnParas._len_z_main)
			new_ly = max(new_ly_y,new_lz_y)
			self._magnParas._len_y_main = new_ly
			self._magnParas._len_z_main = new_lz

	def mirror(self,coord='x') :
		if coord == 'x' :  
			self._p_center._x = - self._p_center._x
			for ind, pnt in enumerate( self._pnts) : 
				self._pnts[ind]._x = -pnt._x
		elif coord == 'y' :
			self._magnParas._frac_y = 1/self._magnParas._frac_y
			self._p_center._y = - self._p_center._y
			for ind, pnt in enumerate( self._pnts) : 
				self._pnts[ind]._y = -pnt._y
		elif coord == 'z' :
			self._magnParas._frac_z = 1/self._magnParas._frac_z
			self._p_center._z = - self._p_center._z
			for ind, pnt in enumerate( self._pnts) : 
				self._pnts[ind]._z = -pnt._z

	def add_to_clc(self, magns_ignore = None,api=None) : 
		if api is None :
			api=self._api
		if self._inactive :
			return 
		if not (magns_ignore is None):
			for ignore in magns_ignore:
				if self._mother==ignore :
					return 
		my_txt = self.create_clc_txt()
		clc_txt=api._prog_paras.in_file_clc_lines.get()
		if len(clc_txt) < 1 :
			api.load_clc_raw()
			clc_txt=api._prog_paras.in_file_clc_lines.get()

		ind_end = -1
		for ind, line in enumerate(clc_txt) :
			if line.find('*PER END') >= 0 :
				ind_end = ind-1
				break
		if ind_end < 0 : 
			return
		clc_txt[ind_end:ind_end] = my_txt
		api._prog_paras.in_file_clc_lines.set(clc_txt)

	def create_clc_txt(self) : 
		chamfer = False
		block_type = 'BlockChamf'
		if self._magnParas._chamf is None :
			corners = True
			block_type = 'Corners'
		else :
			corners = False
		if self._magnParas._material == 'mag' :
			material = 'Magnet'
			materialInd=1
			color = 'ColorMag_Hybrid'
		elif self._magnParas._material == 'pol' :
			material = 'Pole'
			materialInd=2
			color = 'ColorPol_Hybrid'
		elif self._magnParas._material == 'NiCuFoil' :
			material = 'Pole'
			materialInd=3
			color = 'ColorFoil'
		txt = []

		txt.append(f'& {material}\n')
		txt.append(f'{block_type} {self._name} {self._mother} ${color}  !key, name, mother, color\n')
		# txt.append(f'0.0 0.0 0.0\n')
		txt.append(f'{self._p_center._x} {self._p_center._y} {self._p_center._z}\n')
		if material == 'Magnet' :
			txt.append(f'{self._magnParas._magnetization} {self._magnParas._magn_unit_vec._x} {self._magnParas._magn_unit_vec._y} {self._magnParas._magn_unit_vec._z} {materialInd}  !length bc and comp. of magnetization, material index\n')
		elif material == 'Pole' :
			txt.append(f'{materialInd}  !material index\n')
		elif material == 'NiCuFoil' :
			txt.append(f'{materialInd}  !material index\n')
		if not corners :
			txt.append(f'{self._magnParas._len_x_main} {self._magnParas._len_y_main} {self._magnParas._len_z_main} {self._magnParas._chamf:.2f}  !dimensions \n')
		txt.append(f'{self._magnParas._segm_x} {self._magnParas._segm_y} {self._magnParas._segm_z} {self._magnParas._frac_y} {self._magnParas._frac_z}  !segmentation\n')
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

def get_magn_ignore(n_magn_ignore,ignore_level,magns_ignore=[]):
	"""
	n_magn_ignore : number of magnet to be ignored
	ignore_level : 'u'/'l' - upper or lower rows
	"""
	ignore_rowl = ignore_level+'l'
	ignore_rowr = ignore_level+'r'
	magns_ignore_new = [ f'{ignore_rowl}_row_mag{n_magn_ignore}',f'{ignore_rowr}_row_mag{n_magn_ignore}' ]
	if len(magns_ignore) > 0 :
		magns_ignore[0:0]=magns_ignore_new
	else :
		magns_ignore = magns_ignore_new
	return magns_ignore

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

def get_magnet_period(magn_collection, denominator, nperiod, n_el_per_period = 4, offset=0 ) :
	first_magn = offset+nperiod*n_el_per_period
	magn_nums = [ i for i in range(first_magn,first_magn+4) ]
	period = get_magnet_struct(
		magn_collection=magn_collection, 
		denominator=denominator, 
		magn_nums = magn_nums
		)
	return period

def get_magnet_struct(magn_collection, denominator, magn_nums = []) :
	"""
	Returns list of all objects whose name contains one of: 
	{ denominator and str(magn_nums[i]) | forall i}
	""" 
	fnd_list=[]
	if len(magn_nums) > 0 :
		for magn_num in magn_nums :
			if denominator is None :
				denos = [f'obj_{magn_num}' ]
			else:
				denos = [denominator,f'obj_{magn_num}' ]
			magn_collection.find_magn_name(names=denos,fnd_list=fnd_list)
	else :
		denos = [denominator]
		magn_collection.find_magn_name(names=denos,fnd_list=fnd_list)
	magnets_obj = undu_magnets(magnet_blocks=fnd_list)
	return magnets_obj
