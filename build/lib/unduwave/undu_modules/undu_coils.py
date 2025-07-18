"""
Contains the wave_from_b class that incorporates the API from python to WAVE and the function create_wave_instance that 
returns an instance of that class
"""

from unduwave.unduwave_incl import *
from unduwave.undu_modules.undu_parameters import *
import unduwave.helpers.file_folder_helpers as f_h
from unduwave.quantities.quantities import *

class undu_coils : 

	def __init__(self,coils,api=None) : 
		self._coils = coils
		self._api=api

	def move_it(self,vec) : 
		for coil in self._coils :
			coil.move_it(vec=vec)

	def add_to_clc(self,clc_txt) : 
		for coil in self._coils :
			coil.add_to_clc()

	def set_all_apis(self) : 
		for coil in self._coils :
			coil.set_api(api=self._api)

class coil :

	def __init__(self, 
			coil_type, 
			current, 
			center_coords, 
			normal_vec, 
			rot_angle, 
			length, 
			inner_z, 
			outer_z, 
			inner_radius, 
			height, 
			n_vert, 
			n_hor, 
			n_rad, 
			filling = None, 
			n_windings = None,
			api=None,
			) : 

		self._type = coil_type
		self._I = current
		self._center_coords = center_coords
		self._normal_vec = normal_vec
		self._rot_angle = rot_angle
		self._length = length
		self._inner_z = inner_z
		self._outer_z = outer_z
		self._inner_radius = inner_radius
		self._height = height
		self._n_vert = n_vert
		self._n_hor = n_hor
		self._n_rad = n_rad
		self._filling = filling
		self._n_windings = n_windings
		self.set_api(api=api)

	def set_api(self,api) :
		if not (api is None):
			self._api=api

	def add_to_clc(self) : 

		my_txt = self.create_clc_txt()
		clc_txt=self._api._prog_paras.in_file_clc_lines.get()
		if len(clc_txt) < 1 :
			self._api.load_clc_raw()
			clc_txt=self._api._prog_paras.in_file_clc_lines.get()

		my_txt = self.create_clc_txt()
		ind_end = -1
		for ind, line in enumerate(clc_txt) :
			if line.find('*COILS END') >= 0 :
				ind_end = ind-1
				break
		if ind_end < 0 : 
			return
		clc_txt[ind_end:ind_end] = my_txt
		self._api._prog_paras.in_file_clc_lines.set(clc_txt)

	def move_it(self,vec) : 
		p_center = self._center_coords
		p_center._x = p_center._x + vec._x
		p_center._y = p_center._y + vec._y
		p_center._z = p_center._z + vec._z
		self._center_coords = p_center

	def create_clc_txt(self) : 
		txt = []
		if self._type == 'Rectangular' :
			txt.append(f'& Coil\n')
			txt.append("* Current, center, normal vector, rotation-angle, total length,\n")
			txt.append("* inner width, outer width, inner radius, height, vert/hori/rad\n")
			txt.append("* divisions, color-index\n")
			txt.append(f"{self._type}\n")
			txt.append(f"{self._I} {self._center_coords._x} {self._center_coords._y} {self._center_coords._z} ") 
			txt.append(f"{self._normal_vec._x} {self._normal_vec._y} {self._normal_vec._z} {self._rot_angle} {self._length} ")
			txt.append(f"{self._inner_z} {self._outer_z} {self._inner_radius} {self._height} {self._n_vert} ")
			txt.append(f"{self._n_hor} {self._n_rad} {2}")
			txt.append(' \n\n')
		if self._type == 'RectWindings' :
			txt.append(f'& Coil\n')
			txt.append("* $Current $Filling $nWinding $xCoil $ymCoil $zCoil $Vnx $Vny $Vnz $Ang\n")
			txt.append("* $xLenOut $zLenIn $zLenOut $RectRi $Height\n")
			txt.append("* $nDivHeight $nDivWind $nDivArc $nColor\n")
			txt.append(f"{self._type}\n")
			txt.append(f"{self._I} {self._filling} {self._n_windings} {self._center_coords._x} {self._center_coords._y} {self._center_coords._z} ") 
			txt.append(f"{self._normal_vec._x} {self._normal_vec._y} {self._normal_vec._z} {self._rot_angle} {self._length} ")
			txt.append(f"{self._inner_z} {self._outer_z} {self._inner_radius} {self._height} {self._n_vert} ")
			txt.append(f"{self._n_hor} {self._n_rad} {2}")
			txt.append(' \n\n')

		return txt
