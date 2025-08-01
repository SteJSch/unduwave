"""
Defining the functionality to representations of current-carrying coils.
"""

from unduwave.unduwave_incl import *
from unduwave.undu_modules.undu_parameters import *
import unduwave.helpers.file_folder_helpers as f_h
from unduwave.quantities.quantities import *

class undu_coils : 

	def __init__(self,coils,api=None) : 
		"""
		A container class containing a list of coil objects which are to be created and simulated with Undumag.
		:param list coils: A list of coil objects.
		:param api: The undu_api object carrying out the simulations.
		:type api: undu_api or None
		"""
		self._coils = coils
		self._api=api

	def move_it(self,vec) : 
		"""
		Moving the center of the coil by vector vec.

		:param list vec: A list containing the (x,y,z) components of the translation vector vec.
		"""
		for coil in self._coils :
			coil.move_it(vec=vec)

	def add_to_clc(self) : 
		"""
		Add all daughter coils to the Undumag clc file for simulation.
		"""
		for coil in self._coils :
			coil.add_to_clc()

	def set_all_apis(self) : 
		"""
		Sets all daughter coil-api objects to the api object of this collection.
		"""

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
		"""
		Coil class representing a currect carrying coil.
		:param str coil_type: The type of coil. Can be one of the following: \n
			Rectangular: A coil with 90° edges. \n
			RectWindings: A coil with rounded edges. \n
		:param float current: The current in the coil [mA]
		:param undu_magnets.point_coords center_coords: The 3 center coordinates. [mm]
		:param undu_magnets.point_coords normal_vec: Normal vector. The opening in the middle of the coil lies perpendicular to normal_vec. A normal_vec along (0,1,0) is the standard to which the parameter names refer!
		:param float rot_angle: angle by which coil is rotated aroung normal_vec?
		:param float length: Extension in x-direction. [mm]
		:param float inner_z: Extension in z-direction of hole inside coil. [mm]
		:param float outer_z: Extension in z-direction of full coil. The difference of outer and inner z determines the width of the coils running around the central hole. [mm]
		:param float inner_radius: The inner radius of the curvature of the coil edges. [mm]
		:param float height: The extent of the coil in y-direction. [mm]
		:param float n_vert,n_hor,n_rad: Number of segmentations.
		:param filling: The percentage of the cross-section of the coil being filled with current carrying wire
		:type filling: float or None
		:param int n_windings: The number of coils in a cross-section
		:type n_windings: int or None
		:param api: The undu_api object carrying out the simulations.
		:type api: undu_api or None
		"""

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
		"""
		Setting the internal api if it is not None

		:param undu_api api: The api object.
		"""
		if not (api is None):
			self._api=api

	def add_to_clc(self) : 
		"""
		Gets the clc-text for this coil and adds it to the api._prog_paras.in_file_clc_lines list.
		The in_file_clc_lines lines object can be initialized or, if empty, will be initialized via the api.
		"""

		# get the current api clc-text
		clc_txt=self._api._prog_paras.in_file_clc_lines.get()
		if len(clc_txt) < 1 :
			# initialize the clc with the api if nothing there
			self._api.load_clc_raw()
			clc_txt=self._api._prog_paras.in_file_clc_lines.get()

		# get the clc text for this object
		my_txt = self.create_clc_txt()
		ind_end = -1
		for ind, line in enumerate(clc_txt) :
			# find the start of th coils section in clc list
			if line.find('*COILS END') >= 0 :
				ind_end = ind-1
				break
		if ind_end < 0 : 
			return
		# add this clc text to the start of the clc-coil section
		clc_txt[ind_end:ind_end] = my_txt
		# set the new clc text for the api
		self._api._prog_paras.in_file_clc_lines.set(clc_txt)

	def move_it(self,vec) : 
		"""
		Translates the center of the voil by the vector vec. 

		:param undu_magnets.point_coords vec: The translation vector.
		"""
		p_center = self._center_coords
		p_center._x = p_center._x + vec._x
		p_center._y = p_center._y + vec._y
		p_center._z = p_center._z + vec._z
		self._center_coords = p_center

	def create_clc_txt(self) : 
		"""
		Creates the actual clc text from the object's parameters

		:return list: The list containing str representing the content of the undumag clc file
		"""
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
