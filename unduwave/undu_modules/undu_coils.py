"""
Defining the functionality to representations of current-carrying coils.
"""

from unduwave.unduwave_incl import *
from unduwave.undu_modules.undu_parameters import *
import unduwave.helpers.file_folder_helpers as f_h
from unduwave.quantities.quantities import *
import unduwave.undu_modules.undu_blocks as undu_blocks
from unduwave.attribute_classes.attributes import _attribute
from unduwave.attribute_classes.attributes import _attribute_collection

class coil(undu_blocks.undumagBlockObject) :

	def __init__(self, 
			center, 
			coil_type, 
			current, 
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
			name='',
			filling = None, 
			n_windings = None,
			api=None,
			parentName='',
			**kwargs,
		) :
		super(coil, self).__init__(
			center=center,
			magnParas=None,
			pnts=None,
			name=name,
			parentName=parentName,
			api=api,
			**kwargs,
			)

		"""
		Coil class representing a currect carrying coil.
		:param str coil_type: The type of coil. Can be one of the following: \n
			Rectangular: A coil with 90° edges. \n
			RectWindings: A coil with rounded edges. \n
		:param float current: The current in the coil [mA]
		:param undu_magnets.point_coords center: The 3 center coordinates. [mm]
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
		self.coil_type = _attribute(coil_type)
		self.I = _attribute(current)
		self.normal_vec = _attribute(normal_vec)
		self.rot_angle = _attribute(rot_angle)
		self.length = _attribute(length)
		self.inner_z = _attribute(inner_z)
		self.outer_z = _attribute(outer_z)
		self.inner_radius = _attribute(inner_radius)
		self.height = _attribute(height)
		self.n_vert = _attribute(n_vert)
		self.n_hor = _attribute(n_hor)
		self.n_rad = _attribute(n_rad)
		self.filling = _attribute(filling)
		self.n_windings = _attribute(n_windings)

	def getDictInfo(self) : 
		if self._createdWithPoints : 
			pntsD=[ pnts.tolist() for pnts in self._pnts]
		else :
			pntsD=None

		myDict=super().getDictInfo()
		myDict.update({
			'type' : 'coil',
			'coil_type' : self.coil_type(),
			'I':self.I(), 
			'normal_vec':self.normal_vec(), 
			'rot_angle':self.rot_angle(), 
			'length':self.length(),
			'inner_z':self.inner_z(),
			'outer_z':self.outer_z(),
			'inner_radius':inner_radius(),
			'height':self.height(),
			'n_vert':self.n_vert(), #"magnet","pole", 'NiCuFoil'
			'n_hor':self.n_hor(),
			'n_rad':n_rad(),
			'filling':filling(),
			'n_windings' : self.n_windings(),
		})
		return myDict

	@staticmethod
	def loadFromDict(dictParas) : 

		if not (dictParas['type'] == 'coil') :
			return

		thecoil=coil(
			center=dictParas['center'], 
			coil_type=dictParas['coil_type'], 
			current=dictParas['I'], 
			normal_vec=dictParas['normal_vec'], 
			rot_angle=dictParas['rot_angle'], 
			length=dictParas['length'], 
			inner_z=dictParas['inner_z'], 
			outer_z=dictParas['outer_z'], 
			inner_radius=dictParas['inner_radius'], 
			height=dictParas['height'], 
			n_vert=dictParas['n_vert'], 
			n_hor=dictParas['n_hor'], 
			n_rad=dictParas['n_rad'], 
			name=dictParas['name'], 
			filling=dictParas['filling'], 
			n_windings=dictParas['n_windings'], 
			api=None,
			parentName=dictParas['parentName'], 
		)
		return thecoil

	def add_to_clc(self, magns_ignore = None,api=None) : 

		"""
		Gets the clc-text for this coil and adds it to the api._prog_paras.in_file_clc_lines list.
		The in_file_clc_lines lines object can be initialized or, if empty, will be initialized via the api.
		"""

		# get the current api clc-text
		if (api is None) :
			api=self._api
		else :
			self.set_api(api=api)

		if not (magns_ignore is None):
			for ignore in magns_ignore:
				if self.parentName==ignore :
					return 

		clc_txt=api._prog_paras.in_file_clc_lines.get()
		if len(clc_txt) < 1 :
			# initialize the clc with the api if nothing there
			api.load_clc_raw()
			clc_txt=api._prog_paras.in_file_clc_lines.get()

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
		api._prog_paras.in_file_clc_lines.set(clc_txt)

	def create_clc_txt(self) : 
		"""
		Creates the actual clc text from the object's parameters

		:return list: The list containing str representing the content of the undumag clc file
		"""
		txt = []
		if self.coil_type() == 'Rectangular' :
			txt.append(f'& Coil\n')
			txt.append("* Current, center, normal vector, rotation-angle, total length,\n")
			txt.append("* inner width, outer width, inner radius, height, vert/hori/rad\n")
			txt.append("* divisions, color-index\n")
			txt.append(f"{self.coil_type()}\n")
			txt.append(f"{self.I()} {self._center[0]} {self._center[1]} {self._center[2]} ") 
			txt.append(f"{self.normal_vec()[0]} {self.normal_vec()[1]} {self.normal_vec()[2]} {self.rot_angle()} {self.length()} ")
			txt.append(f"{self.inner_z()} {self.outer_z()} {self.inner_radius()} {self.height()} {self.n_vert()} ")
			txt.append(f"{self.n_hor()} {self.n_rad()} {2}")
			txt.append(' \n\n')
		if self.coil_type() == 'RectWindings' :
			txt.append(f'& Coil\n')
			txt.append("* $Current $Filling $nWinding $xCoil $ymCoil $zCoil $Vnx $Vny $Vnz $Ang\n")
			txt.append("* $xLenOut $zLenIn $zLenOut $RectRi $Height\n")
			txt.append("* $nDivHeight $nDivWind $nDivArc $nColor\n")
			txt.append(f"{self.coil_type()}\n")
			txt.append(f"{self.I()} {self.filling()} {self.n_windings()} {self._center[0]} {self._center[1]} {self._center[2]} ") 
			txt.append(f"{self.normal_vec()[0]} {self.normal_vec()[1]} {self.normal_vec()[2]} {self.rot_angle()} {self.length()} ")
			txt.append(f"{self.inner_z()} {self.outer_z()} {self.inner_radius()} {self.height()} {self.n_vert()} ")
			txt.append(f"{self.n_hor()} {self.n_rad()} {2}")
			txt.append(' \n\n')

		return txt

	def get_max_extent(self,maxs=None,mins=None) : 
		# v1=np.cross(self.normal_vec,np.array([1.0,0.0,0.0]))
		# if v1 == np.array([0.0,0.0,0.0]) :
		# 	v1=np.array([0.0,1.0,0.0])
		# 	v2=np.array([0.0,0.0,1.0])
		# else :
		# 	v2=np.cross(v1,self.normal_vec)
		# 	v1=v1/np.linalg.norm(v1)
		# 	v2=v2/np.linalg.norm(v2)
		return self._center(), maxs, mins

	def get_copy(self,name=None,parentName=None) :
		return self