"""
Undu_magnet definitions allowing to construct a 3D magnetic system comprissed of magnet blocks (w/o chamfers)
"""

from unduwave.unduwave_incl import *
from unduwave.undu_modules.undu_parameters import *
import unduwave.helpers.file_folder_helpers as f_h
from unduwave.quantities.quantities import *
from unduwave.attribute_classes.attributes import _attribute
from unduwave.attribute_classes.attributes import _attribute_collection

class unduBlockError(Exception) :
	pass

class magnetic_material(_attribute_collection) : 

	def __init__( self,
			material_id=None, 
			base_material_type=None, # "magnet" or "pole", "NiCuFoil"
			ksi_easy=None,
			ksi_perp=None,
			magnetization_data=None,
			magnetization_file=None,
			loadedMaterialList=None,
			) :


		self._base_material_type=_attribute(base_material_type)
		self._material_id=_attribute(material_id)
		self._ksi_easy=_attribute(ksi_easy)
		self._ksi_perp=_attribute(ksi_perp)
		self._magnetization_file=_attribute(magnetization_file)
		self._magnetization_data=_attribute(magnetization_data)
		self._ini=False
		self.initialize()
		super().__init__()

	def initialize(self) : 
		try:
			if not (self._magnetization_file() is None) :
				self.load_from_material_file()
				self._ini=True
			elif not (self._magnetization_data() is None) :
				self._ini=True
			elif not (self._ksi_easy() is None) :
				if not (self._ksi_perp() is None) :
					self._ini=True
			# if not self._ini :
			# 	raise unduBlockError("You have to set either a material file, material data or parameters.")
			self._ini=True
		except Exception as e:
			print(e)

	def load_from_material_file(self) :

		with open(self._magnetization_file()) as f:
			lines = f.readlines()
		try:
			firstLine=lines[0]
			matType=firstLine.split('type:')[1].strip()
			self._base_material_type.set(matType)
			scndLine=lines[1]
			matID=scndLine.split('id:')[1].strip()
			self._material_id.set(matID)

			if self._base_material_type() == 'magnet' :
				paras=lines[3].split(' ')
				paras=[float(pa) for pa in paras]
				self._ksi_easy.set(paras[0])
				self._ksi_perp.set(paras[1])
			elif self._base_material_type() == 'pole' :
				lines=lines[3:]
				magnetDat=[]
				for line in lines:
					spl=line.strip().split(' ')
					for el in spl: 
						el=el.strip()
						if len(el) > 0 :
							magnetDat.append(float(el))
				magnetDat=np.array(magnetDat)
				magnetDat=np.reshape(magnetDat,(int(len(magnetDat)/2),2))
				self._magnetization_data.set(magnetDat)
		except Exception as e:
			print(e)

	def write_undumag_material_file(self,folder) :
		if not self._ini :
			return
		fullFile=os.path.join(folder,f'{self._material_id}.dat')
		with open(fullFile, 'w') as f:
			if self._base_material_type() == 'magnet' :
				f.write(f'{(self._ksi_easy()+1.0):.4f} {(self._ksi_perp()):.4f} ! mu_Par and ksi_Per')
			elif self._base_material_type() == 'pole' :
				data=self._magnetization_data()
				for dL in data :
					f.write(f'{dL[0]} {dL[1]}\n')

	@staticmethod
	def get_material_from_list(theList,material_id) : 
		for el in theList :
			if el._material_id()==material_id :
				return el

	@staticmethod
	def add_new_to_material_list(theList,other) : 
		fnd=False
		indList=0
		for ind,el in enumerate(theList) :
			if el == other :
				fnd=True
				indList=ind
				break
		if not fnd :
			theList.append(other)
			indList=len(theList)-1
		return theList, indList

	def __eq__(self,other) : 
		"""
		Returns True if two materials share same definition.
		"""
		dataEqual=False
		if not (self._material_id() == other._material_id()) :
			return False
		if self._base_material_type() is None :
			return False
		if other._base_material_type() is None :
			return False
		if self._base_material_type() == other._base_material_type() :
			if not (self._ksi_easy() is None) :

				if not (self._ksi_perp() is None) :
					if (self._ksi_easy() == other._ksi_easy()) and (self._ksi_perp() == other._ksi_perp()) :
						dataEqual=True
					else:
						return False
				else:
					return False

			elif not (self._magnetization_data() is None) :
				if not (other._magnetization_data() is None) :
					if np.array_equal(self._magnetization_data(),other._magnetization_data()) :
						dataEqual=True
				else:
					return False
		# if dataEqual :
		# 	if not (self._mat_id() == other._mat_id()) :
		# 		self._mat_id.set( other._mat_id() )
		return dataEqual

class magParameters(_attribute_collection) : 
	def __init__( self,
			material_id=None,
			len_x_main=None, 
			len_y_main=None, 
			len_z_main=None, 
			magnetization=1.0,
			magn_unit_vec='+x',
			segm_x=1,
			segm_y=1,
			segm_z=1,
			chamf=None,
			frac_y=1,
			frac_z=1,
			**kwargs,
			) :

		super(magParameters, self).__init__(**kwargs)
		self._len_x_main = _attribute(len_x_main)
		self._len_y_main = _attribute(len_y_main)
		self._len_z_main = _attribute(len_z_main)
		self._magnetization = _attribute(magnetization)
		if isinstance(magn_unit_vec,str) :
			magn_unit_vec=create_magnetization_unit_vec(magn_string=magn_unit_vec)
		self._magn_unit_vec = _attribute(magn_unit_vec)
		self._segm_x = _attribute(segm_x)
		self._segm_y = _attribute(segm_y)
		self._segm_z = _attribute(segm_z)
		self._material_id = _attribute(material_id)
		self._chamf = _attribute(chamf)
		self._frac_y = _attribute(frac_y)
		self._frac_z = _attribute(frac_z)
		super().__init__()

	def __eq__(self,other) :
		if self._len_x_main()==other._len_x_main() :
			if self._len_y_main()==other._len_y_main() :
				if self._len_z_main()==other._len_z_main() :
					if self._magnetization()==other._magnetization() :
						if np.array_equal(self._magn_unit_vec(),other._magn_unit_vec()):
							if self._material_id()==other._material_id() :
								return True
		return False

def rotate(pnts,degrees,axis=np.array([0,0,0]),plane='yz'):
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
		yres = pnt[1]*math.cos(degrees_rad)-pnt[2]*math.sin(degrees_rad)
		zres = pnt[1]*math.sin(degrees_rad)+pnt[2]*math.cos(degrees_rad)
		pnt[1] = yres
		pnt[2] = zres
		if not ( (axis[1] == 0) and (axis[2] == 0)) :
			yres = axis[1]*(1+math.cos(degrees_rad))-axis[2]*math.sin(degrees_rad)
			zres = axis[1]*math.sin(degrees_rad)+axis[2]*(1+math.cos(degrees_rad))
			pnt[1] = yres - pnt[1]
			pnt[2] = zres - pnt[2]

	return pnts_tmp

class undumagObjectList : 
	"""
	Collects a group of magnets into a list and implements some functionality on that list, like move,mirror,clc-creation,
	extent determination
	magnet_blocks - list of undu_magnet_block_coords objects
	"""
	def __init__(self,
			magnet_blocks=[],
			name='',
			parentName='',
			api=None,
			center=np.array([0.0,0.0,0.0]),
			includeInStruct=True,
			**kwargs,
			) : 
		self._magnet_blocks = magnet_blocks
		self._api=api
		self._name=name
		self._parentName=parentName
		self._center=center
		self.set_api()
		# self.move_it(vec=self._center,moveSelf=False)
		if includeInStruct :
			self.set_struct_name(parentName=self._parentName)
			self.set_absolute_position(parentPos=np.array([0.0,0.0,0.0]))

	def getDictInfo(self,onlyBase=False) : 
		myDict={
			'name' : self._name,
			'parentName' : self._parentName,
			'center' : self._center.tolist(),
			'type' : 'objList',
		}
		if onlyBase :
			return myDict
		blocks=[]
		for mag in self._magnet_blocks :
			block= mag.getDictInfo()
			blocks.append(block)
		myDict.update({'blocks':blocks})
		return myDict

	@staticmethod
	def loadFromDict(dictParas) : 
		if not (dictParas['type'] == 'objList') :
			return
		magnet_blocks=[]
		for block in dictParas['blocks'] : 
			if block['type'] == 'objList' : 
				nBlock=undumagObjectList.loadFromDict(dictParas=block)
				magnet_blocks.append(nBlock)
			elif block['type'] == 'obj' :
				nBlock=undumagBlockObject.loadFromDict(dictParas=block)
				magnet_blocks.append(nBlock)
		objList=undumagObjectList(
			magnet_blocks=magnet_blocks,
			name=dictParas['name'],
			parentName=dictParas['parentName'],
			api=None,
			center=np.array(dictParas['center']),
			)
		return objList

	@staticmethod
	def getDictFromYAMLFile(file) :
		try:
			with open(file, 'r') as fileYAML:
				unduDict = yaml.safe_load(fileYAML)
		except: 
			pdb.set_trace()
			return []
		return unduDict

	def writeToParameterFile(self,file) : 
		dictLoad=self.getDictFromYAMLFile(file=file)
		myDict=self.getDictInfo()
		fnd=False
		for indd, dictL in enumerate(dictLoad) :
			if dictL['name'] == myDict['name'] :
				dictLoad[indd] = myDict
				fnd=True
				break
		if not fnd :
			dictLoad.append(myDict)
		with open(file,"w") as openFile :
			yaml.dump(
				dictLoad,
				openFile,
				default_flow_style=False,
				indent=4,
				sort_keys=False,
				)

	def set_struct_name(self,parentName='') :
		if len(parentName) > 0 :
			self._parentName=parentName
			self._nameStruct=parentName+'_' + self._name
		else :
			self._nameStruct=self._name
		try:
			for mag in self._magnet_blocks :
				mag.set_struct_name(parentName=self._nameStruct)
		except:
			pdb.set_trace()

	def set_api(self,api=None) :
		if api is None :
			api=self._api 
		self._api=api
		for mag in self._magnet_blocks :
			mag.set_api(api=api)

	def get_copy(self,name=None,parentName=None) :
		if name is None :
			name=self._name
		if parentName is None :
			parentName=self._parentName
		magnet_blocks=[]
		for mag in self._magnet_blocks :
			magnet_blocks.append( mag.get_copy() )
		newObj = undumagObjectList(
			name=name,
			parentName=parentName,
			magnet_blocks=magnet_blocks,
			api=self._api,
			center=copy.deepcopy(self._center),
			)
		return newObj

	# def create_edge_points(self) : 
	# 	for mag in self._magnet_blocks :
	# 		# mag.create_edge_points()

	def move_it(self,vec) : 
		self._center=self._center+vec
		self.set_absolute_position()

	def set_absolute_position(self,parentPos=np.array([0.0,0.0,0.0])) :
		"""
		absolute position is updated for me and all children
		"""
		self._absCenter=copy.deepcopy(parentPos)+copy.deepcopy(self._center)
		for mag in self._magnet_blocks :
			mag.set_absolute_position(parentPos=self._absCenter)

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

		if coord == 'x' :  
			self._absCenter[0] = - self._absCenter[0]
			self._center[0] = - self._center[0]
		elif coord == 'y' :
			self._absCenter[1] = - self._absCenter[1]
			self._center[1] = - self._center[1]
		elif coord == 'z' :
			self._absCenter[2] = - self._absCenter[2]
			self._center[2] = - self._center[2]
		for mag in self._magnet_blocks :
			mag.mirror(coord=coord)

		self.set_absolute_position()

	def change_segm(self,segm_x,segm_y,segm_z,frac_y=None,frac_z=None) :
		for mag in self._magnet_blocks :
			mag.change_segm(segm_x,segm_y,segm_z,frac_y,frac_z)

	def set_inactive(self) :
		for mag in self._magnet_blocks :
			mag.set_inactive()

	def add_to_clc(self, magns_ignore = None, api=None, refresh=False) : 
		if not (api is None) :
			self.set_api(api=api)
			self._api._prog_paras.load_all_magnetic_material_files_from_folders()
		if refresh :
			if not (self._api is None):
				self._api._prog_paras.in_file_clc_lines.set([])
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
			if isinstance(mag,undumagObjectList) :
				centers.append(mag.get_center())
			else: 
				centers.append(mag._center)
		center_res = centers[0]
		for center in centers[1:]:
			center_res[0] = 0.5*(center_res[0]+center[0])
			center_res[1] = 0.5*(center_res[1]+center[1])
			center_res[2] = 0.5*(center_res[2]+center[2])
		return center_res

	def get_period_length(self,name_hints=['ll']) : 
		"""
		n_magn / n_magn_per_per = n_per
		"""
		if self._api is None :
			return -1
		names=self.find_all_names()
		posMagn=[]
		for name in names :
			take=True
			for hint in name_hints :
				if not (name.find(hint) >= 0) :
					take=False
					break
			if not take:
				continue
			magn=self.find_magn_name(names=[name])[-1]
			posMagn.append({'magn':magn,'pos': magn._absCenter[0] })
		posMagn=pd.DataFrame(posMagn).sort_values(by=['pos']).to_dict('records')

		magMaterialsList=self._api._prog_paras.magnetic_materials()

		indStrt=int(len(posMagn)/2)
		while True :

			material=magnetic_material.get_material_from_list(
				theList=magMaterialsList,
				material_id=posMagn[indStrt]['magn']._magnParas._material_id(),
				)

			if material._base_material_type()=='pole' :
				indStrt=indStrt+1
			else :
				break

		magnParas0=posMagn[indStrt]['magn']._magnParas
		pos0=posMagn[indStrt]['pos']
		fnd=False
		indEnd=indStrt
		while True :
			indEnd=indEnd+1
			if indEnd >= len(posMagn) :
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
			period=posMagn[indEnd]['magn']._absCenter[0]-posMagn[indStrt]['magn']._absCenter[0]
		return period

class undumagBlockObject :
	"""
	Implements basic undumag magnet block
	can be moved, mirrored, incorporated into undumag-clc and the extent can be calculated
	"""

	def __init__(self,
				center,
				magnParas=None,
				pnts=None,
				name='name',
				parentName='',
				api=None
				) : 
		"""
		chamf - if some float - chamfer is added
		"""
		self._center = copy.deepcopy(center)
		self.set_absolute_position()
		self._magnParas=copy.deepcopy(magnParas)
		self._name = name
		self._parentName=parentName
		self._mother='mother'
		self.set_struct_name(parentName=self._parentName)
		self._pnts = pnts
		self._createdWithPoints=True
		self._api=api
		self._inactive = False
		if self._pnts is None :
			self._createdWithPoints=False
			self.create_edge_points()

	def writeToParameterFile(self,file) : 
		dictLoad=self.getDictFromYAMLFile(file=file)
		myDict=self.getDictInfo()
		fnd=False
		for indd, dictL in enumerate(dictLoad) :
			if dictL['name'] == myDict['name'] :
				dictLoad[indd] = myDict
				fnd=True
				break
		if not fnd :
			dictLoad.append(myDict)
		with open(file,"w") as openFile :
			yaml.dump(
				dictLoad,
				openFile,
				default_flow_style=False,
				indent=4,
				sort_keys=False,
				)

	@staticmethod
	def getDictFromYAMLFile(file) :

		try:
			with open(file, 'r') as fileYAML:
				unduDict = yaml.safe_load(fileYAML)
		except:
			return []
		return unduDict

	def getDictInfo(self) : 
		if self._createdWithPoints : 
			pntsD=[ pnts.tolist() for pnts in self._pnts]
		else :
			pntsD=None
		myDict={
			'type' : 'obj',
			'name' : self._name,
			'parentName' : self._parentName,
			'center' : self._center.tolist(),
		}
		if not (self._magnParas is None) :
			myDict.update({
				'len_x_main':self._magnParas._len_x_main(), 
				'len_y_main':self._magnParas._len_y_main(), 
				'len_z_main':self._magnParas._len_z_main(), 
				'magnetization':self._magnParas._magnetization(),
				'magn_unit_vec':self._magnParas._magn_unit_vec().tolist(),
				'segm_x':self._magnParas._segm_x(),
				'segm_y':self._magnParas._segm_y(),
				'segm_z':self._magnParas._segm_z(),
				'material_id':self._magnParas._material_id(), #"magnet","pole", 'NiCuFoil'
				'chamf':self._magnParas._chamf(),
				'frac_y':self._magnParas._frac_y(),
				'frac_z':self._magnParas._frac_z(),
				'pnts' : pntsD,
				'createdWithPoints' : self._createdWithPoints,
				})
		return myDict

	@staticmethod
	def loadFromDict(dictParas) : 

		if not (dictParas['type'] == 'obj') :
			return

		magnParasL=magParameters(
			len_x_main=dictParas['len_x_main'], 
			len_y_main=dictParas['len_y_main'], 
			len_z_main=dictParas['len_z_main'], 
			magnetization=dictParas['magnetization'],
			magn_unit_vec=np.array(dictParas['magn_unit_vec']),
			segm_x=dictParas['segm_x'],
			segm_y=dictParas['segm_y'],
			segm_z=dictParas['segm_z'],
			material_id=dictParas['material_id'], #"magnet","pole", 'NiCuFoil'
			chamf=dictParas['chamf'],
			frac_y=dictParas['frac_y'],
			frac_z=dictParas['frac_z'],
		)

		pnts=dictParas['pnts']
		if not (pnts is None) :
			pnts=[np.array(p) for p in pnts]
		block=undumagBlockObject(
			center=np.array(dictParas['center']),
			magnParas=magnParasL,
			pnts=pnts,
			name=dictParas['name'],
			parentName=dictParas['parentName'],
			)
		block._createdWithPoints=dictParas['createdWithPoints']
		return block

	def set_struct_name(self,parentName='') :
		if len(parentName) > 0 :
			self._parentName=parentName
			self._mother = parentName
			self._nameStruct=parentName+'_' + self._name
		else :
			self._nameStruct=self._name

	def set_api(self,api) :
		if not (api is None):
			self._api=api

	def get_copy(self,name=None,parentName=None) :
		if name is None :
			name=self._name
		if parentName is None :
			parentName=self._parentName
		newObj = undumagBlockObject(
			center=copy.deepcopy(self._center),
			magnParas=copy.deepcopy(self._magnParas),
			pnts=copy.deepcopy(self._pnts),
			name=name,
			parentName=parentName,
			api=self._api
			)
		newObj._createdWithPoints=self._createdWithPoints
		return newObj

	def find_all_names(self,names=None) :
		if names is None: 
			names = []
		names.append(self._nameStruct)
		return names

	def find_all_mag_blocks(self, mag_blocks = []) : 
		mag_blocks.append(self)
		return mag_blocks

	def set_inactive(self) :
		self._inactive = True

	def change_segm(self,segm_x,segm_y,segm_z,frac_y=None,frac_z=None) :
		if self._magnParas is None :
			return
		self._magnParas._segm_x.set(segm_x)
		self._magnParas._segm_y.set(segm_y)
		self._magnParas._segm_z.set(segm_z)
		if not (frac_y is None) :
			if self._magnParas._frac_y() < 1 :
				if frac_y > 1 :
					self._magnParas._frac_y.set(1/frac_y)
				else :
					self._magnParas._frac_y.set(frac_y)
			else :
				if frac_y < 1 :
					self._magnParas._frac_y.set(1/frac_y)
				else :
					self._magnParas._frac_y.set(frac_y)
		if not (frac_z is None) :
			if self._magnParas._frac_z() < 1 :
				if frac_z > 1 :
					self._magnParas._frac_z.set(1/frac_z)
				else :
					self._magnParas._frac_z.set(frac_z)
			else :
				if frac_z < 1 :
					self._magnParas._frac_z.set(1/frac_z)
				else :
					self._magnParas._frac_z.set(frac_z)


	def find_magn_name(self,names,fnd_list=None) :
		if fnd_list is None :
			fnd_list=[]
		not_fnd = False
		for name in names :
			if len(name) > 0 :
				if not (self._nameStruct.find(name) >= 0) :
					not_fnd=True
					break
		if not not_fnd :
			fnd_list.append(self)

	def set_magnetization(self,magnetization,magn_unit_vec) : 
		if self._magnParas is None :
			return
		self._magnParas._magnetization.set(magnetization)
		self._magnParas._magn_unit_vec.set(magn_unit_vec)

	def create_edge_points(self) : 
		if self._magnParas is None :
			return
		p1 = np.array([
			-self._magnParas._len_x_main()/2.0,
			-self._magnParas._len_y_main()/2.0,
			-self._magnParas._len_z_main()/2.0,
			])

		p2 = np.array([
			-self._magnParas._len_x_main()/2.0,
			+self._magnParas._len_y_main()/2.0,
			-self._magnParas._len_z_main()/2.0,
			])

		p3 = np.array([
			-self._magnParas._len_x_main()/2.0,
			-self._magnParas._len_y_main()/2.0,
			+self._magnParas._len_z_main()/2.0,
			])

		p4 = np.array([
			-self._magnParas._len_x_main()/2.0,
			+self._magnParas._len_y_main()/2.0,
			+self._magnParas._len_z_main()/2.0,
			])

		p5 = np.array([
			self._magnParas._len_x_main()/2.0,
			-self._magnParas._len_y_main()/2.0,
			-self._magnParas._len_z_main()/2.0,
			])

		p6 = np.array([
			self._magnParas._len_x_main()/2.0,
			+self._magnParas._len_y_main()/2.0,
			-self._magnParas._len_z_main()/2.0,
			])

		p7 = np.array([
			self._magnParas._len_x_main()/2.0,
			-self._magnParas._len_y_main()/2.0,
			+self._magnParas._len_z_main()/2.0,
			])

		p8 = np.array([
			self._magnParas._len_x_main()/2.0,
			+self._magnParas._len_y_main()/2.0,
			+self._magnParas._len_z_main()/2.0,
			])

		self._pnts = [p1,p2,p3,p4,p5,p6,p7,p8]

	def move_it(self,vec) : 
		"""
		changes center, i.e. relative coordinates
		also changes its own absolute position
		"""
		p_center = self._center+vec
		self._center = p_center
		self.set_absolute_position()

	def set_absolute_position(self,parentPos=np.array([0.0,0.0,0.0])) :
		"""
		for a singal magnet block, relative and absolute coordinates 
		are defined to coincide, makes life easier
		"""
		self._absCenter=copy.deepcopy(parentPos)+copy.deepcopy(self._center)

	def rotate(self,degrees,axis,plane='yz'):
		if self._magnParas is None :
			return
		self._pnts = rotate(pnts=self._pnts,degrees=degrees)
		new_center = rotate(
			pnts=[self._center],
			degrees=degrees,
			axis=axis,
			plane=plane
			)[-1]
		# pdb.set_trace()
		# new_abs_center = rotate(
		# 	pnts=[self._absCenter],
		# 	degrees=degrees,
		# 	# axis=axis,
		# 	plane=plane
		# 	)[-1]
		# self._absCenter = new_abs_center
		self._center=new_center
		# need to update the lengths here in order to get the right dimensions later on
		if degrees == 180 :
			self._magnParas._frac_y.set(1/self._magnParas._frac_y())
			self._magnParas._frac_z.set(1/self._magnParas._frac_z())
		degrees_rad = 2*math.pi*degrees/360
		if ( not (self._magnParas._len_y_main() is None)) and ( not (self._magnParas._len_z_main() is None)) :
			new_ly_z = abs(math.sin(degrees_rad)*self._magnParas._len_y_main())
			new_lz_z = abs(math.cos(degrees_rad)*self._magnParas._len_z_main())
			new_lz = max(new_ly_z,new_lz_z)

			new_ly_y = abs(math.cos(degrees_rad)*self._magnParas._len_y_main())
			new_lz_y = abs(math.sin(degrees_rad)*self._magnParas._len_z_main())
			new_ly = max(new_ly_y,new_lz_y)
			self._magnParas._len_y_main.set(new_ly)
			self._magnParas._len_z_main.set(new_lz)

	def mirror(self,coord='x') :
		if self._magnParas is None :
			return
		if coord == 'x' :  
			self._absCenter[0] = - self._absCenter[0]
			self._center[0] = - self._center[0]
			for ind, pnt in enumerate( self._pnts) : 
				self._pnts[ind][0] = -pnt[0]
		elif coord == 'y' :
			self._magnParas._frac_y.set(1/self._magnParas._frac_y())
			self._absCenter[1] = - self._absCenter[1]
			self._center[1] = - self._center[1]
			newMagn=copy.deepcopy(self._magnParas._magn_unit_vec())
			newMagn[0]=-newMagn[0]
			self._magnParas._magn_unit_vec.set(newMagn)
			for ind, pnt in enumerate( self._pnts) : 
				self._pnts[ind][1] = -pnt[1]
		elif coord == 'z' :
			self._magnParas._frac_z.set(1/self._magnParas._frac_z())
			self._absCenter[2] = - self._absCenter[2]
			self._center[2] = - self._center[2]
			for ind, pnt in enumerate( self._pnts) : 
				self._pnts[ind][2] = -pnt[2]
		self.set_absolute_position()

	def add_to_clc(self, magns_ignore = None,api=None) : 
		if api is None :
			api=self._api
			self._api._prog_paras.load_all_magnetic_material_files_from_folders()
		if self._inactive :
			return 
		if not (magns_ignore is None):
			for ignore in magns_ignore:
				if self._mother==ignore :
					return 

		materialClass=magnetic_material.get_material_from_list(
			theList=api._prog_paras.magnetic_materials(),
			material_id=self._magnParas._material_id(),
			)
		newList,materialIndex=magnetic_material.add_new_to_material_list(
				theList=api._prog_paras.magnetic_materials(),
				other=materialClass,
				)
		api._prog_paras.magnetic_materials.set(newList)

		my_txt = self.create_clc_txt(materialIndex=materialIndex,materialClass=materialClass)
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

	def create_clc_txt(self,materialIndex,materialClass) : 
		chamfer = False
		block_type = 'BlockChamf'
		if self._magnParas._chamf() is None :
			corners = True
			block_type = 'Corners'
		else :
			corners = False
		if materialClass._base_material_type() == 'magnet' :
			material = 'Magnet'
			color = 'ColorMag_Hybrid'
		elif materialClass._base_material_type() == 'pole' :
			material = 'Pole'
			color = 'ColorPol_Hybrid'
		elif materialClass._base_material_type() == 'NiCuFoil' :
			material = 'Pole'
			color = 'ColorFoil'
		txt = []

		txt.append(f'& {material}\n')
		txt.append(f'{block_type} {self._nameStruct} {self._mother} ${color}  !key, name, mother, color\n')
		# txt.append(f'0.0 0.0 0.0\n')
		txt.append(f'{self._absCenter[0]} {self._absCenter[1]} {self._absCenter[2]}\n')
		if materialClass._base_material_type() == 'magnet' :
			txt.append(f'{self._magnParas._magnetization()} {self._magnParas._magn_unit_vec()[0]} {self._magnParas._magn_unit_vec()[1]} {self._magnParas._magn_unit_vec()[2]} {materialIndex+1}  !length bc and comp. of magnetization, material index\n')
		elif materialClass._base_material_type() == 'pole' :
			txt.append(f'{materialIndex+1}  !material index\n')
		elif materialClass._base_material_type() == 'NiCuFoil' :
			txt.append(f'{materialIndex+1}  !material index\n')
		if not corners :
			txt.append(f'{self._magnParas._len_x_main()} {self._magnParas._len_y_main()} {self._magnParas._len_z_main()} {self._magnParas._chamf():.2f}  !dimensions \n')
		txt.append(f'{self._magnParas._segm_x()} {self._magnParas._segm_y()} {self._magnParas._segm_z()} {self._magnParas._frac_y()} {self._magnParas._frac_z()}  !segmentation\n')
		if corners :
			txt.append(f'{len(self._pnts)}\n')
			for pnt in self._pnts : 
				txt.append(f'{pnt[0]} {pnt[1]} {pnt[2]}\n')
		txt.append(' \n')
		return txt

	def get_max_extent(self,maxs=None,mins=None) : 
		"""
		determines the extent of this magnet block and compares to maxs and mins vals given, returns
		the max and min vals
		"""
		max_x = self._pnts[0][0]
		min_x = self._pnts[0][0]
		max_y = self._pnts[0][1]
		min_y = self._pnts[0][1]
		max_z = self._pnts[0][2]
		min_z = self._pnts[0][2]
		for pnt in self._pnts[1:] :
			if pnt[0] > max_x :
				max_x = pnt[0]
			elif pnt[0] < min_x :
				min_x = pnt[0]
			if pnt[1] > max_y :
				max_y = pnt[1]
			elif pnt[1] < min_y :
				min_y = pnt[1]
			if pnt[2] > max_z :
				max_z = pnt[2]
			elif pnt[2] < min_z :
				min_z = pnt[2]

		max_x = self._absCenter[0] + max_x
		min_x = self._absCenter[0] + min_x
		max_y = self._absCenter[1] + max_y
		min_y = self._absCenter[1] + min_y
		max_z = self._absCenter[2] + max_z
		min_z = self._absCenter[2] + min_z
		# max_x = self._center[0] + self._len_x/2.0
		# min_x = self._center[0] - self._len_x/2.0
		# max_y = self._center[1] + self._len_y/2.0
		# min_y = self._center[1] - self._len_y/2.0
		# max_z = self._center[2] + self._len_z/2.0
		# min_z = self._center[2] - self._len_z/2.0
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
		p_center = np.array([
			0.5*(maxs[0]+mins[0]),
			0.5*(maxs[1]+mins[1]),
			0.5*(maxs[2]+mins[2])
			])
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
	magn_unit_vec = np.array([1.0,0.0,0.0])
	if magn_string == '+x' : 
		magn_unit_vec = np.array([1.0,0.0,0.0])
	elif magn_string == '-x' : 
		magn_unit_vec = np.array([-1.0,0.0,0.0])
	if magn_string == '+y' : 
		magn_unit_vec = np.array([0.0,1.0,0.0])
	elif magn_string == '-y' : 
		magn_unit_vec = np.array([0.0,-1.0,0.0])
	if magn_string == '+z' : 
		magn_unit_vec = np.array([0.0,0.0,1.0])
	elif magn_string == '-z' : 
		magn_unit_vec = np.array([0.0,0.0,-1.0])
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
	magnets_obj = undumagObjectList(magnet_blocks=fnd_list)
	return magnets_obj
