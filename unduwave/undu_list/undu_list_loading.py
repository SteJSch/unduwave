"""

"""

from unduwave.unduwave_incl import *
from unduwave.undu_modules.undu_parameters import *
import unduwave.helpers.file_folder_helpers as f_h
from unduwave.quantities.quantities import *
# from unduwave import undu_blocks
# from unduwave import undulatorComponents

import unduwave.undu_modules.undu_blocks as undu_blocks
import unduwave.undu_modules.undu_magneticObjectGeometries as magneticObjectGeometries
import unduwave.undu_modules.undu_undulatorComponents as undulatorComponents

class unduListError(Exception) :
	pass

def rotate_xy(vecIn,angle,negative=False) :

	vecOut=np.array([
		math.cos(angle)*vecIn[0]-math.sin(angle)*vecIn[1],
		math.sin(angle)*vecIn[0]+math.cos(angle)*vecIn[1],
		vecIn[2]
		])

	return vecOut

def from_el_sequence_get_obj_list_length(
		sequence,
		objectDict,
		unduDict,
		center=np.array([0.0,0.0,0.0]),
		name_core='',
		) :
	objlist_tmps=[]
	all_objlist=[]
	len_objs=0.0
	all_lens=[]

	for ind_sq, el in enumerate(sequence) : 
		if isinstance(el,float) : 
			all_lens.append(el)
			all_objlist.append('dist')
		else :
			try : 
				el=float(el)
				all_lens.append(el)
				all_objlist.append('dist')
			except :
				if isinstance(el,str) : 
					if el in objectDict.keys() : 
						el=objectDict[el]['magn_object']
						p_center, maxs, mins=el.get_max_extent()
						el_len=maxs[0]-mins[0]	
						all_lens.append(el_len)
						objlist_tmps.append(el)
						all_objlist.append(el)
					elif el in unduDict['undulator'].keys() :
						el_val=unduDict['undulator'][el]
						if isinstance(el_val,float) :
							all_lens.append(el_val)
							all_objlist.append('dist')
						else :
							raise unduListError(f"The entry {el} of undulator dict, {el_val}, is not suitable to construct end_struct")
					else :
						raise unduListError(f"The entry {el} not found for constructing end_struct")

	len_objs=sum(all_lens)
	# creation and positioning
	onlyObjcs=[]
	numEl=0
	for ind, obj in enumerate(all_objlist) :
		if isinstance(obj,str) : 
			continue
		numEl=numEl+1
		lenBefore=sum(all_lens[0:ind])
		thisLen=all_lens[ind]
		realObj=obj.get_copy(name=f'{name_core}{numEl}')
		realObj.move_it(vec=center+np.array([
			-len_objs/2.0+lenBefore+thisLen/2.0,
			0.0,
			0.0,
			]))
		onlyObjcs.append(realObj)

	return len_objs, onlyObjcs

def get_set_std_val(thedict,stdVal,name=None,names=[]) :
	if not (name is None) :
		names=[name]
	vals=()
	for name in names :
		val=stdVal
		if not name in thedict.keys() : 
			thedict.update({name:stdVal})
		else:
			val=thedict[name]
		vals=vals+(val,)
	return vals

class undu_list_element :

	def __init__(self,
		listFile,
		listFolder=None,
		**kwargs,
		) :
		ht = os.path.split(listFile)
		self._listFile=ht[1]
		if listFolder is None :
			self._listFolder=ht[0]
		else:
			self._listFolder=listFolder
		self._material_file_folder=None

		self._fullFile=os.path.join(self._listFolder, self._listFile)

	def loadYAML(self,file=None) :
		if file is None :
			file = self._fullFile
		with open(file, 'r') as fileYAML:
			unduDict = yaml.safe_load(fileYAML)
		return unduDict

	def constructUndulator(self, 
			unduDict=None,
			gap=None,
			shift=None,
			center=None,
			useSymmetries=True,
			) :
		undulator=None
		try:
			if unduDict is None : 
				unduDict=self.loadYAML()

			unduType=unduDict['undulator']['type']

			if unduType == 'Planar-Hybrid' :
				undulator=self.constructPlanarHybridUndulator(
					unduDict=unduDict,
					gap=gap,
					center=center,
					useSymmetries=useSymmetries,
					)

			if unduType == 'APPLE' :
				undulator=self.constructAppleUndulator(
					unduDict=unduDict,
					gap=gap,
					shift=shift,
					center=center,
					useSymmetries=useSymmetries,
					)

			"""
			1) Check type - planar / elliptical
			Planar)
				1) If Hybrid, check for pole and magnet definitions
				2) Load the pole and magnet definitions, construct the elements

			"""
		except unduListError as e: 
			undulator=None
			print(e)

		# pdb.set_trace()
		return undulator

	def constructPlanarHybridUndulator(self,
			unduDict,
			gap=None,
			center=None,
			useSymmetries=True,
			) : 
		undulator=None
		try:
			if center is None :
				if not 'origin' in unduDict['general'].keys() : 
					center=np.array([0.0,0.0,0.0])
				else :
					center=unduDict['general']['origin']
					center=np.array([
						center[2],
						center[0],
						center[1],
						])

			if not 'name' in unduDict['undulator'].keys() : 
				raise unduListError("The name entry is missing from the list")
			name=unduDict['undulator']['name']

			if not 'accelerator' in unduDict['undulator'].keys() : 
				raise unduListError("The accelerator entry is missing from the list")
			accelerator=unduDict['undulator']['accelerator']

			if not 'periods' in unduDict['undulator'].keys() : 
				raise unduListError("The periods entry is missing from the list")
			nperiods=unduDict['undulator']['periods']

			if not 'minimum_gap' in unduDict['undulator'].keys() : 
				raise unduListError("The minimum_gap entry is missing from the list")
			minimum_gap=unduDict['undulator']['minimum_gap']
			if gap is None :
				gap=minimum_gap
			gap_range=[minimum_gap,100.0]

			if not 'symmetries' in unduDict['general'].keys() : 
				raise unduListError("The symmetries entry is missing from the list")
			symmetries=[]
			if useSymmetries :
				symmetries=unduDict['general']['symmetries']
			for ind,el in enumerate(symmetries) :
				if el == 'X' :
					symmetries[ind]='z'
				elif el=='Y':
					symmetries[ind]='y'
			if not 'period_length' in unduDict['undulator'].keys() : 
				raise unduListError("The period_length entry is missing from the list")
			period_length=unduDict['undulator']['period_length']

			(mag_pole_slit,)=get_set_std_val(
				name='shim_mp',
				thedict=unduDict['undulator'],
				stdVal=0.0
				)

			(pole_mag_slit,)=get_set_std_val(
				name='shim_pm',
				thedict=unduDict['undulator'],
				stdVal=0.0
				)

			if not 'magnetization_first_magnet' in unduDict['undulator'].keys() : 
				raise unduListError("The magnetization_first_magnet entry is missing from the list")
			magnetization_first_magnet=np.array(unduDict['undulator']['magnetization_first_magnet'])
			magnetization_first_magnet=np.array([
				magnetization_first_magnet[2],
				magnetization_first_magnet[0],
				magnetization_first_magnet[1],
				])
			magnetization_first_magnet=magnetization_first_magnet/np.linalg.norm(magnetization_first_magnet)

			magnetization_scnd_magnet=rotate_xy(
				vecIn=magnetization_first_magnet,
				angle=math.pi,
				negative=False
				)

			periodicMagnetizationSequence=np.array([
				np.array([ 0.0, 1.0, 0.0 ]),
				magnetization_first_magnet,
				np.array([ 0.0, -1.0,0.0 ]),
				magnetization_scnd_magnet,
				])

			if not 'geometry_info_file' in unduDict.keys() : 
				raise unduListError("The geometry_info_file entry is missing from the list")

			if 'material_file_folder' in unduDict.keys() : 
				material_file_folder=unduDict['material_file_folder']
				if len(material_file_folder) > 0 :
					self._material_file_folder=material_file_folder

			geometry_info_file = unduDict['geometry_info_file']
			fullGeometryFile=os.path.join(self._listFolder, geometry_info_file)
			geometry_info=self.loadYAML(file=fullGeometryFile)

			obj_names={ 
				'periodic_magnet',
				'end_magnet_1',
				'end_magnet_2',
				'periodic_pole',
				'end_pole',
			}

			objectDict=self.find_construct_magnetic_objects_from_dict(
				obj_names=obj_names,
				geometry_info=geometry_info,
				unduDict=unduDict,
				)

			magnetizationNorm=unduDict['periodic_magnet']['remanence']
			periodicMagnetizationSequence=[el*magnetizationNorm for el in periodicMagnetizationSequence]

			period_seq=['periodic_pole','shim_pm','periodic_magnet','shim_mp','periodic_pole','shim_pm','periodic_magnet','shim_mp']
			len_period, period_objects = from_el_sequence_get_obj_list_length(
				sequence=period_seq,
				objectDict=objectDict,
				unduDict=unduDict,
				center=np.array([0.0,0.0,0.0])
				)

			if not (len_period==period_length) :
				raise unduListError(f"Reconstructed period does not match given period! reconstructed={len_period}, given={period_length}")

			period=undulatorComponents.period(
				period_length=period_length/1e3,
				objects=period_objects, # list of objects making up one period
				center=np.array([0.0,0.0,0.0]),
				name='period1'
			)

			(end_struct_dict,)=get_set_std_val(
				name='end_struct',
				thedict=unduDict,
				stdVal={'type':'standard','sequence':[]}
				)

			mytype=end_struct_dict['type']
			if mytype == 'standard' : 
				sequence_up=['end_magnet_2','shim_mp','end_pole','shim_pm','end_magnet_1','shim_mp']
				sequence_down=['shim_mp','periodic_pole','shim_pm','end_magnet_1','shim_mp','end_pole','shim_pm','end_magnet_2']
			elif mytype == 'user' : 
				if 'sequence' in end_struct_dict.keys() : 
					sequence=end_struct_dict['sequence']
					sequence_down=sequence
					sequence_up=sequence
				elif 'sequence_up' in end_struct_dict.keys() : 
					sequence_up=end_struct_dict['sequence_up']
					if 'sequence_down' in end_struct_dict.keys() : 
						sequence_down=end_struct_dict['sequence_down']
					else :
						raise unduListError("The sequence_down entry is missing from end_struct list")
				else :
					raise unduListError("The sequence/up/down entries are wrong or missing from end_struct list")
			else :
				raise unduListError("The type entry from the end_struct element is unknown")

			len_upStrm, upstream_end_objects = from_el_sequence_get_obj_list_length(
				sequence=sequence_up,
				objectDict=objectDict,
				unduDict=unduDict,
				center=np.array([0.0,0.0,0.0])
				)
			len_dwnStrm, downstream_end_objects = from_el_sequence_get_obj_list_length(
				sequence=sequence_down,
				objectDict=objectDict,
				unduDict=unduDict,
				center=np.array([0.0,0.0,0.0])
				)

			endUpstream=undulatorComponents.endConfig(
				dist_period=0.0,
				objects=upstream_end_objects, # list of objects making up one period
				center=np.array([0.0,0.0,0.0]),
				name='endUS'
			)

			endDownstream=undulatorComponents.endConfig(
				dist_period=0.0,
				objects=downstream_end_objects, # list of objects making up one period
				center=np.array([0.0,0.0,0.0]),
				name='endDS'
			)

			row=undulatorComponents.row(
				pos='ll',
				period=period,
				downstream_end=endDownstream,
				upstream_end=endUpstream,
				period_length=period_length,
				center=np.array([0.0,0.0,0.0]),
				name='row1',
				nperiods=nperiods,
				periodicMagnetizationSequence=periodicMagnetizationSequence,
				api=None,
				)
			undulator=undulatorComponents.undulator(
				name=name,
				accelerator=accelerator,
				period_length=period_length,
				rows=[row],
				center=center,
				gap=gap,
				gap_range=gap_range,
				shift=0.0,
				symmetries=symmetries,
				construct_quadrants=['ll'],
				)

		except unduListError as e :
			undulator=None
			print(e)
		return undulator

	def constructAppleUndulator(self,
			unduDict,
			gap=None,
			shift=None,
			center=None,
			useSymmetries=True,
			) : 
		undulator=None
		try:
			if center is None :
				if not 'origin' in unduDict['general'].keys() : 
					center=np.array([0.0,0.0,0.0])
				else :
					center=unduDict['general']['origin']
					center=np.array([
						center[2],
						center[0],
						center[1],
						])

			if not 'name' in unduDict['undulator'].keys() : 
				raise unduListError("The name entry is missing from the list")
			name=unduDict['undulator']['name']

			(compensation,)=get_set_std_val(
				name='compensation',
				thedict=unduDict['undulator'],
				stdVal=False
				)

			if compensation:
				if not 'dist_fm_cm' in unduDict['undulator'].keys() : 
					raise unduListError("The dist_fm_cm entry is missing from the list")
				dist_fm_cm=unduDict['undulator']['dist_fm_cm']

				if not 'drop_cm' in unduDict['undulator'].keys() : 
					raise unduListError("The drop_cm entry is missing from the list")
				drop_cm=unduDict['undulator']['drop_cm']

			if not 'accelerator' in unduDict['undulator'].keys() : 
				raise unduListError("The accelerator entry is missing from the list")
			accelerator=unduDict['undulator']['accelerator']

			if not 'periods' in unduDict['undulator'].keys() : 
				raise unduListError("The periods entry is missing from the list")
			nperiods=unduDict['undulator']['periods']

			if not 'minimum_gap' in unduDict['undulator'].keys() : 
				raise unduListError("The minimum_gap entry is missing from the list")
			minimum_gap=unduDict['undulator']['minimum_gap']
			if gap is None :
				gap=minimum_gap
			gap_range=[minimum_gap,100.0]

			if not 'shift' in unduDict['undulator'].keys() : 
				shift=shift
			else:
				shift=unduDict['undulator']['shift']

			if not 'symmetries' in unduDict['general'].keys() : 
				raise unduListError("The symmetries entry is missing from the list")
			symmetries=[]
			if useSymmetries :
				symmetries=unduDict['general']['symmetries']
			for ind,el in enumerate(symmetries) :
				if el == 'X' :
					symmetries[ind]='z'
				elif el=='Y':
					symmetries[ind]='y'
			if not 'period_length' in unduDict['undulator'].keys() : 
				raise unduListError("The period_length entry is missing from the list")
			period_length=unduDict['undulator']['period_length']

			(shim_m,row_slit,keeper_slit,)=get_set_std_val(
				names=['shim_m','row_slit','keeper_slit'],
				thedict=unduDict['undulator'],
				stdVal=False
				)

			"""
			getting rows, finding three keepers if have to
			"""
			(rows_dicts,)=get_set_std_val(
				names=['rows'],
				thedict=unduDict,
				stdVal={
					'll' : {
						'sequence' : 
							['upstream_end','periodic','downstream_end']
							}
						},
				)

			# Create magnetization vector for ll-row

			if not 'magnetization_first_magnet' in unduDict['undulator'].keys() : 
				raise unduListError("The magnetization_first_magnet entry is missing from the list")
			magnetization_first_magnet=np.array(unduDict['undulator']['magnetization_first_magnet'])
			magnetization_first_magnet=np.array([
				magnetization_first_magnet[2],
				magnetization_first_magnet[0],
				magnetization_first_magnet[1],
				])
			magnetization_first_magnet=magnetization_first_magnet/np.linalg.norm(magnetization_first_magnet)

			magnetization_scnd_magnet=rotate_xy(
				vecIn=magnetization_first_magnet,
				angle=math.pi/2.0,
				)
			magnetization_thrd_magnet=rotate_xy(
				vecIn=magnetization_scnd_magnet,
				angle=math.pi/2.0,
				)
			magnetization_frth_magnet=rotate_xy(
				vecIn=magnetization_thrd_magnet,
				angle=math.pi/2.0,
				)

			periodicMagnetizationSequence=np.array([
				magnetization_first_magnet,
				magnetization_scnd_magnet,
				magnetization_thrd_magnet,
				magnetization_frth_magnet,
				])

			if not 'geometry_info_file' in unduDict.keys() : 
				raise unduListError("The geometry_info_file entry is missing from the list")

			if 'material_file_folder' in unduDict.keys() : 
				material_file_folder=unduDict['material_file_folder']
				if len(material_file_folder) > 0 :
					self._material_file_folder=material_file_folder

			geometry_info_file = unduDict['geometry_info_file']
			fullGeometryFile=os.path.join(self._listFolder, geometry_info_file)
			geometry_info=self.loadYAML(file=fullGeometryFile)

			"""create the end-magnet info"""

			unduDict.update({
				'end_mag_1' : copy.deepcopy(unduDict['periodic_magnet'])
				})
			unduDict['end_mag_1']['geometry']='end_mag_1'
			geometry_info.update({
				'end_mag_1' : copy.deepcopy(geometry_info['periodic_magnet'])
				})
			geometry_info['end_mag_1']['main_block_dimensions'][2] = \
				geometry_info['end_mag_1']['main_block_dimensions'][2]*0.75
			geometry_info['end_mag_1']['clamp_cutout'][2] = \
				geometry_info['end_mag_1']['clamp_cutout'][2]*0.75

			unduDict.update({
				'end_mag_2' : copy.deepcopy(unduDict['periodic_magnet'])
				})
			unduDict['end_mag_2']['geometry']='end_mag_2'
			geometry_info.update({
				'end_mag_2' : copy.deepcopy(geometry_info['periodic_magnet'])
				})
			geometry_info['end_mag_2']['main_block_dimensions'][2] = \
				geometry_info['end_mag_2']['main_block_dimensions'][2]*0.5
			geometry_info['end_mag_2']['clamp_cutout'][2] = \
				geometry_info['end_mag_2']['clamp_cutout'][2]*0.5

			unduDict.update({
				'end_mag_3' : copy.deepcopy(unduDict['periodic_magnet'])
				})
			geometry_info.update({
				'end_mag_3' : copy.deepcopy(geometry_info['periodic_magnet'])
				})
			unduDict['end_mag_3']['geometry']='end_mag_3'
			geometry_info['end_mag_3']['main_block_dimensions'][2] = \
				geometry_info['end_mag_3']['main_block_dimensions'][2]*0.25
			geometry_info['end_mag_3']['clamp_cutout'][2] = \
				geometry_info['end_mag_3']['clamp_cutout'][2]*0.25

			obj_names={
				'periodic_magnet',
				'compensation_magnet',
				'end_mag_1',
				'end_mag_2',
				'end_mag_3',
			}

			objectDict=self.find_construct_magnetic_objects_from_dict(
				obj_names=obj_names,
				geometry_info=geometry_info,
				unduDict=unduDict,
				)
			magnLen0=objectDict['periodic_magnet']['geo_info']['main_block_dimensions'][2]

			magnetizationNorm=unduDict['periodic_magnet']['remanence']
			if not isinstance(magnetizationNorm,list) :
				magnetizationNorm=[magnetizationNorm]
			lenMagnNorm=len(magnetizationNorm)

			periodicMagnetizationSequence=[el*magnetizationNorm[ind%lenMagnNorm] for ind,el in enumerate(periodicMagnetizationSequence)]

			three_seq=['periodic_magnet',shim_m,'periodic_magnet',shim_m,'periodic_magnet']
			period_seq=['periodic_magnet',shim_m,'periodic_magnet',keeper_slit,'periodic_magnet',shim_m,'periodic_magnet',keeper_slit]
			len_period, period_objects = from_el_sequence_get_obj_list_length(
				sequence=period_seq,
				objectDict=objectDict,
				unduDict=unduDict,
				center=np.array([0.0,0.0,0.0])
				)

			if not (len_period==period_length) :
				raise unduListError(f"Reconstructed period does not match given period! reconstructed={len_period}, given={period_length}")

			period=undulatorComponents.period(
				period_length=period_length/1e3,
				objects=period_objects, # list of objects making up one period
				center=np.array([0.0,0.0,0.0]),
				name='p1'
			)

			# # creating compensation parts if needed
			if compensation : 

				compMagnLen0=objectDict['compensation_magnet']['geo_info']['main_block_dimensions'][2]
				keeper_slit_comp=period_length/2.0-(2.0*compMagnLen0+shim_m)
				period_seq_c=['compensation_magnet',shim_m,'compensation_magnet',keeper_slit_comp,'compensation_magnet',shim_m,'compensation_magnet',keeper_slit_comp]
				len_period_comp, period_objects_comp = from_el_sequence_get_obj_list_length(
					sequence=period_seq_c,
					objectDict=objectDict,
					unduDict=unduDict,
					center=np.array([0.0,0.0,0.0])
					)

				if not (len_period_comp==period_length) :
					raise unduListError(f"Reconstructed compens. period does not match given period! reconstructed={len_period_comp}, given={period_length}")

				period_c=undulatorComponents.period(
					period_length=period_length/1e3,
					objects=period_objects_comp, # list of objects making up one period
					center=np.array([0.0,0.0,0.0]),
					name='pc1'
				)

				three_seq_c=['compensation_magnet',shim_m,'compensation_magnet',shim_m,'compensation_magnet']

			# creating the ends

			(end_struct_dict,)=get_set_std_val(
				name='end_struct',
				thedict=unduDict,
				stdVal={'type':'standard','sequence':[]}
				)

			mytype=end_struct_dict['type']
			if mytype == 'standard' : 
				sequence_up=['end_mag_3',0.5*magnLen0,'end_mag_2',0.5*magnLen0,'end_mag_1']
				# minus keeper_slit in sequence_down is for closing the keeper-slit of last period
				sequence_down=['end_mag_1',0.5*magnLen0,'end_mag_2',0.5*magnLen0,'end_mag_3']
			elif mytype == 'user' : 
				if 'sequence' in end_struct_dict.keys() : 
					sequence=end_struct_dict['sequence']
					sequence_down=sequence
					sequence_up=sequence
				elif 'sequence_up' in end_struct_dict.keys() : 
					sequence_up=end_struct_dict['sequence_up']
					if 'sequence_down' in end_struct_dict.keys() : 
						sequence_down=end_struct_dict['sequence_down']
					else :
						raise unduListError("The sequence_down entry is missing from end_struct list")
				else :
					raise unduListError("The sequence/up/down entries are wrong or missing from end_struct list")
			elif mytype == 'none' : 
				pass
			else :
				raise unduListError("The type entry from the end_struct element is unknown")

			rows=[]
			for key, item in rows_dicts.items() : 

				rowPos=key
				rowSeq=item['sequence']
				objcSeqUp=[]
				objcSeqDown=[]
				indPeriod=0
				try:
					indPeriod=rowSeq.index('periodic')
				except ValueError as e:
					print("No period found")

				for indEl,elSeq in enumerate(rowSeq) :
					if elSeq=='upstream_end' :
						elSeq=sequence_up
					elif elSeq=='downstream_end' :
						elSeq=sequence_down
					elif elSeq=='three_keeper' :
						elSeq=three_seq
					else :
						elSeq=[elSeq]
					if indEl<indPeriod:
						if len(objcSeqUp)>0:
							elSeq=[keeper_slit]+elSeq
						objcSeqUp=objcSeqUp+elSeq
					if indEl>indPeriod:
						if len(objcSeqDown)>0:
							elSeq=[keeper_slit]+elSeq
						objcSeqDown=objcSeqDown+elSeq

				len_upStrm, upstream_end_objects = from_el_sequence_get_obj_list_length(
					sequence=objcSeqUp,
					objectDict=objectDict,
					unduDict=unduDict,
					center=np.array([0.0,0.0,0.0]),
					)

				len_dwnStrm, downstream_end_objects = from_el_sequence_get_obj_list_length(
					sequence=objcSeqDown,
					objectDict=objectDict,
					unduDict=unduDict,
					center=np.array([0.0,0.0,0.0])
					)

				endUpstream=undulatorComponents.endConfig(
					dist_period=keeper_slit,
					objects=upstream_end_objects, # list of objects making up one period
					center=np.array([0.0,0.0,0.0]),
					name='endUS'
				)

				endDownstream=undulatorComponents.endConfig(
					dist_period=0.0,
					objects=downstream_end_objects, # list of objects making up one period
					center=np.array([0.0,0.0,0.0]),
					name='endDS'
				)

				row=undulatorComponents.row(
					pos=rowPos,
					period=period,
					downstream_end=endDownstream,
					upstream_end=endUpstream,
					period_length=period_length,
					center=np.array([0.0,0.0,0.0]),
					name='',
					nperiods=nperiods,
					periodicMagnetizationSequence=periodicMagnetizationSequence,
					api=None,
					)
				row.move_it(vec=np.array([0.0,0.0,-row_slit/2.0]))

				if rowPos[0]=='u' : # mirror at zx plane
					row.mirror(coord='y')
				if rowPos[1]=='r' : # mirror at yx plane
					row.mirror(coord='z')

				rows.append(row)

				if compensation :

					per_center, per_maxs, per_mins=period._magnet_blocks[0].get_max_extent()
					perLenY=per_maxs[1]-per_mins[1]
					perLenZ=per_maxs[2]-per_mins[2]
					com_center, com_maxs, com_mins=period_c._magnet_blocks[0].get_max_extent()
					comLenY=com_maxs[1]-com_mins[1]
					comLenZ=com_maxs[2]-com_mins[2]

					upSeq_c=[]
					downSeq_c=[]
					startIndPer=None
					for indEl,elSeq in enumerate(rowSeq) :

						if elSeq.find('_end') >= 0 :
							continue
						if elSeq=='three_keeper' :
							elSeq=three_seq_c
						elif elSeq.find('periodic') >= 0 :
							if startIndPer is None :
								startIndPer=indEl
							continue
						else :
							elSeq=[elSeq]
						if startIndPer is None :
							if len(upSeq_c)>0:
								elSeq=[keeper_slit_comp]+[elSeq]
							upSeq_c=upSeq_c+elSeq
						else:
							if len(downSeq_c)>0:
								elSeq=[keeper_slit_comp]+[elSeq]
							downSeq_c=downSeq_c+elSeq

					len_upStrm_c, upstream_end_objects_c1 = from_el_sequence_get_obj_list_length(
						sequence=upSeq_c,
						objectDict=objectDict,
						unduDict=unduDict,
						center=np.array([0.0,0.0,0.0]),
						)
					if len(upstream_end_objects_c1) > 0 :
						endUpstream_c1=undulatorComponents.endConfig(
							dist_period=keeper_slit,
							objects=upstream_end_objects_c1, # list of objects making up one period
							center=np.array([0.0,0.0,0.0]),
							name='endUS'
						)
					else :
						endUpstream_c1=None

					len_dwnStrm_c, downstream_end_objects_c1 = from_el_sequence_get_obj_list_length(
						sequence=downSeq_c,
						objectDict=objectDict,
						unduDict=unduDict,
						center=np.array([0.0,0.0,0.0])
						)

					if len(downstream_end_objects_c1) > 0 :
						endDownstream_c1=undulatorComponents.endConfig(
							dist_period=0.0,
							objects=downstream_end_objects_c1, # list of objects making up one period
							center=np.array([0.0,0.0,0.0]),
							name='endDS'
						)
					else :
						endDownstream_c1=None

					row_c_dwn=undulatorComponents.row(
						pos=rowPos,
						period=period_c,
						downstream_end=endDownstream_c1,
						upstream_end=endUpstream_c1,
						period_length=period_length,
						center=np.array([0.0,0.0,0.0]),
						name='cd',
						nperiods=nperiods,
						periodicMagnetizationSequence=periodicMagnetizationSequence,
						api=None,
						)

					my_p_center, my_maxs, my_mins=row_c_dwn.get_max_extent()
					row_c_dwn.rotate(
						degrees=-90,
						axis=my_p_center,
						plane='yz'
						)
					row_c_dwn.move_it(vec=-my_p_center)
					row_c_dwn.move_it(vec=np.array([0.0,-comLenZ/2.0,-comLenY/2.0]))

					row_c_dwn.move_it(vec=np.array([
						0.0,
						-perLenY-dist_fm_cm,
						row_slit/2.0 - drop_cm,
						]))

					len_upStrm_c, upstream_end_objects_c2 = from_el_sequence_get_obj_list_length(
						sequence=upSeq_c,
						objectDict=objectDict,
						unduDict=unduDict,
						center=np.array([0.0,0.0,0.0]),
						)
					if len(upstream_end_objects_c2) > 0 :
						endUpstream_c2=undulatorComponents.endConfig(
							dist_period=keeper_slit,
							objects=upstream_end_objects_c2, # list of objects making up one period
							center=np.array([0.0,0.0,0.0]),
							name='endUS'
						)
					else :
						endUpstream_c2=None

					len_dwnStrm_c, downstream_end_objects_c2 = from_el_sequence_get_obj_list_length(
						sequence=downSeq_c,
						objectDict=objectDict,
						unduDict=unduDict,
						center=np.array([0.0,0.0,0.0])
						)

					if len(downstream_end_objects_c2) > 0 :
						endDownstream_c2=undulatorComponents.endConfig(
							dist_period=0.0,
							objects=downstream_end_objects_c2, # list of objects making up one period
							center=np.array([0.0,0.0,0.0]),
							name='endDS'
						)
					else :
						endDownstream_c2=None

					row_c_lft=undulatorComponents.row(
						pos=rowPos,
						period=period_c,
						downstream_end=endDownstream_c2,
						upstream_end=endUpstream_c2,
						period_length=period_length,
						center=np.array([0.0,0.0,0.0]),
						name='cl',
						nperiods=nperiods,
						periodicMagnetizationSequence=periodicMagnetizationSequence,
						api=None,
						)

					row_c_lft.move_it(vec=np.array([
						0.0,
						- drop_cm,
						row_slit/2.0 - perLenZ - dist_fm_cm,
						]))

					if rowPos[0]=='u' : # mirror at zx plane
						row_c_dwn.mirror(coord='y')
						row_c_lft.mirror(coord='y')
					if rowPos[1]=='r' : # mirror at yx plane
						row_c_dwn.mirror(coord='z')
						row_c_lft.mirror(coord='z')

					rows.append(row_c_dwn)
					rows.append(row_c_lft)

			undulator=undulatorComponents.undulator(
				name=name,
				accelerator=accelerator,
				period_length=period_length,
				rows=rows,
				center=center,
				gap=gap,
				gap_range=gap_range,
				shift=shift,
				symmetries=symmetries,
				construct_quadrants=['ll','lr','ul','ur'],
				# construct_quadrants=['ul'],
				)

		except unduListError as e :
			undulator=None
			print(e)
		return undulator

	def copy_magnetic_object_dict(self,copy_dict,magn_obj_name=None) :

		new_dict={}
		for key, el in copy_dict.items() :
			if isinstance(el,dict) : 
				new_dict.update({ 
					key : self.copy_magnetic_object_dict(copy_dict=el) 
					})
			elif key == 'magn_object' :
				if magn_obj_name is None :
					magn_obj_name=el._name
				new_dict.update({key:el.get_copy(name=magn_obj_name)})
			elif key == 'geometry_name' :
				if magn_obj_name is None :
					magn_obj_name=el._name
				new_dict.update({key:magn_obj_name})
			else :
				new_dict.update({key:copy.deepcopy(el)})
		return new_dict

	def find_construct_magnetic_objects_from_dict(self,
			obj_names,
			geometry_info,
			unduDict
			) :

		resDict={}
		for indObj, objName in enumerate(obj_names) :
			key=objName
			resDict.update({key:{}})
			item=resDict[key]
			if not key in unduDict.keys() : 
				raise unduListError(f"The {key} entry is missing from the file")

			if not 'geometry' in unduDict[key].keys() : 
				raise unduListError(f"The geometry entry is missing from the {key}")

			item.update({'geometry_name' : unduDict[key]['geometry']})

			(y_drop,)=get_set_std_val(
				thedict=unduDict[key],
				stdVal=0.0,
				name='y_drop',
				)
			(y_alignment_0,)=get_set_std_val(
				thedict=unduDict[key],
				stdVal='gap',
				name='y_alignment_0',
				)

			item.update({'y_drop' : y_drop})
			item.update({'y_alignment_0' : y_alignment_0})

			if not item['geometry_name'] in geometry_info.keys() : 
				raise unduListError(f"The {item['geometry_name']} entry is missing from the geometry info file")

			geo_info = geometry_info[item['geometry_name']]
			item.update({'geo_info':geo_info})

			magn_object=self.constructMagneticObject(
					objectDict=item['geo_info'],
					material_info=unduDict[key],
					)
			item.update({'magn_object':magn_object})

		# Alignment after dict is complete
		for key, item in resDict.items() :

			aly=item['y_alignment_0']
			alignWith=None
			alignBottom=False
			if aly.find('bottom_') >= 0 :
				alignWith=aly.split('bottom_')[-1] 
				alignBottom=True
			elif aly.find('top_') >= 0 : 
				alignWith=aly.split('top_')[-1] 
			if not (alignWith is None) :

				if not (alignWith in resDict.keys()) :
					raise unduListError(f"Trying alignment of {info['name']}, but target {alignWith} is missing from unduDict")
				my_p_center, my_maxs, my_mins=item['magn_object'].get_max_extent()
				o_p_center, o_maxs, o_mins=resDict[alignWith]['magn_object'].get_max_extent()
				if alignBottom : 
					moveBy=np.array([0.0,o_mins[1]-my_mins[1],0.0])
				else :
					moveBy=np.array([0.0,o_maxs[1]-my_maxs[1],0.0])
				item['magn_object'].move_it(vec=moveBy)

		# finally, apply the drops
		for key, item in resDict.items() :
			moveBy=np.array([0.0,-item['y_drop'],0.0])			
			item['magn_object'].move_it(vec=moveBy)

		return resDict

	def constructMagneticObject(self,objectDict,material_info) :
		magnObject=None
		try:
			if not 'form_type' in objectDict.keys() : 
				raise unduListError("constructMagneticObject: form_type entry has to be present in geometry definition")
			if not 'material_id' in material_info.keys() : 
				raise unduListError("constructMagneticObject: material_id entry has to be present in definition")

			if objectDict['form_type'] == 'oneSidedClamps' : 
				magnObject=self.construct_oneSidedClamps_object(
					objectDict=objectDict,
					material_info=material_info,
					)
			elif objectDict['form_type'] == 'square' : 
				magnObject=self.construct_square_object(
					objectDict=objectDict,
					material_info=material_info,
					)
			elif objectDict['form_type'] == 'twoSidedClamps_a' : 
				magnObject=self.twoSidedClamps_a(
					objectDict=objectDict,
					material_info=material_info,
					)
			elif objectDict['form_type'] == 'twoSidedClamps_b' : 
				magnObject=self.twoSidedClamps_b(
					objectDict=objectDict,
					material_info=material_info,
					)
			elif objectDict['form_type'] == 'cpmuStdPole' : 
				magnObject=self.construct_cpmuStdPole_object(
					objectDict=objectDict,
					material_info=material_info,
					)
			else :
				raise unduListError("form_type not defined")
		except unduListError as e :
			magnObject=None
			print(e)
		return magnObject

	def twoSidedClamps_a(self,objectDict,material_info) : 

		obj=None
		"""
		check if material_info is all there and objectDict has all necessary 
		info
		"""
		try:
			if not 'main_block_dimensions' in objectDict.keys() : 
				raise unduListError("The main_block_dimensions entry is missing from the geometry info")
			if not 'clamp_cutout' in objectDict.keys() : 
				raise unduListError("The clamp_cutout entry is missing from the geometry info")
			if not 'fracs' in objectDict.keys() : 
				raise unduListError("The fracs entry is missing from the geometry info")
			if not 'segms' in objectDict.keys() : 
				raise unduListError("The segms entry is missing from the geometry info")
			if not 'chamf' in objectDict.keys() : 
				raise unduListError("The chamf entry is missing from the geometry info")

			fullDims=objectDict['main_block_dimensions']
			clampCut=objectDict['clamp_cutout']
			segms=objectDict['segms']
			fracs=objectDict['fracs']
			chamf=objectDict['chamf']
			if chamf == 0.0 : 
				chamf=None
			frac_z=fracs[0]
			if frac_z < 1 :
				frac_z=1

			magnParasMain=undu_blocks.magParameters(
				len_x_main=fullDims[2], 
				len_y_main=fullDims[1],
				len_z_main=fullDims[0], 
				segm_x=segms[2],
				segm_y=segms[1],
				segm_z=int(segms[0]),
				frac_y=fracs[1],
				frac_z=frac_z,
				chamf=chamf,
				material_id=material_info['material_id'],
			)

			magnParasClamp=undu_blocks.magParameters(
				len_x_main=clampCut[2], 
				len_y_main=clampCut[1],
				len_z_main=clampCut[0], 
				segm_x=2,
				segm_y=2,
				segm_z=2,
				frac_y=fracs[1],
				frac_z=frac_z,
				chamf=chamf,
				material_id=material_info['material_id'],
			)

			magnObject=magneticObjectGeometries.create_twoSidedClamps_a( 
				center=np.array([
					0.0,
					-magnParasMain._len_y_main()/2.0,
					-magnParasMain._len_z_main()/2.0,
					]), 
				magn_paras=magnParasMain,
				magn_paras_clamp=magnParasClamp,
				name=''
				)

			obj=magnObject

		except unduListError as e :
			obj=None
			print(e)
		return obj

	def twoSidedClamps_b(self,objectDict,material_info) : 

		obj=None
		"""
		check if material_info is all there and objectDict has all necessary 
		info
		"""
		try:
			if not 'main_block_dimensions' in objectDict.keys() : 
				raise unduListError("The main_block_dimensions entry is missing from the geometry info")
			if not 'clamp_cutout' in objectDict.keys() : 
				raise unduListError("The clamp_cutout entry is missing from the geometry info")
			if not 'fracs' in objectDict.keys() : 
				raise unduListError("The fracs entry is missing from the geometry info")
			if not 'segms' in objectDict.keys() : 
				raise unduListError("The segms entry is missing from the geometry info")
			if not 'chamf' in objectDict.keys() : 
				raise unduListError("The chamf entry is missing from the geometry info")

			fullDims=objectDict['main_block_dimensions']
			clampCut=objectDict['clamp_cutout']
			segms=objectDict['segms']
			fracs=objectDict['fracs']
			frac_z=fracs[0]
			if frac_z < 1 :
				frac_z=1
			chamf=objectDict['chamf']
			if chamf == 0.0 : 
				chamf=None

			magnParasMain=undu_blocks.magParameters(
				len_x_main=fullDims[2], 
				len_y_main=fullDims[1],
				len_z_main=fullDims[0], 
				segm_x=segms[2],
				segm_y=segms[1],
				segm_z=int(segms[0]),
				frac_y=fracs[1],
				frac_z=frac_z,
				chamf=chamf,
				material_id=material_info['material_id'],
			)

			magnParasClamp=undu_blocks.magParameters(
				len_x_main=clampCut[2], 
				len_y_main=clampCut[1],
				len_z_main=clampCut[0], 
				segm_x=2,
				segm_y=2,
				segm_z=2,
				frac_y=fracs[1],
				frac_z=frac_z,
				chamf=chamf,
				material_id=material_info['material_id'],
			)

			magnObject=magneticObjectGeometries.create_twoSidedClamps_b( 
				center=np.array([
					0.0,
					-magnParasMain._len_y_main()/2.0,
					-magnParasMain._len_z_main()/2.0,
					]), 
				magn_paras=magnParasMain,
				magn_paras_clamp=magnParasClamp,
				name=''
				)

			obj=magnObject

		except unduListError as e :
			obj=None
			print(e)
		return obj

	def construct_cpmuStdPole_object(self,objectDict,material_info) :

		obj=None
		"""
		check if material_info is all there and objectDict has all necessary 
		info
		"""
		try:
			if not 'main_block_dimensions' in objectDict.keys() : 
				raise unduListError("The main_block_dimensions entry is missing from the geometry info")
			if not 'lPolMSY' in objectDict.keys() : 
				raise unduListError("The lPolMSY entry is missing from the geometry info")
			if not 'lPolSZ' in objectDict.keys() : 
				raise unduListError("The lPolSZ entry is missing from the geometry info")
			if not 'lPolLSY' in objectDict.keys() : 
				raise unduListError("The lPolLSY entry is missing from the geometry info")
			if not 'lPolMChamf' in objectDict.keys() : 
				raise unduListError("The lPolMChamf entry is missing from the geometry info")
			if not 'fracs' in objectDict.keys() : 
				raise unduListError("The fracs entry is missing from the geometry info")
			if not 'segms' in objectDict.keys() : 
				raise unduListError("The segms entry is missing from the geometry info")
			if not 'chamf' in objectDict.keys() : 
				raise unduListError("The chamf entry is missing from the geometry info")

			fullDims=objectDict['main_block_dimensions']
			segms=objectDict['segms']
			fracs=objectDict['fracs']
			frac_z=int(fracs[0]/2)
			chamf=objectDict['chamf']
			if chamf == 0.0 : 
				chamf=None

			if frac_z < 1 :
				frac_z=1
			poleParasMain=undu_blocks.magParameters(
				len_x_main=fullDims[2], 
				len_y_main=fullDims[1],
				len_z_main=fullDims[0]/2.0, 
				segm_x=segms[2],
				segm_y=segms[1],
				segm_z=int(segms[0]/2),
				frac_y=fracs[1],
				frac_z=frac_z,
				chamf=chamf,
				material_id=material_info['material_id'],
			)
			poleParasLowerSide=undu_blocks.magParameters(
				len_x_main=poleParasMain._len_x_main(), 
				len_y_main=objectDict['lPolLSY'],
				len_z_main=objectDict['lPolSZ'], 
				segm_x=poleParasMain._segm_x(),
				material_id=material_info['material_id'],
			)
			poleParasUpperSide=undu_blocks.magParameters(
				len_x_main=poleParasMain._len_x_main(), 
				len_y_main=objectDict['lPolMSY'],
				len_z_main=objectDict['lPolSZ'], 
				segm_x=poleParasMain._segm_x(),
				material_id=material_info['material_id'],
			)

			magnObject=magneticObjectGeometries.create_cpmuStdPole_geometry(
				center=np.array([
					0.0,
					-poleParasMain._len_y_main()/2.0,
					- poleParasMain._len_z_main()/2.0
					]),
				poleParasMain=poleParasMain,
				poleParasLowerSide=poleParasLowerSide,
				poleParasUpperSide=poleParasUpperSide,
				cliffSide=objectDict['lPolMChamf'],
				name=''
				)
			# obj=[magnObject, poleParasMain, poleParasLowerSide, poleParasUpperSide]
			obj=magnObject

		except unduListError as e :
			obj=None
			print(e)
		return obj

	def construct_oneSidedClamps_object(self,objectDict,material_info) :

		obj=None
		"""
		check if material_info is all there and objectDict has all necessary 
		info
		"""
		try:
			if not 'main_block_dimensions' in objectDict.keys() : 
				raise unduListError("The main_block_dimensions entry is missing from the geometry info")
			if not 'clamp_cutout' in objectDict.keys() : 
				raise unduListError("The clamp_cutout entry is missing from the geometry info")
			if not 'fracs' in objectDict.keys() : 
				raise unduListError("The fracs entry is missing from the geometry info")
			if not 'segms' in objectDict.keys() : 
				raise unduListError("The segms entry is missing from the geometry info")
			if not 'chamf' in objectDict.keys() : 
				raise unduListError("The chamf entry is missing from the geometry info")

			fullDims=objectDict['main_block_dimensions']
			cutOut=objectDict['clamp_cutout']
			segms=objectDict['segms']
			fracs=objectDict['fracs']
			frac_z=int(fracs[0]/2)
			if frac_z < 1 :
				frac_z=1			
			chamf=objectDict['chamf']
			if chamf == 0.0 : 
				chamf=None			

			magnParasMain=undu_blocks.magParameters(
				len_x_main=fullDims[2], 
				len_y_main=fullDims[1],
				len_z_main=fullDims[0]/2.0-cutOut[0], 
				segm_x=segms[2],
				segm_y=segms[1],
				segm_z=int(segms[0]/2),
				frac_y=fracs[1],
				frac_z=frac_z,
				chamf=chamf,
				material_id=material_info['material_id'],
			)

			magnParasSide=undu_blocks.magParameters(
				len_x_main=cutOut[2], 
				len_y_main=fullDims[1]-cutOut[1],
				len_z_main=cutOut[0], 
				segm_x=segms[2],
				segm_y=1,
				segm_z=1,
				frac_y=1,
				frac_z=1,
				chamf=chamf,
				material_id=material_info['material_id'],
			)

			magnObject=magneticObjectGeometries.create_cpmuStdMagnet_geometry(
				center=np.array([
					0.0,
					-magnParasMain._len_y_main()/2.0,
					- magnParasMain._len_z_main()/2.0,
					]),
				magnParasMain=magnParasMain,
				magnParasSide=magnParasSide,
				name=''
				)

			# obj=[magnObject, magnParasMain, magnParasSide]
			obj=magnObject

		except unduListError as e :
			obj=None
			print(e)
		return obj

	def construct_square_object(self,objectDict,material_info) :

		obj=None
		"""
		check if material_info is all there and objectDict has all necessary 
		info
		"""
		try:
			if not 'main_block_dimensions' in objectDict.keys() : 
				raise unduListError("The main_block_dimensions entry is missing from the geometry info")
			if not 'fracs' in objectDict.keys() : 
				raise unduListError("The fracs entry is missing from the geometry info")
			if not 'segms' in objectDict.keys() : 
				raise unduListError("The segms entry is missing from the geometry info")
			if not 'chamf' in objectDict.keys() : 
				raise unduListError("The chamf entry is missing from the geometry info")

			fullDims=objectDict['main_block_dimensions']
			segms=objectDict['segms']
			fracs=objectDict['fracs']
			frac_z=int(fracs[0]/2)
			if frac_z < 1 :
				frac_z=1
			chamf=objectDict['chamf']
			if chamf == 0.0 : 
				chamf=None			

			magnParasMain=undu_blocks.magParameters(
				len_x_main=fullDims[2], 
				len_y_main=fullDims[1],
				len_z_main=fullDims[0]/2.0, 
				segm_x=segms[2],
				segm_y=segms[1],
				segm_z=int(segms[0]/2),
				frac_y=fracs[1],
				frac_z=frac_z,
				chamf=chamf,
				material_id=material_info['material_id'],
			)

			magnObject = undu_blocks.undumagBlockObject(
				center=np.array([
					0.0,
					-magnParasMain._len_y_main()/2.0,
					- magnParasMain._len_z_main()/2.0,
					]),
				magnParas=magnParasMain,
				name='',
				parentName='',
				)
			# obj=[magnObject,magnParasMain]
			obj=magnObject

		except unduListError as e :
			obj=None
			print(e)
		return obj
