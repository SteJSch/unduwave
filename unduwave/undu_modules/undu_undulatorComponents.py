"""

"""

from unduwave.unduwave_incl import *
from unduwave.undu_modules.undu_parameters import *
import unduwave.helpers.file_folder_helpers as f_h
from unduwave.quantities.quantities import *
import unduwave.undu_modules.undu_blocks as undu_blocks

class unduComponentsError(Exception) :
	pass

class period(undu_blocks.undumagObjectList) :

	def __init__(self,
		period_length,
		objects, # list of objects making up one period
		center,
		api=None,
		**kwargs,
		) :
		super(period, self).__init__(
			magnet_blocks=objects,
			api=api,
			center=center,
			**kwargs,
			)
		self._period_length=period_length

	def getDictInfo(self) : 

		periodDict={
		'period_length' : self._period_length, # distance btwn center and begin/end of periodic part
		}
		periodDict.update(super(period, self).getDictInfo())
		return periodDict

	@staticmethod
	def loadFromDict(dictParas) : 

		periodList=undu_blocks.undumagObjectList.loadFromDict(dictParas=dictParas)
		nperiod=period(
			name=dictParas['name'],
			parentName=dictParas['parentName'],
			center=np.array(dictParas['center']),
			period_length=dictParas['period_length'],
			objects=periodList._magnet_blocks, # list of objects making up one period
			)
		return nperiod

class endConfig(undu_blocks.undumagObjectList) :

	def __init__(self,
		dist_period,
		objects, # list of objects making up one period
		center,
		api=None,
		**kwargs,
		) :
		super(endConfig, self).__init__(
			magnet_blocks=objects,
			api=api,
			center=center,
			**kwargs,
			)
		self._dist_period=dist_period

	def getDictInfo(self) : 

		endConfigDict={
		'dist_period' : self._dist_period, # distance btwn center and begin/end of periodic part
		}
		endConfigDict.update(super(endConfig, self).getDictInfo())
		return endConfigDict

	@staticmethod
	def loadFromDict(dictParas) : 

		endBlcks=undu_blocks.undumagObjectList.loadFromDict(dictParas=dictParas)
		end=endConfig(
			name=dictParas['name'],
			parentName=dictParas['parentName'],
			center=np.array(dictParas['center']),
			dist_period=dictParas['dist_period'],
			objects=endBlcks._magnet_blocks, # list of objects making up one period
			)
		return end

class row(undu_blocks.undumagObjectList) :

	def __init__(self,
		period,
		center,
		period_length,
		nperiods,
		periodicMagnetizationSequence,
		pos='',
		downstream_end=None,
		upstream_end=None,
		api=None,
		**kwargs,
		) :
	
		self._pos=pos
		self._period=period
		self._downstream_end=downstream_end
		self._upstream_end=upstream_end
		self._period_length=period_length 
		self._nperiods=nperiods
		self._periodicMagnetizationSequence=periodicMagnetizationSequence
		self._periodicMagnetizationSequenceNorm=[]

		self._periodicMagnetizationSequenceNorm=self.magnNormUnitVecFromVecs(
				periodicMagnetizationSequence=self._periodicMagnetizationSequence
				)

		blocks=self.positionMagnetizeEnds(center=center)

		super(row, self).__init__(
			magnet_blocks=blocks,
			api=api,
			center=center,
			**kwargs,
			)
		self._center=copy.deepcopy(center)
		self._absCenter=copy.deepcopy(center)

	def positionMagnetizeEnds(self,center) : 

		blocks=[]
		startMagnNum=0

		lenPeriodic=self._period_length*self._nperiods
		lenUpX=0.0
		lenDownX=0.0
		if not ( self._upstream_end is None ) :
			centerUp, minsUp, maxsUp = self._upstream_end.get_max_extent()
			lenUpX=abs(maxsUp[0]-minsUp[0])+self._upstream_end._dist_period
		if not ( self._downstream_end is None ) :
			centerDown, minsDown, maxsDown = self._downstream_end.get_max_extent()
			lenDownX=abs(maxsDown[0]-minsDown[0])+self._downstream_end._dist_period
		totalLen=lenUpX+lenDownX+lenPeriodic

		endPeriodUpX=center[0]-totalLen/2.0
		numMagntznElements=len(self._periodicMagnetizationSequence)
		if not ( self._upstream_end is None ) :

			# numUpstreamEndObjcts=len(self._upstream_end._magnet_blocks)
			# startMagnNum=numUpstreamEndObjcts % numMagntznElements
			# startMagnNum=numMagntznElements-startMagnNum
			startMagnNum=0

			endPeriodUpX=endPeriodUpX+lenUpX/2.0
			self._upstream_end.move_it(vec=np.array([endPeriodUpX-self._upstream_end._center[0],0.0,0.0]))

			for obj in self._upstream_end._magnet_blocks :
				magnNow=self._periodicMagnetizationSequenceNorm[startMagnNum%numMagntznElements]
				obj.set_magnetization(magnetization=magnNow['magnLen'],magn_unit_vec=magnNow['magnUnitVec'])
				startMagnNum=startMagnNum+1

			blocks.append(self._upstream_end)

		center_period=endPeriodUpX+lenUpX/2.0+lenPeriodic/2.0
		center_vec_period=copy.deepcopy(center)
		center_vec_period[0]=center_period
		periods, startMagnNum=self.createNperiods(center=center_vec_period,startMagnNum=startMagnNum)
		for periodT in periods:
			blocks.append(periodT) 

		if not ( self._downstream_end is None ) :

			endPeriodDownX=center_period+lenPeriodic/2.0+lenDownX/2.0
			self._downstream_end.move_it(vec=np.array([endPeriodDownX-self._downstream_end._center[0],0.0,0.0]))

			for obj in self._downstream_end._magnet_blocks :
				magnNow=self._periodicMagnetizationSequenceNorm[startMagnNum%numMagntznElements]
				obj.set_magnetization(magnetization=magnNow['magnLen'],magn_unit_vec=magnNow['magnUnitVec'])
				startMagnNum=startMagnNum+1

			blocks.append(self._downstream_end)

		return blocks

	def magnNormUnitVecFromVecs(self,periodicMagnetizationSequence) : 
		periodicMagnetizationSequenceNorm=[]
		for magnetization in periodicMagnetizationSequence :
			norm=0
			for el in magnetization :
				norm=norm+el**2
			norm=math.sqrt(norm)
			normMagnetization=magnetization/norm
			periodicMagnetizationSequenceNorm.append({'magnLen':norm,'magnUnitVec':normMagnetization})
		return periodicMagnetizationSequenceNorm

	def createNperiods(self,center=None,startMagnNum=0) :
		periodList=[]
		if center is None :
			center=self._center
		totalLength=self._nperiods*self._period_length
		numMagntznElements=len(self._periodicMagnetizationSequence)

		for i in range(self._nperiods) :
			period=self._period.get_copy(name=f'period{i+1}')

			for obj in period._magnet_blocks :
				magnNow=self._periodicMagnetizationSequenceNorm[startMagnNum%numMagntznElements]
				obj.set_magnetization(magnetization=magnNow['magnLen'],magn_unit_vec=magnNow['magnUnitVec'])
				startMagnNum=startMagnNum+1

			posX=center[0]-totalLength/2.0+self._period_length/2.0+\
				self._period_length*i
			period.move_it(vec=np.array([
				posX-period._center[0],
				0.0,
				0.0,
				]))
			periodList.append(period)
		return periodList, startMagnNum

	def getDictInfo(self) : 

		rowDict=super(row, self).getDictInfo(onlyBase=True)

		downstream_end_dict={}
		if not ( self._downstream_end is None ) :
			downstream_end_dict=self._downstream_end.getDictInfo()

		upstream_end_dict={}
		if not ( self._upstream_end is None ) :
			upstream_end_dict=self._upstream_end.getDictInfo()

		period_dict=self._period.getDictInfo()

		rowDict.update({
			'pos' : self._pos,
			'period_length' : self._period_length,
			'nperiods' : self._nperiods,
			'periodicMagnetizationSequence' : [el.tolist() for el in self._periodicMagnetizationSequence],
			'periodic_part' : period_dict,
			'downstream_end' : downstream_end_dict,
			'upstream_end' : upstream_end_dict,
		})
		return rowDict

	@staticmethod
	def loadFromDict(dictParas,nperiods=None) : 
		nperiod=period.loadFromDict(dictParas['periodic_part'])
		downstream_end=endConfig.loadFromDict(dictParas['downstream_end'])
		upstream_end=endConfig.loadFromDict(dictParas['upstream_end'])
		if nperiods is None :
			nperiods=dictParas['nperiods']

		nrow = row(
			pos=dictParas['pos'],
			name=dictParas['name'],
			parentName=dictParas['parentName'],
			period=nperiod,
			downstream_end=downstream_end,
			upstream_end=upstream_end,
			center=np.array(dictParas['center']),
			period_length=dictParas['period_length'],
			nperiods=nperiods,
			periodicMagnetizationSequence=[np.array(el) for el in dictParas['periodicMagnetizationSequence']],
		)
		return nrow

class undulator(undu_blocks.undumagObjectList) :

	def __init__(self,
		rows,
		gap,
		gap_range=None,
		shift=0.0,
		shifts=None,
		beff=None,
		shift_range=np.array([0.0,0.0]),
		center=np.array([0.0,0.0,0.0]),
		symmetries=[],
		construct_quadrants=['ll','lr','ul','ur'],
		api=None,
		name='',
		accelerator='',
		period_length=None,
		nperiods=None,
		**kwargs,
		) :


		try:

			self._rows=rows
			self._beff=beff
			if gap_range is None :
				gap_range=np.array([gap,10*gap])
			self._gap_range=gap_range
			self._gap=gap
			self._shift=shift
			self._symmetries=symmetries
			self._accelerator=accelerator
			self._shift_range=shift_range
			self._period_length=period_length
			self._nperiods=nperiods

			period_lengths=[]
			nperiods=[]

			if not (shifts is None) : 
				shfts=shifts
				shift=shifts[int(len(shfts)/2)]
			else :
				shfts = [0.0,0.0,0.0,0.0]
				shift=shift*self._period_length
				if shift < 0 :
					shfts = [-shift,0.0,0.0,shift]
				elif shift > 0 :
					shfts = [shift,0.0,0.0,shift]

			quadrants={
					'll_rows':[],'ll_obj' : None,'ll_shifts' : np.array([shfts[2],0.0,0.0]),
					'lr_rows':[],'lr_obj' : None,'lr_shifts' : np.array([shfts[3],0.0,0.0]),
					'ul_rows':[],'ul_obj' : None,'ul_shifts' : np.array([shfts[0],0.0,0.0]),
					'ur_rows':[],'ur_obj' : None,'ur_shifts' : np.array([shfts[1],0.0,0.0]),
				}
			for row in self._rows :
				pos=row._pos
				quadrants[pos+'_rows'].append(row)
			for key, quadrant in quadrants.items() :
				if key.find('_rows')>=0:
					continue
				if key.find('_shifts')>=0:
					continue
				pos=key.split('_obj')[0]
				if len(quadrants[pos+'_rows']) < 1 :
					continue
				quadrants[pos+'_obj']=undu_blocks.undumagObjectList(
						magnet_blocks=quadrants[pos+'_rows'],
						api=api,
						center=center,
						**kwargs,
						)

			blocks=[]

			for quad_pos in construct_quadrants :

				quadr_key=quad_pos+'_obj'
				if not (quadrants[quadr_key] is None) :
					thequadrant=quadrants[quadr_key].get_copy(name=quad_pos+'R')
					if quad_pos[0] == 'u' :
						thequadrant.move_it(vec=np.array([0.0,self._gap/2.0,0.0]))
					elif quad_pos[0] == 'l' :
						thequadrant.move_it(vec=np.array([0.0,-self._gap/2.0,0.0]))
				else :
					base_quadr_key='ll_obj'
					if not (quadrants[base_quadr_key] is None) :
						thequadrant=quadrants[base_quadr_key].get_copy(name=quad_pos+'R')
					else:
						raise unduComponentsError(f"construct undulator: quadrant base {base_quadr_key} is not defined.")

					thequadrant.move_it(vec=np.array([0.0,-self._gap/2.0,0.0]))
					if quad_pos[0]=='u' : # mirror at zx plane
						thequadrant.mirror(coord='y')
					if quad_pos[1]=='r' : # mirror at yx plane
						thequadrant.mirror(coord='z')

				if not (shift == 0.0) :
					thequadrant.move_it(vec=quadrants[quad_pos+'_shifts'])

				blocks.append(thequadrant)

			super(undulator, self).__init__(
				magnet_blocks=blocks,
				api=api,
				center=center,
				name=name,
				**kwargs,
				)
		except unduComponentsError as e: 
			print(e)

	def createUndulator(self) :
		periodList=[]
		return periodList

	def getDictInfo(self) : 

		unduDict=super(undulator, self).getDictInfo(onlyBase=True)

		rowDicts=[]
		for row in self._rows :
			rowDict=row.getDictInfo()
			rowDicts.append(rowDict)

		unduDict.update({
			'accelerator' : self._accelerator,
			'gap_range' : self._gap_range,
			'period_lengths' : self._period_length,
			'nperiods' : self._nperiods,
			'shift_range' : self._shift_range,
			'symmetries': self._symmetries,
			'beff' : self._beff,
			'geometry' : {
				'rows' : rowDicts,
			},
		})
		return unduDict

	@staticmethod
	def loadFromDict(name,file=None,dictParas=None,gap=None,nperiods=None,shift=None,shifts=None,symmetries=None,api=None) :
		if file is None : 
			if dictParas is None :
				return
		else :
			dictParas=undulator.getDictFromYAMLFile(file)

		fnd=False
		for el in dictParas :
			if el['name'] == name :
				dictParas=el
				fnd=True
				break
		if not fnd:
			return
		if gap is None :
			gap=dictParas['gap_range'][0]
		if shift is None :
			if shifts is None :
				shift=dictParas['shift_range'][0]
			else :
				shift=0.0
		if nperiods is None :
			nperiods=dictParas['nperiods']

		geo=dictParas['geometry']
		rowDicts=geo['rows']
		rows=[]
		for rowDict in rowDicts :
			rowl=row.loadFromDict(
				dictParas=rowDict,
				nperiods=nperiods,
				)
			rows.append(rowl)
		if (symmetries is None) : 
			symmetries=dictParas['symmetries']

		undulatorLoad=undulator(
			rows=rows,
			period_length=dictParas['period_lengths'],
			beff=dictParas['beff'],
			gap=gap,
			gap_range=dictParas['gap_range'],
			shift=shift,
			shift_range=dictParas['shift_range'],
			center=np.array(dictParas['center']),
			symmetries=symmetries,
			api=api,
			name=dictParas['name'],
			parentName=dictParas['parentName'],
			accelerator=dictParas['accelerator'],
			shifts=shifts,
			# nperiods=nperiods,
		)
		return undulatorLoad
