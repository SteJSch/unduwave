"""
Contains the wave_from_b class that incorporates the API from python to WAVE and the function create_wave_instance that 
returns an instance of that class
"""

from unduwave.unduwave_incl import *
from unduwave.undu_modules.undu_parameters import *
import unduwave.helpers.file_folder_helpers as f_h
from unduwave.quantities.quantities import *
from unduwave import undu_magnets

class magnetRepresentation : 
	def __init__(self,
			center,
			magnParas,
			magnFun=None
			) :
		self._magnFun=magnFun
		self._magnParas=magnParas
		self._center=center

	def createMagnet(self,
			center=None,
			magnParas=None,
			name='mag'
			) :
		if center is None :
			center=self._center
		if magnParas is None :
			magnParas=self._magnParas
		centerTmp = copy.deepcopy(center)
		magnParasTmp = copy.deepcopy(magnParas)
		if self._magnFun is None : 
			return self.createEasyMagnet(magnParasTmp,magnParasTmp,name=name)
		return self._magnFun(center,magnParasTmp,name=name)

	def createEasyMagnet(self,
			center,
			magnParas,
			name='mag'
			) : 
		name_main = combine_strings(str_add='main',str_add_to=name)
		main_block = undu_magnets.undu_magnet_block_coords(
			p_center=p_main,
			name=name_main,
			magnetParameters=magnParas,
			mother=name,
			)
		return undu_magnets.undu_magnets(magnet_blocks=[main_block])

class unduRepresentation() : 
	def __init__( self,
			gap = 0.0, 
			shift=0.0,
			glueSlit = 0.0, 
			keeperSlit=0.0, 
			rowSlit=0.0, 
			nperiods =1,
			unduType='planar', # "planar", 'ellipt',
			magnetsPerKeeper=[],
			endMagnets=[],
			nKeeperPerPeriod=2,
			) : 
		"""
		Defines the standard parameters of the ivue32 constellation
		gap - gap
		glueSlit - slit between 2 magnets on same keeper
		keeperSlit - slit between 2 keepers on same row
		rowSlit - slit between 2 rows on same horizontal plane
		"""
		self._gap = gap
		self._glueSlit = glueSlit
		self._keeperSlit = keeperSlit
		self._rowSlit = rowSlit
		self._unduType=unduType
		self._nperiods=nperiods
		self._shift=shift
		self._magnetsPerKeeper=magnetsPerKeeper # defines the createMagnet funs to use in order
		self._nKeeperPerPeriod=nKeeperPerPeriod
		self._endMagnets=endMagnets

		# if magnFuns==[] :
		# 	magFuns=[self.create_easy_magnet]

	def getPeriodLength(self) : 
		length_keeper= self.getKeeperLength()
		per_l=self._nKeeperPerPeriod*length_keeper
		return per_l

	def getKeeperLength(self) : 
		length_keeper=0.0
		for magnC in self._magnetsPerKeeper : 
			length_keeper=length_keeper+magnC._magnParas._len_x_main
		length_keeper=length_keeper+(len(self._magnetsPerKeeper)-1)*self._glueSlit+self._keeperSlit
		return length_keeper

	def get_end_magn_paras(self,endNum,magnParas):
		"""
		Definition is such that that func. magn: (ll-row)
		 ---   
		| |  ->main
		|---  ->side1 (|)
		|-    ->side1 (|)

		 ^ - side2
		with y:|^ and z: ->
		"""
		endMagnParas=copy.deepcopy(magnParas)
		shrink=1.0
		if endNum==1:
			shrink=3.0/4.0
		elif endNum==2:
			shrink=1.0/2.0
		elif endNum==3:
			shrink=1.0/4.0

		endMagnParas._len_x_main=endMagnParas._len_x_main*shrink

		return endMagnParas

	def createKeeper(self,
			center, 
			magnetizations, 
			magnSeq, 
			startNum=0,
			nameCore='keeper_1',
			) :
		"""
		returns undu_magnets object with the magnets ordered
		"""
		length_keeper= self.getKeeperLength()
		centerTmp=copy.deepcopy(center)

		pX = centerTmp._x-length_keeper/2.0+self._keeperSlit/2.0
		pY0 = centerTmp._y
		pZ0 = centerTmp._z

		magnet_blocks=[]

		for indMagn, magnetC in enumerate( self._magnetsPerKeeper ):
			pX=pX+magnetC._magnParas._len_x_main/2.0
			pY=pY0+magnetC._center._y
			pZ=pZ0+magnetC._center._z
			magnCenter = undu_magnets.point_coords(x=pX,y=pY,z=pZ)
			name = combine_strings(str_add=f'obj_{indMagn+1}',str_add_to=nameCore)
			magnetC._magnParas._magnetization=magnetizations[(startNum+indMagn)%len(magnetizations)]
			magnetC._magnParas._magn_unit_vec=undu_magnets.create_magnetization_unit_vec(
				magn_string=magnSeq[(startNum+indMagn)%len(magnSeq)],
				)
			magnet = magnetC.createMagnet(center=magnCenter,name=name)
			magnet_blocks.append(magnet)
			pX=pX+magnetC._magnParas._len_x_main/2.0+self._glueSlit
			# print(f"The startNum {startNum} indMagn {indMagn} magnetization vec is {magnSeq[(startNum+indMagn)%len(magnSeq)]}")
			# pdb.set_trace()
		return undu_magnets.undu_magnets(magnet_blocks=magnet_blocks)

	def createPeriod(self,
			center, 
			magnetizations, 
			magnSeq, 
			nameCore='period_1',
			startNum=0,
			) :

		length_keeper= self.getKeeperLength()
		per_l= self.getPeriodLength()

		# pX=p_center._x-length_keeper/2.0+magnParas._len_x_main/2.0+undulatorP._keeperSlit/2.0

		centerTmp=copy.deepcopy(center)
		pX=centerTmp._x-per_l/2.0+length_keeper/2.0
		pY = centerTmp._y
		pZ = centerTmp._z
		nObjPerKeeper=len(self._magnetsPerKeeper)
		magnet_blocks=[]

		for keeperNum in range(1,self._nKeeperPerPeriod+1) :

			keepCenter = undu_magnets.point_coords(x=pX,y=pY,z=pZ)
			name = combine_strings(str_add=f'keeper_{keeperNum}',str_add_to=nameCore)
			keeper=self.createKeeper( 
					center=keepCenter, 
					magnetizations=magnetizations, 
					magnSeq=magnSeq,
					startNum=startNum,
					nameCore=name,
					)
			magnet_blocks.append(keeper)
			startNum=startNum+nObjPerKeeper
			pX=pX+length_keeper
		return undu_magnets.undu_magnets(magnet_blocks=magnet_blocks)

	def createNPeriods(self,
			center, 
			magnetizations, 
			magnSeq, 
			nameCore='period_1',
			) :
		per_l= self.getPeriodLength()

		p_center_per = copy.deepcopy(center)
		center_x_0 = center._x - per_l*self._nperiods/2.0 + per_l/2.0
		half_per = self._nperiods - int(self._nperiods)
		npers = self._nperiods
		half = False
		if half_per == 0.5 :
			half = True
			npers = int(npers)
		blocks = []
		last_n = 0
		startNum=0
		for n_ind in range(1,npers+1) : 
			p_center_per._x = center_x_0 + (n_ind-1)*per_l
			name = combine_strings(str_add=f'period_{n_ind}',str_add_to=nameCore)
			per=self.createPeriod(
					center=p_center_per, 
					magnetizations=magnetizations, 
					magnSeq=magnSeq, 
					nameCore=name,
					startNum=startNum,
					)
			blocks = blocks + per._magnet_blocks
			last_n = n_ind
			startNum=startNum+len(self._magnetsPerKeeper)+self._nKeeperPerPeriod
		if half :
			p_center_per._x = p_center_per._x + per_l*0.75 
			name = combine_strings(str_add=f'period_{npers+1}_keeper_1',str_add_to=nameCore)
			per=createKeeper( 
					center=p_center_per, 
					magnetizations=magnetizations, 
					magnSeq=magnSeq,
					nameCore=name,
					startNum=startNum,
					)
			blocks = blocks + per._magnet_blocks
		return undu_magnets.undu_magnets(magnet_blocks=blocks)

	def createEndStructElliptStd(self,
			center, 
			location, 
			magnets,
			magnetizations, 
			startNum,
			magnSeq, 
			nameCore='',
			) :
		"""
		creates the endstructure
		type='us' or 'ds'
		"""

		pX = center._x
		pY0 = center._y
		pZ0 = center._z

		magnet_blocks=[]
		if location=='us':
			facts=[0.3,0.3,0.75]
			magnet1=magnets[startNum%len(magnets)]
			magnet2=magnets[(startNum+1)%len(magnets)]
			magnet3=magnets[(startNum+2)%len(magnets)]
			pX=pX-magnet2._magnParas._len_x_main/2.0-self._glueSlit-magnet1._magnParas._len_x_main
		if location=='ds':
			facts=[0.75,0.75,0.3,0.3]
			magnet1=magnets[startNum%len(magnets)]
			magnet2=magnets[(startNum+1)%len(magnets)]
			pX=pX-magnet2._magnParas._len_x_main-1.5*self._glueSlit-magnet1._magnParas._len_x_main
		newStart=startNum
		for indMagn, fact in enumerate( facts ):
			magnet=magnets[(startNum+indMagn)%len(magnets)]
			magnetC=copy.deepcopy(magnet)
			magnetL=magnetC._magnParas._len_x_main
			pX=pX+magnetL/2.0
			pY=pY0+magnetC._center._y
			pZ=pZ0+magnetC._center._z
			magnCenter = undu_magnets.point_coords(x=pX,y=pY,z=pZ)
			name = combine_strings(str_add=f'obj_{indMagn+1}',str_add_to=nameCore)
			magnetC._magnParas._magnetization=magnetizations[(startNum+indMagn)%len(magnetizations)]
			magnetC._magnParas._magn_unit_vec=undu_magnets.create_magnetization_unit_vec(
				magn_string=magnSeq[(startNum+indMagn)%len(magnSeq)],
				)
			magnetC._magnParas._len_x_main=magnetL*fact
			magnet = magnetC.createMagnet(center=magnCenter,name=name)
			magnet_blocks.append(magnet)
			pX=pX+magnetL/2.0+self._glueSlit
			newStart=newStart+1

		return undu_magnets.undu_magnets(magnet_blocks=magnet_blocks)

	def createRow(self,
			center, 
			magnetizations, 
			magnSeq, 
			nameCore='',
			pos='ll',
			) :

		nperiods=self.createNPeriods(
				center=center, 
				magnetizations=magnetizations, 
				magnSeq=magnSeq, 
				nameCore=nameCore,
				)

		endMag2=self._endMagnets[1%len(self._endMagnets)]
		endMag3=self._endMagnets[2%len(self._endMagnets)]
		endMagnLengthL=endMag2._magnParas._len_x_main*0.5+\
					endMag3._magnParas._len_x_main
		endStructLCenter=copy.deepcopy(center)
		endStructLCenter._x=center._x-\
				self._nperiods*self.getPeriodLength()/2.0-\
				self._keeperSlit-endMagnLengthL-self._glueSlit
		reversedEndMagn=list(reversed(self._endMagnets))
		reversedEndMagn=self._endMagnets
		endMagR1=reversedEndMagn[0%len(self._endMagnets)]
		endMagR2=reversedEndMagn[1%len(self._endMagnets)]
		endMagnLengthR=endMagR1._magnParas._len_x_main+\
					endMagR2._magnParas._len_x_main
		endStructRCenter=copy.deepcopy(center)
		endStructRCenter._x=center._x+\
				self._nperiods*self.getPeriodLength()/2.0+\
				self._keeperSlit+endMagnLengthR+self._glueSlit*1.5
		name = combine_strings(str_add=f'endL',str_add_to=nameCore)
		endStructL=self.createEndStructElliptStd(
			center=endStructLCenter, 
			location='us', 
			magnets=self._endMagnets,
			magnetizations=magnetizations, 
			magnSeq=magnSeq, 
			nameCore=name,
			startNum=1,
			)
		name = combine_strings(str_add=f'endR',str_add_to=nameCore)
		endStructR=self.createEndStructElliptStd(
			center=endStructRCenter, 
			location='ds', 
			magnets=reversedEndMagn,
			magnetizations=magnetizations, 
			magnSeq=magnSeq, 
			nameCore=name,
			startNum=0,
			)
		row1=undu_magnets.undu_magnets(magnet_blocks=[endStructL,nperiods,endStructR])

		axis_vec = undu_magnets.point_coords(0,0,0)
		# if self._unduType == 'planar' :
		# 	if (pos == 'ur') or (pos == 'ul') or (pos == 'upper') :
		# 		row1.mirror(coord='y')
		# elif self._unduType=='ellipt' :
		if pos == 'lr' :
			row1.mirror(coord='z')
		elif pos == 'ur' :
			row1.rotate(
				degrees=180,
				axis=axis_vec,
				)
		elif pos == 'ul' :
			row1.mirror(coord='y')

		return undu_magnets.undu_magnets(magnet_blocks=[endStructL,nperiods,endStructR])	

	def reflectMagnetizations(self,magnSeq,mirrorPlane='xz') :
		newMagnSeq=[]
		if mirrorPlane == 'xz' :
			for magn in magnSeq :
				if magn=='+x' :
					newMagnSeq.append('-x')
				elif magn=='-x' :
					newMagnSeq.append('+x')
				else:
					newMagnSeq.append(magn)					
		return newMagnSeq

	def createUndulator(self,
			center, 
			magnetizations, 
			magnSeq, 
			nameCore='',
			onlyLL=False,
			) :

		magsY=[]
		magsZ=[]
		for mag in self._magnetsPerKeeper :
			magsY.append(mag._magnParas._len_y_main)
			magsZ.append(mag._magnParas._len_z_main)
		maxY=max(magsY)
		maxZ=max(magsZ)
		center_ll=undu_magnets.point_coords(
			x=center._x,
			y=center._y-maxY/2.0,
			z=center._z-maxZ/2.0
			)
		shifts=[0.0,0.0,0.0,0.0]
		rowSlit=0.0
		if self._unduType=='ellipt' :
			rowSlit=self._rowSlit
			if self._shift >= 0 :
				shifts = [self._shift,0.0,0.0,self._shift]
			else :
				shifts = [-self._shift,0.0,0.0,self._shift]

		llRow=self.createRow(
			center=center_ll, 
			magnetizations=magnetizations, 
			magnSeq=magnSeq, 
			nameCore='ll',
			pos='ll'
			)
		v_move_ll = undu_magnets.point_coords(
			x=shifts[2],
			y=-self._gap/2.0,
			z=-rowSlit/2.0
			)
		llRow.move_it(vec=v_move_ll)
		if onlyLL or (self._unduType=='planar') :
			blocks = [llRow]
			return undu_magnets.undu_magnets(magnet_blocks=blocks)

		lrRow=self.createRow(
			center=center_ll, 
			magnetizations=magnetizations, 
			magnSeq=magnSeq, 
			nameCore='lr',
			pos='lr'
			)
		magnSeqMir=self.reflectMagnetizations(magnSeq=magnSeq,mirrorPlane='xz')
		ulRow=self.createRow(
			center=center_ll, 
			magnetizations=magnetizations, 
			magnSeq=magnSeqMir, 
			nameCore='ul',
			pos='ul'
			)

		urRow=self.createRow(
			center=center_ll, 
			magnetizations=magnetizations, 
			magnSeq=magnSeqMir, 
			nameCore='ur',
			pos='ur'
			)
		v_move_ul = undu_magnets.point_coords(
			x=shifts[0],
			y=self._gap/2.0,
			z=-rowSlit/2.0
			)

		v_move_ur = undu_magnets.point_coords(
			x=shifts[1],
			y=self._gap/2.0,
			z=rowSlit/2.0
			)

		v_move_lr = undu_magnets.point_coords(
			x=shifts[3],
			y=-self._gap/2.0,
			z=rowSlit/2.0
			)
		ulRow.move_it(vec=v_move_ul)
		urRow.move_it(vec=v_move_ur)
		lrRow.move_it(vec=v_move_lr)
		blocks = [ulRow,urRow,llRow,lrRow]
		# blocks = [llRow]

		return undu_magnets.undu_magnets(magnet_blocks=blocks)

	# def create_easy_endStruct(self,
	# 		center, 
	# 		location, 
	# 		magnParas, 
	# 		magnetizations, 
	# 		magnSeq, 
	# 		nameCore='',
	# 		) :
	# 	"""
	# 	creates the endstructure
	# 	type='us' or 'ds'
	# 	"""
	# 	blocks=[]

	# 	mlen=magnParas._len_x_main
	# 	if type=='us' :
	# 		name1 = combine_strings(str_add=f'end_us_obj_1',str_add_to=nameCore)
	# 		name2 = combine_strings(str_add=f'end_us_obj_2',str_add_to=nameCore)
	# 		name3 = combine_strings(str_add=f'end_us_obj_3',str_add_to=nameCore)

	# 		p_center_1=copy.deepcopy(p_center)
	# 		p_center_2=copy.deepcopy(p_center)
	# 		p_center_3=copy.deepcopy(p_center)
	# 		p_center_1._x=p_center_1._x-mlen/2.0-mlen/4.0-mlen/8.0
	# 		p_center_3._x=p_center_3._x+mlen/4.0+mlen/2.0+mlen*3.0/8.0
	# 		endMagParas1 = get_end_magn_paras(endNum=3,magnParas=magnParas)
	# 		endMagParas2 = get_end_magn_paras(endNum=2,magnParas=magnParas)
	# 		endMagParas3 = get_end_magn_paras(endNum=1,magnParas=magnParas)
	# 	elif type=='ds' :
	# 		name1 = combine_strings(str_add=f'end_ds_obj_1',str_add_to=nameCore)
	# 		name2 = combine_strings(str_add=f'end_ds_obj_2',str_add_to=nameCore)
	# 		name3 = combine_strings(str_add=f'end_ds_obj_3',str_add_to=nameCore)

	# 		p_center_1=copy.deepcopy(p_center)
	# 		p_center_2=copy.deepcopy(p_center)
	# 		p_center_3=copy.deepcopy(p_center)
	# 		p_center_1._x=p_center_1._x-mlen/4.0-mlen/2.0-mlen*3.0/8.0
	# 		p_center_3._x=p_center_3._x+mlen/2.0+mlen/4.0+mlen/8.0
	# 		endMagParas1 = get_end_magn_paras(endNum=1,magnParas=magnParas)
	# 		endMagParas2 = get_end_magn_paras(endNum=2,magnParas=magnParas)
	# 		endMagParas3 = get_end_magn_paras(endNum=3,magnParas=magnParas)

	# 	endMagParas1._magnetization = magnetizations[0%len(magnetizations)]
	# 	endMagParas1._magn_unit_vec = undu_magnets.create_magnetization_unit_vec(magn_string=magnSeq[0%len(magnSeq)])
	# 	endMagParas2._magnetization = magnetizations[1%len(magnetizations)]
	# 	endMagParas2._magn_unit_vec = undu_magnets.create_magnetization_unit_vec(magn_string=magnSeq[1%len(magnSeq)])
	# 	endMagParas3._magnetization = magnetizations[2%len(magnetizations)]
	# 	endMagParas3._magn_unit_vec = undu_magnets.create_magnetization_unit_vec(magn_string=magnSeq[2%len(magnSeq)])

	# 	magnet1 = create_easy_magnet( 
	# 		p_center=p_center_1, 
	# 		magn_para=endMagParas1,
	# 		name=name1
	# 		)
	# 	magnet2 = create_easy_magnet( 
	# 		p_center=p_center_2, 
	# 		magn_para=endMagParas2,
	# 		name=name2
	# 		)
	# 	magnet3 = create_easy_magnet( 
	# 		p_center=p_center_3, 
	# 		magn_para=endMagParas3,
	# 		name=name3
	# 		)
	# 	blocks = blocks + magnet1._magnet_blocks + magnet2._magnet_blocks + magnet3._magnet_blocks

	# 	return undu_magnets.undu_magnets(magnet_blocks=blocks)
