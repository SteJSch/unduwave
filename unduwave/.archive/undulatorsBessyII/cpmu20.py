"""
Contains the wave_from_b class that incorporates the API from python to WAVE and the function create_wave_instance that 
returns an instance of that class
"""

# import unduwave as unduwave
from unduwave.unduwave_incl import *
from unduwave.attribute_classes.attributes import _attribute
from unduwave.attribute_classes.attributes import _attribute_collection
import unduwave.helpers.file_folder_helpers as f_h
from unduwave import undu_representation
import unduwave.constants as uc
from unduwave import bfield
from unduwave import undu_blocks

class cpmuParas(_attribute_collection) : 
	def __init__(self,
			file=None,
			center=None,
			magn_seqs_ul=None,
			periodicMagnParas=None,
			periodicPolParas=None,
			endMagn1Paras=None,
			endMagn2Paras=None,
			endPolParas=None,
			gap = 0.0, 
			nperiods=1,
			magnetization=None,
			polDrop=0.0,
			magDrop=0.0,
			end_pole_drop = 0.0,
			mag_pole_slit = 0.0, 
			pole_mag_slit = 0.0, 
			endmag2_pole_l_slit=0.0,
			pole_endmag1_l_slit=0.0,
			endmag1_pole_l_slit=0.0,
			pole_endmag1_r_slit=0.0,
			endmag1_pole_r_slit=0.0,
			pole_endmag2_r_slit=0.0,
			endMagn1Drop=0.0,
			endMagn2Drop=0.0,
			coating=0.003,
			coatingMaterial='TiN',
			):
		if not (file is None) :
			obj=self.fromFile(file=file)
			self.__dict__=obj.__dict__
			super().__init__()
			return
		self._periodicMagnParas=periodicMagnParas
		self._periodicPolParas=periodicPolParas
		self._endMagn1Paras=endMagn1Paras
		self._endMagn2Paras=endMagn2Paras
		self._endPolParas=endPolParas
		self._gap=_attribute(gap)
		self._nperiods=_attribute(nperiods)
		self._magnetization=_attribute(magnetization)
		self._end_pole_drop=_attribute(end_pole_drop)
		self._mag_pole_slit=_attribute(mag_pole_slit)
		self._pole_mag_slit=_attribute(pole_mag_slit)
		self._endmag2_pole_l_slit=_attribute(endmag2_pole_l_slit)
		self._pole_endmag1_l_slit=_attribute(pole_endmag1_l_slit)
		self._endmag1_pole_l_slit=_attribute(endmag1_pole_l_slit)
		self._pole_endmag1_r_slit=_attribute(pole_endmag1_r_slit)
		self._endmag1_pole_r_slit=_attribute(endmag1_pole_r_slit)
		self._pole_endmag2_r_slit=_attribute(pole_endmag2_r_slit)
		self._magn_seqs_ul=_attribute(magn_seqs_ul)
		self._polDrop=_attribute(polDrop)
		self._magDrop=_attribute(magDrop)
		self._endMagn1Drop=_attribute(endMagn1Drop)
		self._endMagn2Drop=_attribute(endMagn2Drop)
		self._coating=coating
		self._coatingMaterial=coatingMaterial
		if center is None:
			center=np.array([0.0,0.0,0.0])
		self._center=_attribute(center)
		super().__init__()

	def writeValsJSON(self,file) :

		myTable={
			'gap' : self._gap(),
			'nperiods' : self._nperiods(),
			'lMag' : self._periodicMagnParas._len_x_main(),
			'end magnet 1 drop' : self._endMagn1Drop(),
			'end magnet 2 drop' : self._endMagn2Drop(),
			'chamfers' : self._periodicMagnParas._chamf(),
			'lPolMZ' : self._periodicPolParas._len_z_main(),
			'lPolLSY' : self._periodicPolParas._len_y_side_low(),
			'lPolMSY' : self._periodicPolParas._len_y_side_middle(),
			'lPolSZ' : self._periodicPolParas._len_z_side(),
			'lPolMChamf' : 0.5,
			'lPolMY' : self._periodicPolParas._len_y_main(),
			'lPolMZ_EP' : self._endPolParas._len_z_main(),
			'lPole_EP' : self._endPolParas._len_x_main(),
			'lPolMY_EP' : self._endPolParas._len_y_main(),
			'lMagMY' : self._periodicMagnParas._len_y_main(),
			'lMagMZ' : self._periodicMagnParas._len_z_main(),
			'lMagSY' : self._periodicMagnParas._len_y_side(),
			'lMagSZ' : self._periodicMagnParas._len_z_side(),
			'lMagE1MY' : self._endMagn1Paras._len_y_main(),
			'lMagE1MZ' : self._endMagn1Paras._len_z_main(),
			'lMagE1SY' : self._endMagn1Paras._len_y_side(),
			'lMagE1SZ' : self._endMagn1Paras._len_z_side(),
			'lMagE2MY' : self._endMagn2Paras._len_y_main(),
			'lMagE2MZ' : self._endMagn2Paras._len_z_main(),
			'lMagE1' : self._endMagn1Paras._len_x_main(),
			'lMagE2' : self._endMagn2Paras._len_x_main(),
			'end pole drop' : self._end_pole_drop(),
			'mag_pole_slit' : self._mag_pole_slit(),
			'pole_mag_slit' : self._pole_mag_slit(),
			'endmag2_pole_l_slit' : self._endmag2_pole_l_slit(),
			'pole_endmag1_l_slit' : self._pole_endmag1_l_slit(),
			'endmag1_pole_l_slit' : self._endmag1_pole_l_slit(),
			'pole_endmag1_r_slit' : self._pole_endmag1_r_slit(),
			'endmag1_pole_r_slit' : self._endmag1_pole_r_slit(),
			'pole_endmag2_r_slit' : self._pole_endmag2_r_slit(),
		}
		for key, el in myTable.items() :
			myTable[key]=float(el)

		files = [(file,  {"indent": 4}) ]
		for filename, options in files:
			with open(filename, "w") as f:
				try:
					json.dump(myTable, f, **options)
				except:
					pdb.set_trace()
			return myTable


def createPole(
		center, # coordinates such that center is upper right corner of main magnet (for lower left row)
		poleParasMain,
		poleParasLowerSide,
		poleParasUpperSide,
		cliffSide,
		name='pol'
		) : 

	pnts_side1 = [
		np.array([
			-poleParasMain._len_x_main()/2.0,
			poleParasUpperSide._len_y_main()/2.0-cliffSide,
			-poleParasUpperSide._len_z_main()/2.0
			]),
		np.array([
			-poleParasMain._len_x_main()/2.0,
			poleParasUpperSide._len_y_main()/2.0,
			-poleParasUpperSide._len_z_main()/2.0+cliffSide,
			]),
		np.array([
			-poleParasMain._len_x_main()/2.0,
			poleParasUpperSide._len_y_main()/2.0,
			+poleParasUpperSide._len_z_main()/2.0
			]),
		np.array([
			-poleParasMain._len_x_main()/2.0,
			-poleParasUpperSide._len_y_main()/2.0,
			+poleParasUpperSide._len_z_main()/2.0
			]),
		np.array([
			-poleParasMain._len_x_main()/2.0,
			-poleParasUpperSide._len_y_main()/2.0,
			-poleParasUpperSide._len_z_main()/2.0
			]),
		np.array([
			poleParasMain._len_x_main()/2.0,
			poleParasUpperSide._len_y_main()/2.0-cliffSide,
			-poleParasUpperSide._len_z_main()/2.0
			]),
		np.array([
			poleParasMain._len_x_main()/2.0,
			poleParasUpperSide._len_y_main()/2.0,
			-poleParasUpperSide._len_z_main()/2.0+cliffSide
			]),
		np.array([
			poleParasMain._len_x_main()/2.0,
			poleParasUpperSide._len_y_main()/2.0,
			+poleParasUpperSide._len_z_main()/2.0
			]),
		np.array([
			poleParasMain._len_x_main()/2.0,
			-poleParasUpperSide._len_y_main()/2.0,
			+poleParasUpperSide._len_z_main()/2.0
			]),
		np.array([
			poleParasMain._len_x_main()/2.0,
			-poleParasUpperSide._len_y_main()/2.0,
			-poleParasUpperSide._len_z_main()/2.0
			]),
	]

	pnts_side2 = [
		np.array([
			-poleParasMain._len_x_main()/2.0,
			-poleParasLowerSide._len_y_main()/2.0,
			poleParasLowerSide._len_z_main()/2.0
			]),
		np.array([
			0.0,
			-poleParasLowerSide._len_y_main()/2.0,
			-poleParasLowerSide._len_z_main()/2.0
			]),
		np.array([
			poleParasMain._len_x_main()/2.0,
			-poleParasLowerSide._len_y_main()/2.0,
			poleParasLowerSide._len_z_main()/2.0
			]),
		np.array([
			poleParasMain._len_x_main()/2.0,
			poleParasLowerSide._len_y_main()/2.0,
			poleParasLowerSide._len_z_main()/2.0
			]),
		np.array([
			0.0,
			poleParasLowerSide._len_y_main()/2.0,
			-poleParasLowerSide._len_z_main()/2.0
			]),
		np.array([
			-poleParasMain._len_x_main()/2.0,
			poleParasLowerSide._len_y_main()/2.0,
			poleParasLowerSide._len_z_main()/2.0
			]),
	]

	name_main = combine_strings(str_add='main',str_add_to=name)
	p_main = copy.deepcopy(center)

	# p_main[1]=p_main[1]-poleParasMain._len_y_main()/2.0-poleDrop
	# p_main[2] = p_main[2] - poleParasMain._len_z_main()/2.0
	main_block = undu_blocks.undumagBlockObject(
		center=p_main,
		magnParas=poleParasMain,
		pnts=None,
		name=name_main,
		parentName=name,
		api=None
		)

	name_u_side_left = combine_strings(str_add='uside_l',str_add_to=name)
	p_usl = copy.deepcopy(center)
	p_usl[1] = p_main[1] - poleParasMain._len_y_main()/2.0 + poleParasLowerSide._len_y_main() + poleParasUpperSide._len_y_main()/2.0
	p_usl[2] = p_main[2] - poleParasMain._len_z_main()/2.0 - poleParasLowerSide._len_z_main()/2.0
	upper_side_block_left = undu_blocks.undumagBlockObject(
		center=p_usl,
		magnParas=poleParasUpperSide,
		pnts=copy.deepcopy(pnts_side1),
		name=name_u_side_left,
		parentName=name,
		api=None
		)

	name_l_side_left = combine_strings(str_add='lside_l',str_add_to=name)
	p_lsl = copy.deepcopy(center)
	p_lsl[1] = p_main[1] - poleParasMain._len_y_main()/2.0 + poleParasLowerSide._len_y_main()/2.0
	p_lsl[2] = p_main[2] - poleParasMain._len_z_main()/2.0 - poleParasLowerSide._len_z_main()/2.0
	lower_side_block_left = undu_blocks.undumagBlockObject(
		center=p_lsl,
		magnParas=poleParasLowerSide,
		pnts=copy.deepcopy(pnts_side2),
		name=name_l_side_left,
		parentName=name,
		api=None
		)

	for ind,pnt in enumerate(pnts_side1):
		pnts_side1[ind][2] = - pnt[2]

	name_u_side_right = combine_strings(str_add='uside_r',str_add_to=name)
	p_usr = copy.deepcopy(center)
	p_usr[1] = p_main[1] - poleParasMain._len_y_main()/2.0 + poleParasLowerSide._len_y_main() + poleParasUpperSide._len_y_main()/2.0
	p_usr[2] = p_main[2] + poleParasMain._len_z_main()/2.0 + poleParasUpperSide._len_z_main()/2.0
	upper_side_block_right = undu_blocks.undumagBlockObject(
		center=p_usr,
		magnParas=copy.deepcopy(poleParasUpperSide),
		pnts=copy.deepcopy(pnts_side1),
		name=name_u_side_right,
		parentName=name,
		api=None
		)

	for ind,pnt in enumerate(pnts_side2):
		pnts_side2[ind][2] = - pnt[2]

	name_l_side_right = combine_strings(str_add='lside_r',str_add_to=name)
	p_lsr = copy.deepcopy(center)
	p_lsr[1] = p_main[1] - poleParasMain._len_y_main()/2.0 + poleParasLowerSide._len_y_main()/2.0
	p_lsr[2] = p_main[2] + poleParasMain._len_z_main()/2.0 + poleParasLowerSide._len_z_main()/2.0
	lower_side_block_right = undu_blocks.undumagBlockObject(
		center=p_lsr,
		magnParas=copy.deepcopy(poleParasLowerSide),
		pnts=copy.deepcopy(pnts_side2),
		name=name_l_side_right,
		parentName=name,
		api=None
		)

	# return undu_blocks.undumagObjectList(
	# 	magnet_blocks=[main_block]
	# 	)
	return undu_blocks.undumagObjectList(
		magnet_blocks=[main_block,upper_side_block_left,lower_side_block_left]
		)
	# return undu_blocks.undumagObjectList(
	# 	magnet_blocks=[main_block,upper_side_block_left,lower_side_block_left,upper_side_block_right,lower_side_block_right]
	# 	)

def createMagnet(
		center,
		magnParasMain,
		magnParasSide,
		name
		) :

	name_main = combine_strings(str_add='main',str_add_to=name)
	p_main = copy.deepcopy(center)

	# p_main[2] = p_main[2] - magnParasMain._len_z_main()/2.0
	main_block = undu_blocks.undumagBlockObject(
		center=p_main,
		magnParas=magnParasMain,
		pnts=None,
		name=name_main,
		parentName=name,
		api=None
		)

	name_s1 = combine_strings(str_add='sideL',str_add_to=name)
	p_side_left = copy.deepcopy(center)
	p_side_left[1] = p_side_left[1] - magnParasMain._len_y_main()/2.0 + \
		magnParasSide._len_y_main()/2.0
	p_side_left[2] = p_side_left[2]-magnParasMain._len_z_main()/2.0-magnParasSide._len_z_main()/2.0

	side_left = undu_blocks.undumagBlockObject(
		center=p_side_left,
		magnParas=magnParasSide,
		pnts=None,
		name=name_s1,
		parentName=name,
		api=None
		)
	name_s2 = combine_strings(str_add='sideR',str_add_to=name)
	p_side_right = copy.deepcopy(center)
	p_side_right[1] = p_side_right[1] - magnParasMain._len_y_main()/2.0 + \
		magnParasSide._len_y_main()/2.0
	p_side_right[2] = p_side_right[2]+magnParasMain._len_z_main()/2.0+magnParasSide._len_z_main()/2.0
	side_right = undu_blocks.undumagBlockObject(
		center=p_side_right,
		magnParas=magnParasSide,
		pnts=None,
		name=name_s2,
		parentName=name,
		api=None
		)

	# return undu_blocks.undumagObjectListmagnet_blocks=[main_block])
	return undu_blocks.undumagObjectList(magnet_blocks=[main_block,side_left])
	# return undu_blocks.undumagObjectList(magnet_blocks=[main_block,side_left,side_right])

def getCPMU20Paras_OldOptm() :
	magnParas=undu_blocks.magParameters(
		len_x_main=6.1, 
		len_y_main=37,
		len_z_main=49/2.0, 
		segm_x=4,
		segm_y=2,
		segm_z=2,
		frac_y=1,
		frac_z=1,
		chamf=0.3,
		material="mag", #"mag","pole", 'NiCuFoil'
	)
	magnParas._len_z_side=_attribute(5)
	magnParas._len_y_side=_attribute(22.7)
	endMagn1Paras=undu_blocks.magParameters(
		len_x_main=6.3, 
		len_y_main=32.8,
		len_z_main=37.8/2.0, 
		segm_x=4,
		segm_y=2,
		segm_z=2,
		frac_y=1,
		frac_z=1,
		chamf=0.3,
		material="mag", #"mag","pole", 'NiCuFoil'
	)
	endMagn1Paras._len_z_side=_attribute(5)
	endMagn1Paras._len_y_side=_attribute(19.9)
	endMagn2Paras=undu_blocks.magParameters(
		len_x_main=6.5, 
		len_y_main=12.8,
		len_z_main=48.4/2.0, 
		segm_x=4,
		segm_y=2,
		segm_z=2,
		frac_y=1,
		frac_z=1,
		chamf=0.3,
		material="mag", #"mag","pole", 'NiCuFoil'
	)
	polParas=undu_blocks.magParameters(
		len_x_main=3.4, 
		len_y_main=28.5,
		len_z_main=25.2/2.0, 
		segm_x=3,
		segm_y=12,
		segm_z=3,
		frac_y=12,
		frac_z=3,
		chamf=0.3,
		material="pole", #"mag","pol", 'NiCuFoil'
	)
	polParas._len_z_side=_attribute(3)
	polParas._len_y_side_low=_attribute(21.9)
	polParas._len_y_side_middle=_attribute(5.4)
	polParas._cliff_side=_attribute(0.5)
	endPolParas=copy.deepcopy(polParas)
	endPolParas._len_z_main.set(34.6/2.0)
	endPolParas._len_y_main.set(35.8)
	cpmuPs=cpmuParas(
		file=None,
		center=None,
		magn_seqs_ul=[],
		periodicMagnParas=magnParas,
		periodicPolParas=polParas,
		endMagn1Paras=endMagn1Paras,
		endMagn2Paras=endMagn2Paras,
		endPolParas=endPolParas,
		gap = 6.0, # 5.5 
		nperiods=12, #73
		magnetization=1.62, # 1.58
		polDrop=0.0,
		magDrop=-0.07,
		end_pole_drop = 0.4,
		mag_pole_slit = 0.0, #0.0
		pole_mag_slit = 0.5, #0.5 
		endmag2_pole_l_slit=0.1,
		pole_endmag1_l_slit=0.7,
		endmag1_pole_l_slit=0.7,
		pole_endmag1_r_slit=0.7,
		endmag1_pole_r_slit=0.1,
		pole_endmag2_r_slit=0.1,
		endMagn1Drop=9.2,
		endMagn2Drop=14.1,
		)
	return cpmuPs

def getCPMU17Paras() :
	magnParas=undu_blocks.magParameters(
		len_x_main=5.3, 
		len_y_main=35,
		len_z_main=40/2.0, 
		segm_x=4,
		segm_y=2,
		segm_z=2,
		frac_y=1,
		frac_z=1,
		chamf=0.3,
		material="mag", #"mag","pole", 'NiCuFoil'
	)
	magnParas._len_z_side=_attribute(5)
	magnParas._len_y_side=_attribute(22.7)
	endMagn1Paras=undu_blocks.magParameters(
		len_x_main=5.3, 
		len_y_main=32.2,
		len_z_main=40/2.0, 
		segm_x=4,
		segm_y=2,
		segm_z=2,
		frac_y=1,
		frac_z=1,
		chamf=0.3,
		material="mag", #"mag","pole", 'NiCuFoil'
	)
	endMagn1Paras._len_z_side=_attribute(5)
	endMagn1Paras._len_y_side=_attribute(19.9)
	endMagn2Paras=undu_blocks.magParameters(
		len_x_main=5.3, 
		len_y_main=12.5,
		len_z_main=50/2.0, 
		segm_x=4,
		segm_y=2,
		segm_z=2,
		frac_y=1,
		frac_z=1,
		chamf=0.3,
		material="mag", #"mag","pole", 'NiCuFoil'
	)
	polParas=undu_blocks.magParameters(
		len_x_main=2.7, 
		len_y_main=30,
		len_z_main=34/2.0, 
		segm_x=3,
		segm_y=12,
		segm_z=3,
		frac_y=12,
		frac_z=3,
		chamf=0.3,
		material="pole", #"mag","pole", 'NiCuFoil'
	)
	polParas._len_z_side=_attribute(3)
	polParas._len_y_side_low=_attribute(18)
	polParas._len_y_side_middle=_attribute(2.84)
	polParas._cliff_side=_attribute(0.5)
	endPolParas=copy.deepcopy(polParas)
	cpmuPs=cpmuParas(
		file=None,
		center=None,
		magn_seqs_ul=[],
		periodicMagnParas=magnParas,
		periodicPolParas=polParas,
		endMagn1Paras=endMagn1Paras,
		endMagn2Paras=endMagn2Paras,
		endPolParas=endPolParas,
		gap = 5.5, 
		nperiods=1,
		magnetization=1.58,
		polDrop=0.0,
		magDrop=-0.07,
		end_pole_drop = -0.07,
		mag_pole_slit = 0.0, 
		pole_mag_slit = 0.5, 
		endmag2_pole_l_slit=0.0,
		pole_endmag1_l_slit=0.5,
		endmag1_pole_l_slit=0.0,
		pole_endmag1_r_slit=0.5,
		endmag1_pole_r_slit=0.0,
		pole_endmag2_r_slit=0.5,
		endMagn1Drop=1.81,
		endMagn2Drop=4.66,
		)
	return cpmuPs

def getCPMU20Paras() :
	width_z=40
	magnParas=undu_blocks.magParameters(
		len_x_main=6.1, 
		len_y_main=35,
		len_z_main=width_z/2.0, 
		segm_x=4,
		segm_y=2,
		segm_z=2,
		frac_y=1,
		frac_z=1,
		chamf=0.3,
		material="mag", #"mag","pole", 'NiCuFoil'
	)
	magnParas._len_z_side=_attribute(5)
	magnParas._len_y_side=_attribute(22.7)
	endMagn1Paras=undu_blocks.magParameters(
		len_x_main=7.9, 
		len_y_main=29.0,
		len_z_main=width_z/2.0, 
		segm_x=4,
		segm_y=2,
		segm_z=2,
		frac_y=1,
		frac_z=1,
		chamf=0.3,
		material="mag", #"mag","pole", 'NiCuFoil'
	)
	endMagn1Paras._len_z_side=_attribute(5)
	endMagn1Paras._len_y_side=_attribute(19.9)
	endMagn2Paras=undu_blocks.magParameters(
		len_x_main=6.5, 
		len_y_main=18.1,
		len_z_main=50/2.0, 
		segm_x=4,
		segm_y=2,
		segm_z=2,
		frac_y=1,
		frac_z=1,
		chamf=0.3,
		material="mag", #"mag","pole", 'NiCuFoil'
	)
	polParas=undu_blocks.magParameters(
		len_x_main=3.4, 
		len_y_main=30,
		len_z_main=34/2.0, 
		segm_x=4,
		segm_y=12,
		segm_z=3,
		frac_y=12,
		frac_z=3,
		chamf=0.3,
		material="pole", #"mag","pole", 'NiCuFoil'
	)
	polParas._len_z_side=_attribute(3)
	polParas._len_y_side_low=_attribute(18)
	polParas._len_y_side_middle=_attribute(2.84)
	polParas._cliff_side=_attribute(0.5)
	endPolParas=copy.deepcopy(polParas)
	endPolParas._len_x_main.set(2.2)
	endPolParas._len_z_main.set(34.0/2.0)
	endPolParas._len_y_main.set(35.0)
	cpmuPs=cpmuParas(
		file=None,
		center=None,
		magn_seqs_ul=[],
		periodicMagnParas=magnParas,
		periodicPolParas=polParas,
		endMagn1Paras=endMagn1Paras,
		endMagn2Paras=endMagn2Paras,
		endPolParas=endPolParas,
		gap = 5.5, # 5.5 
		nperiods=10, #73
		magnetization=1.58, # 1.58
		polDrop=0.0,
		magDrop=-0.07,
		end_pole_drop = 0.4,
		mag_pole_slit = 0.0, 
		pole_mag_slit = 0.5, 
		endmag2_pole_l_slit=0.0,
		pole_endmag1_l_slit=0.5,
		endmag1_pole_l_slit=0.0,
		pole_endmag1_r_slit=0.5,
		endmag1_pole_r_slit=0.0,
		pole_endmag2_r_slit=0.5,
		endMagn1Drop=9.2,
		endMagn2Drop=14.1,
		)
	return cpmuPs

def getCPMURepresentation(cpmuPs,undu,center) :

	magnetCenter=np.array([0.0,0.0,0.0])
	magnetCenter[1]=-cpmuPs._periodicMagnParas._len_y_main()/2.0-cpmuPs._gap()/2.0-cpmuPs._magDrop()
	myMagnetRepr=undu_representation.magnetRepresentation(
		center=magnetCenter,
		magnParas=cpmuPs._periodicMagnParas,
		magnFun=cpmuMagnetFun,
		)

	polCenter=np.array([0.0,0.0,0.0])
	polCenter[1]=-cpmuPs._periodicPolParas._len_y_main()/2.0-cpmuPs._gap()/2.0-cpmuPs._polDrop()
	myPolRepr=undu_representation.magnetRepresentation(
		center=polCenter,
		magnParas=cpmuPs._periodicPolParas,
		magnFun=cpmuPolFun,
		)

	endMagnet1Center=np.array([0.0,0.0,0.0])
	# endMagnet1Center[1]=-cpmuPs._endMagn1Paras._len_y_main()/2.0-cpmuPs._gap()/2.0-cpmuPs._endMagn1Drop()
	endMagnet1Center[1]=-cpmuPs._periodicMagnParas._len_y_main()+cpmuPs._endMagn1Paras._len_y_main()/2.0-cpmuPs._gap()/2.0
	endMagnet1Repr=undu_representation.magnetRepresentation(
		center=endMagnet1Center,
		magnParas=cpmuPs._endMagn1Paras,
		magnFun=cpmuMagnetFun,
		)
	endMagnet2Center=np.array([0.0,0.0,0.0])
	# endMagnet2Center[1]=-cpmuPs._endMagn2Paras._len_y_main()/2.0-cpmuPs._gap()/2.0-cpmuPs._endMagn2Drop()
	endMagnet2Center[1]=-cpmuPs._periodicMagnParas._len_y_main()+cpmuPs._endMagn2Paras._len_y_main()/2.0-cpmuPs._gap()/2.0
	endMagnet2Center[2]=-cpmuPs._endMagn2Paras._len_z_main()/2.0
	endMagnet2Repr=undu_representation.magnetRepresentation(
		center=endMagnet2Center,
		magnParas=cpmuPs._endMagn2Paras,
		magnFun=None,
		)
	endPolCenter=np.array([0.0,0.0,0.0])
	# endPolCenter[1]=-cpmuPs._endPolParas._len_y_main()/2.0-cpmuPs._gap()/2.0-cpmuPs._end_pole_drop()
	endPolCenter[1]=-cpmuPs._periodicMagnParas._len_y_main()+cpmuPs._endPolParas._len_y_main()/2.0-cpmuPs._gap()/2.0
	endPolRepr=undu_representation.magnetRepresentation(
		center=endPolCenter,
		magnParas=cpmuPs._endPolParas,
		magnFun=cpmuPolFun,
		)

	undulatorMagPol=undu_representation.unduRepresentation(
		gap = cpmuPs._gap(), 
		shift=0.0,
		glueSlit = cpmuPs._pole_mag_slit(), #0.0, 
		keeperSlit= cpmuPs._mag_pole_slit(), #0.0, 
		rowSlit=0.0, 
		nperiods = cpmuPs._nperiods(),
		unduType='planar', # "planar", 'ellipt' 
		magnetsPerKeeper=[myPolRepr,myMagnetRepr],
		endStructType='std', # - std just puts em there, ellipt applies shrinking
		endMagnets=[
			[endMagnet2Repr,endPolRepr,endMagnet1Repr],
			[myPolRepr,endMagnet1Repr,endPolRepr,endMagnet2Repr]
			],
		endSlits=[
			[
				cpmuPs._endmag2_pole_l_slit(),
				cpmuPs._pole_endmag1_l_slit(),
				cpmuPs._endmag1_pole_l_slit(),
			],
			[
				cpmuPs._mag_pole_slit(),
				cpmuPs._pole_endmag1_r_slit(),
				cpmuPs._endmag1_pole_r_slit(),
				cpmuPs._pole_endmag2_r_slit(),
			],
			],
		)
	llRow=undulatorMagPol.createRow(
			center=center, 
			magnetizations=[cpmuPs._magnetization], 
			magnSeq=['-x','-y','+x','+y'], 
			nameCore='ll',
			)
	magnetObjects=llRow

	# undulator=undulatorMagPol.createUndulator(
	# 	center=np.array([0.0,0.0,0.0]), 
	# 	magnetizations=[1.22,1.33],
	# 	magnSeq=['-x','-y','+x','+y'], 
	# 	nameCore='',
	# 	onlyLL=False,
	# 	)
	# undulator.add_to_clc(api=undu)

	return undulatorMagPol, magnetObjects 