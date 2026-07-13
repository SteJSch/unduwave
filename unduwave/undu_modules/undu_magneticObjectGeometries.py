"""

"""

from unduwave.unduwave_incl import *
from unduwave.undu_modules.undu_parameters import *
import unduwave.helpers.file_folder_helpers as f_h
from unduwave.quantities.quantities import *

import unduwave.undu_modules.undu_blocks as undu_blocks

def create_cpmuStdPole_geometry(
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

def create_cpmuStdMagnet_geometry(
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

def create_twoSidedClamps_a( 
		center, 
		magn_paras,
		magn_paras_clamp,
		name='mag'
		) : 
	"""
	Two-sided Clamps - a
	    _____
       |     |
	  _|     |
	 |    ___|
	 |___|
	"""

	magn_paras_main=copy.deepcopy(magn_paras)

	name_main = combine_strings(str_add='main',str_add_to=name)

	heightAll=magn_paras._len_y_main()
	widthAll=magn_paras._len_z_main()

	heightMain=heightAll-magn_paras_clamp._len_z_main()
	widthMain=widthAll-magn_paras_clamp._len_z_main()

	heightS1=heightAll-magn_paras_clamp._len_z_main()-magn_paras_clamp._len_y_main()
	widthS1=magn_paras_clamp._len_z_main()

	heightS2=magn_paras_clamp._len_z_main()
	widthS2=widthAll-magn_paras_clamp._len_y_main()

	magn_paras_main._len_z_main.set(widthMain)
	magn_paras_main._len_y_main.set(heightMain)

	p_main = copy.deepcopy(center)
	p_main[1]=p_main[1]-0.5*heightAll+heightS2\
					+0.5*heightMain
	p_main[2]=p_main[2]-0.5*widthAll+widthS1+0.5*widthMain

	main_block = undu_blocks.undumagBlockObject(
		center=p_main,
		magnParas=magn_paras_main,
		pnts=None,
		name=name_main,
		parentName=name,
		api=None
		)

	# 1st one lies vertical, to the lower left of main block
	name_side = combine_strings(str_add='side',str_add_to=name)

	magn_paras_side1=undu_blocks.magParameters(
		material_id=magn_paras_main._material_id(),
		len_x_main=magn_paras_clamp._len_x_main(), 
		len_y_main=heightS1, 
		len_z_main=widthS1, 
		magnetization=magn_paras_main._magnetization(),
		magn_unit_vec=magn_paras_main._magn_unit_vec(),
		segm_x=magn_paras_clamp._segm_x(),
		segm_y=magn_paras_clamp._segm_y(),
		segm_z=magn_paras_clamp._segm_z(),
		chamf=magn_paras_clamp._chamf(),
		frac_y=magn_paras_clamp._frac_y(),
		frac_z=magn_paras_clamp._frac_z(),
		)

	p_side = copy.deepcopy(center)
	p_side[1] = p_side[1]-0.5*heightAll+ \
			heightS2+0.5*heightS1
	p_side[2] = p_side[2]-0.5*widthAll+0.5*widthS1

	side_block = undu_blocks.undumagBlockObject(
		center=p_side,
		magnParas=magn_paras_side1,
		pnts=None,
		name=name_side,
		parentName=name,
		api=None
		)

	# 2nd one lies horizontally, to the lower left of main block
	name_side2 = combine_strings(str_add='side2',str_add_to=name)

	magn_paras_side2=undu_blocks.magParameters(
		material_id=magn_paras_main._material_id(),
		len_x_main=magn_paras_clamp._len_x_main(), 
		len_y_main=heightS2,
		len_z_main=widthS2,
		magnetization=magn_paras_main._magnetization(),
		magn_unit_vec=magn_paras_main._magn_unit_vec(),
		segm_x=magn_paras_clamp._segm_x(),
		segm_y=magn_paras_clamp._segm_z(),
		segm_z=magn_paras_clamp._segm_y(),
		chamf=magn_paras_clamp._chamf(),
		frac_y=magn_paras_clamp._frac_z(),
		frac_z=magn_paras_clamp._frac_y(),
		)

	p_side2 = copy.deepcopy(center)
	p_side2[1] = p_side2[1]-0.5*heightAll+ \
			0.5*heightS2
	p_side2[2] = p_side2[2]-0.5*widthAll+0.5*widthS2

	side_block2 = undu_blocks.undumagBlockObject(
		center=p_side2,
		magnParas=magn_paras_side2,
		pnts=None,
		name=name_side2,
		parentName=name,
		api=None
		)

	# return undu_blocks.undumagObjectList(magnet_blocks=[main_block])
	return undu_blocks.undumagObjectList(magnet_blocks=[main_block,side_block,side_block2])

def create_twoSidedClamps_b( 
		center, 
		magn_paras,
		magn_paras_clamp,
		name='mag'
		) : 
	"""
	Two-sided Clamps - b
	    ____
	 __|    |__
	|__________|
	"""

	name_main = combine_strings(str_add='main',str_add_to=name)

	magn_paras_main=copy.deepcopy(magn_paras)
	magn_paras_main._len_y_main.set(magn_paras_main._len_y_main()-magn_paras_clamp._len_y_main())

	p_main = copy.deepcopy(center)
	p_main[1]=p_main[1]-0.5*magn_paras._len_y_main()+0.5*magn_paras_main._len_y_main()

	main_block = undu_blocks.undumagBlockObject(
		center=p_main,
		magnParas=magn_paras_main,
		pnts=None,
		name=name_main,
		parentName=name,
		api=None
		)

	name_side = combine_strings(str_add='side',str_add_to=name)

	magn_paras_side=undu_blocks.magParameters(
		material_id=magn_paras_main._material_id(),
		len_x_main=magn_paras_clamp._len_x_main(), 
		len_y_main=magn_paras._len_y_main()-magn_paras_main._len_y_main(), 
		len_z_main=magn_paras_main._len_z_main()-2*magn_paras_clamp._len_z_main(), 
		magnetization=magn_paras_main._magnetization(),
		magn_unit_vec=magn_paras_main._magn_unit_vec(),
		segm_x=magn_paras_clamp._segm_x(),
		segm_y=magn_paras_clamp._segm_y(),
		segm_z=magn_paras_clamp._segm_z(),
		chamf=magn_paras_clamp._chamf(),
		frac_y=magn_paras_clamp._frac_y(),
		frac_z=magn_paras_clamp._frac_z(),
		)

	p_side = copy.deepcopy(center)
	p_side[1] = p_side[1]+0.5*magn_paras._len_y_main()-0.5*magn_paras_clamp._len_y_main()

	side_block = undu_blocks.undumagBlockObject(
		center=p_side,
		magnParas=magn_paras_side,
		pnts=None,
		name=name_side,
		parentName=name,
		api=None
		)

	return undu_blocks.undumagObjectList(magnet_blocks=[main_block,side_block])

	# p_main = copy.deepcopy(p_center)
	# p_main._y = p_main._y-0.5*(corr_magn_paras._len_y_side)

	# name_main = ana.combine_strings(str_add='main',str_add_to=name)
	# main_block = undu_magnets.undu_magnet_block_coords(p_center=p_main,len_x=corr_magn_paras._len_x_main,len_y=corr_magn_paras._len_y_main,
	# 	len_z=corr_magn_paras._len_z_main,magnetization=corr_magn_paras._magnetization,magn_unit_vec=corr_magn_paras._magn_unit_vec,
	# 	name=name_main,mother=name,
	# 	segm_x = corr_magn_paras._segm_x, segm_y = corr_magn_paras._segm_y, 
	# 	segm_z = corr_magn_paras._segm_z, material = "mag")

	# x_block1 = p_main._x
	# y_block1 = (p_main._y+corr_magn_paras._len_y_main/2.0)+corr_magn_paras._len_y_side/2.0
	# z_block1 = p_main._z
	# p_center_block1 = undu_magnets.point_coords(x=x_block1,y=y_block1,z=z_block1)

	# name_s = ana.combine_strings(str_add='side',str_add_to=name)
	# side_block_1 = undu_magnets.undu_magnet_block_coords(p_center=p_center_block1,
	# 	len_x=corr_magn_paras._len_x_side,len_y=corr_magn_paras._len_y_side,
	# 	len_z=corr_magn_paras._len_z_side,magnetization=corr_magn_paras._magnetization,magn_unit_vec=corr_magn_paras._magn_unit_vec,
	# 	name=name_s,mother=name,
	# 	segm_x = corr_magn_paras._segm_x, segm_y = corr_magn_paras._segm_y, 
	# 	segm_z = corr_magn_paras._segm_z, material = "mag")

	# return undu_magnets.undu_magnets(magnet_blocks=[main_block,side_block_1])
