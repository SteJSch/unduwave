
import pdb
import sys
import os
sys.path.insert(0, '../../../')

try :
	# works when calling script with python3 script_file
	dir_path = os.path.dirname(os.path.realpath(__file__))
except:
	# works when calling script with exec from python console
	dir_path = os.getcwd()

import unduwave as uw
from unduwave import undu_blocks
import numpy as np

res_folder='res'
res_folder_full=dir_path+f'/{res_folder}/'

undu = uw.undu(undu_mode='from_undu_magns')
undu_prog_paras = undu._prog_paras
undu_prog_paras.res_folder.set(res_folder+'/')
undu_prog_paras.plotGeometry.set(1)
undu_prog_paras.create_z_sym.set(0)
undu_prog_paras.bmap_nz.set(10)
undu_prog_paras.bmap_nx.set(10)
undu_prog_paras.center_magnet_struct.set(0)

pos_magn=np.array([ 0.0, -15.0, 0.0 ])
magn_paras=undu_blocks.magParameters(
	len_x_main=10, 
	len_y_main=20,
	len_z_main=30, 
	segm_x=2,
	segm_y=2,
	segm_z=2,
	frac_y=1,
	frac_z=1,
	chamf=0.3,
	material_id="pm_rec_77K",
)

magnet = undu_blocks.undumagBlockObject(
	center=pos_magn,
	magnParas=magn_paras,
	name='MyMagnet',
	parentName='',
	)

pos_pole=np.array([ 12.0, -15.0, 0.0 ])
pole_paras=undu_blocks.magParameters(
	len_x_main=5, 
	len_y_main=15,
	len_z_main=25, 
	segm_x=3,
	segm_y=9,
	segm_z=5,
	frac_y=9,
	frac_z=5,
	chamf=0.3,
	material_id="fm_vanadium_permendur",
)

pole = undu_blocks.undumagBlockObject(
	center=pos_pole,
	magnParas=pole_paras,
	name='MyPole',
	parentName='',
	)

objs=undu_blocks.undumagObjectList(
	magnet_blocks=[magnet,pole],
	name='magnet_pole',
	parentName='',
	center=np.array([0.0,0.0,0.0]),
	)

undu.set_magnet_objects(magn_objects=objs)
undu.run()

results = undu.get_results()

trajx = results.get_result(which='trajx')
by = results.get_result(which='by')
bz = results.get_result(which='bz')
intBy = results.get_result(which='intBy')
intBz = results.get_result(which='intBz')
intBy2 = results.get_result(which='intBy2')
intBz2 = results.get_result(which='intBz2')

zProfile = results.get_result(which='profz')
ByProfile = results.get_result(which='profBy')
BzProfile = results.get_result(which='profBz')

bmap = results.get_result(which='bmap')

nfig=0
by.plot_over(
	x_quant=trajx,
	nfig=nfig,
	nosave=True,
	clear=True,
	)
nfig=bz.plot_over(
	x_quant=trajx,
	nfig=nfig,
	file_name=f'byz.png',
	plot=False,
	title=f'B$_y$ and B$_z$',
	)

intBy.plot_over(
	x_quant=trajx,
	nfig=nfig,
	nosave=True,
	clear=True,
	)
nfig=intBz.plot_over(
	x_quant=trajx,
	nfig=nfig,
	file_name=f'intByz.png',
	plot=False,
	title=f'1. Integral of \n B$_y$ and B$_z$',
	)

intBy2.plot_over(
	x_quant=trajx,
	nfig=nfig,
	nosave=True,
	clear=True,
	)
nfig=intBz2.plot_over(
	x_quant=trajx,
	nfig=nfig,
	file_name=f'int2Byz.png',
	plot=False,
	title=f'2. Integral of \n B$_y$ and B$_z$',
	)

ByProfile.plot_over(
	x_quant=zProfile,
	nfig=nfig,
	nosave=True,
	clear=True,
	)
nfig=BzProfile.plot_over(
	x_quant=zProfile,
	nfig=nfig,
	file_name=f'profileBz.png',
	plot=False,
	title=f'Profile of B$_y$ and B$_z$',
	)

bmap.plot_fld_map(
	bWhat='By', # "Bx", "By" or "Bz"
	xPos=None,
	yPos=0.0,
	zPos=None,
	nfig=nfig,
	filename=res_folder_full+"bymap.png",
	title="By Map",
	)

pdb.set_trace()

