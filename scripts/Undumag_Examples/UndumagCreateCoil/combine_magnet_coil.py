
import pdb
import sys
import os
import math
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

undu = uw.undu()
undu_prog_paras = undu._prog_paras
undu_prog_paras.res_folder.set(res_folder+'/')
undu_prog_paras.plotGeometry.set(1)
undu_prog_paras.create_z_sym.set(0)
undu_prog_paras.bmap_nz.set(10)
undu_prog_paras.bmap_nx.set(10)

pos_magnet=np.array([ 0.0, -15.0, 0.0 ])
center_coil = np.array([0.0,-20.0,0.0])
normal_vec = np.array([1.0,0.0,0.0])

"""
standard orientation for interpretation of the variable names is in 
(0,1,0) direction, so the coils lying in the x-z plane. 

coil_len_x - extent of coil in x-direction
coil_thick - thickness of coil in x-z-plane
outer_z - outer length of coil in z-direction
inner_z - inner length of coil in z-direction
inner_r - the inner curvature of the coils at the edges, the higher, the more circle-like
height - the extent in y-direction
segm_v - segmentations
segm_h=10
segm_r=10
filling-the percentage of coil querschnitt filled by coils
wire_diameter - diameter of one wire
n_windings=int(height*coil_thick/(math.pi*wire_diameter**2))
"""

coil_len_x=200 
coil_thick=20
outer_z=90
inner_z=outer_z-2*coil_thick
inner_r=20
height=5
segm_v=10
segm_h=10
segm_r=10
filling=0.6
wire_diameter=1
n_windings=int(height*coil_thick/(math.pi*wire_diameter**2))

coil = uw.undu_coils.coil(
	coil_type='RectWindings', 
	current=1.0,
	center=center_coil,
	normal_vec=normal_vec, 
	rot_angle=0.0,
	length=coil_len_x,
	inner_z=inner_z,
	outer_z=outer_z,
	inner_radius=inner_r, 
	height=height,
	n_vert=segm_v,
	n_hor=segm_h,
	n_rad=segm_r,
	filling = filling, 
	n_windings = n_windings,
	)

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

magnObject = undu_blocks.undumagBlockObject(
	center=pos_magnet,
	magnParas=magn_paras,
	name='MyMagnet',
	parentName='',
	)

objs=undu_blocks.undumagObjectList(
	magnet_blocks=[coil,magnObject],
	name='magnet_coil',
	parentName='system',
	api=None,
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