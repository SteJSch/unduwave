
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

res_folder='res'
res_folder_full=dir_path+f'/{res_folder}/'

undu = uw.undu(undu_mode='from_undu_magns')
undu_prog_paras = undu._prog_paras
undu_prog_paras.res_folder.set(res_folder+'/')
undu_prog_paras.plotGeometry.set(1)
undu_prog_paras.create_z_sym.set(0)

center_coil = uw.undu_magnets.point_coords(x=0.0,y=-20.0,z=0.0)
normal_vec = uw.undu_magnets.point_coords(x=1.0,y=0.0,z=0.0)

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
		center_coords=center_coil,
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
		api=undu
		)
coil.add_to_clc()

undu.run()
results = undu.get_results()

trajx = results.get_result(which='trajx')
by = results.get_result(which='by')
bz = results.get_result(which='bz')
intBy = results.get_result(which='intBy')
intBz = results.get_result(which='intBz')
intBy2 = results.get_result(which='intBy2')
intBz2 = results.get_result(which='intBz2')

add=''

nfig=0
by.plot_over(
	x_quant=trajx,
	nfig=nfig,
	nosave=True,
	clear=True,
	dataFile=f'/{res_folder_full}by{add}.txt'
	)
nfig=bz.plot_over(
	x_quant=trajx,
	nfig=nfig,
	file_name=f'byz{add}.png',
	plot=False,
	title=f'B$_y$ and B$_z$',
	dataFile=	f'/{res_folder_full}bz{add}.txt'
	)

intBy.plot_over(
	x_quant=trajx,
	nfig=nfig,
	nosave=True,
	clear=True,
	dataFile=f'/{res_folder_full}intBy{add}.txt'
	)
nfig=intBz.plot_over(
	x_quant=trajx,
	nfig=nfig,
	file_name=f'intByz{add}.png',
	plot=False,
	title=f'1. Integral of \n B$_y$ and B$_z$',
	dataFile=	f'/{res_folder_full}intBz{add}.txt'
	)

intBy2.plot_over(
	x_quant=trajx,
	nfig=nfig,
	nosave=True,
	clear=True,
	dataFile=f'/{res_folder_full}int2By{add}.txt'
	)
nfig=intBz2.plot_over(
	x_quant=trajx,
	nfig=nfig,
	file_name=f'int2Byz{add}.png',
	plot=False,
	title=f'2. Integral of \n B$_y$ and B$_z$',
	dataFile=	f'/{res_folder_full}int2Bz{add}.txt'
	)

pdb.set_trace()

