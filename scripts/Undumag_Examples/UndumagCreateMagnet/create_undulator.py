
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
from unduwave import undulatorComponents

import numpy as np

res_folder='res'
res_folder_full=dir_path+f'/{res_folder}/'

undu = uw.undu(undu_mode='from_undu_magns')
undu_prog_paras = undu._prog_paras
undu_prog_paras.res_folder.set(res_folder+'/')
undu_prog_paras.plotGeometry.set(1)
undu_prog_paras.create_y_sym.set(1)
undu_prog_paras.bmap_nz.set(10)
undu_prog_paras.bmap_nx.set(10)
undu_prog_paras.center_magnet_struct.set(0)

period_length=20
len_y_main=20

magn_pos=np.array([ 0.0, -len_y_main/2, 0.0 ])
magn_paras=undu_blocks.magParameters(
	len_x_main=period_length/4, 
	len_y_main=len_y_main,
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
	center=magn_pos,
	magnParas=magn_paras,
	name='magnet',
	parentName='',
	)

mag1=magnet.get_copy(name='m1')
mag1.move_it(vec=np.array([
	-period_length/2.0+magn_paras._len_x_main()/2.0,
	0.0,
	0.0,
	]))

mag2=magnet.get_copy(name='m2')
mag2.move_it(vec=np.array([
	mag1._center[0]+magn_paras._len_x_main(),
	0.0,
	0.0,
	]))

mag3=magnet.get_copy(name='m3')
mag3.move_it(vec=np.array([
	mag2._center[0]+magn_paras._len_x_main(),
	0.0,
	0.0,
	]))

mag4=magnet.get_copy(name='m4')
mag4.move_it(vec=np.array([
	mag3._center[0]+magn_paras._len_x_main(),
	0.0,
	0.0,
	]))

period_objects=[mag1,mag2,mag3,mag4]

period=undulatorComponents.period(
	period_length=period_length*1e-3,
	objects=period_objects, # list of objects making up one period
	center=np.array([0.0,0.0,0.0]),
	name='period1'
)

magnetizations=[1.1,1.2]
periodicMagnetizationSequence=[
	np.array([ 0.0, magnetizations[1],0.0 ]),
	np.array([ -magnetizations[0], 0.0,0.0 ]),
	np.array([ 0.0, -magnetizations[1], 0.0 ]),
	np.array([ magnetizations[0], 0.0,0.0 ]),
	]
row=undulatorComponents.row(
	period=period,
	period_length=period_length,
	center=np.array([0.0,0.0,0.0]),
	name='row1',
	nperiods=3,
	periodicMagnetizationSequence=periodicMagnetizationSequence,
	)

undulator=undulatorComponents.undulator(
	name='ivue20',
	accelerator='',
	rows=[row],
	center=np.array([0.0,0.0,0.0]),
	gap=10.0,
	symmetries=['y','z'],
	period_length=period_length,
	)

undu.set_magnet_objects(magn_objects=undulator)
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

