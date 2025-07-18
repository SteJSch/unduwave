
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
from unduwave import undu_magnets

res_folder='res'
res_folder_full=dir_path+f'/{res_folder}/'

undu = uw.undu(undu_mode='from_clc_file')
undu_prog_paras = undu._prog_paras
undu_prog_paras.res_folder.set(res_folder+'/')

clc_file='undumag_example.clc'
clc_folder=dir_path+'/'

undu_prog_paras.copy_clc_folder.set(clc_folder)
undu_prog_paras.in_file_clc.set(clc_file)

# undu_center = undu_magnets.point_coords(x=0.0,y=-20,z=0)
# magn_unit_vec=undu_magnets.create_magnetization_unit_vec(magn_string='+x')

# main_block = undu_magnets.undu_magnet_block_coords(
# 	p_center=undu_center,
# 	len_x=10,
# 	len_y=20,
# 	len_z=20,
# 	magnetization=1.12,
# 	magn_unit_vec=magn_unit_vec,
# 	name="mag",
# 	mother="mom",
# 	segm_x = 2, 
# 	segm_y = 2, 
# 	segm_z = 2, 
# 	frac_y=1,
# 	frac_z=1,
# 	material = 'mag',
# 	chamf = 0.3,
# 	api=undu
# 	)
# undu.add_element(element=main_block)

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

