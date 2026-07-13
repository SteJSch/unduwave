
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

res_folder='res'
res_folder_full=dir_path+f'/{res_folder}/'

undu = uw.undu(undu_mode='from_clc_file')
undu_prog_paras = undu._prog_paras
undu_prog_paras.res_folder.set(res_folder+'/')
undu_prog_paras.plotGeometry.set(1)

clc_file='undumag_example.clc'
clc_folder=dir_path+'/'

undu_prog_paras.copy_clc_folder.set(clc_folder)
undu_prog_paras.in_file_clc.set(clc_file)

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