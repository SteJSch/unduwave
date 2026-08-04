import unduwave as uw
from unduwave.unduwave_incl import *
from unduwave import undu_list_loading

try :
	# works when calling script with python3 script_file
	dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
except:
	# works when calling script with exec from python console
	dir_path = Path(os.getcwd())

add=""

"""
Leave one of the following 3-line blocks uncommented to 
simulate the specific file
"""

# res_folder='u24'
# official_yaml_folder=dir_path/'U24_YAML/'
# official_yaml_file='u24.yaml'

# res_folder='ue30'
# official_yaml_folder=dir_path/'UE30_YAML/'
# official_yaml_file='ue30.yaml'

res_folder='hyu24'
official_yaml_folder=dir_path/'HYU24_YAML/'
official_yaml_file='hyu24.yaml'

gap=6
shift=0.0
nperiods=3

res_folder_full=dir_path/f'{res_folder}/'

undu = uw.undu(undu_mode='from_undu_magns')
undu_prog_paras = undu._prog_paras
undu_prog_paras.res_folder.set(res_folder_full)
undu_prog_paras.plotGeometry.set(1)
undu_prog_paras.bmap_nz.set(1)
undu_prog_paras.bmap_ny.set(1)
undu_prog_paras.bmap_nx.set(1)
undu_prog_paras.center_magnet_struct.set(0)
undu_prog_paras.nthreads.set(6)

ule=undu_list_loading.undu_list_element(
	listFile=official_yaml_file,
	listFolder=official_yaml_folder,
	)

undulator=ule.constructUndulator(
	gap=gap,
	shift=shift,
	center=None,#np.array([-200.0,0.0,0.0]),
	nperiods=nperiods,
	api=undu,
	)

periodLength=undulator.get_period_length(name_hints=['ll'])

undu_prog_paras.periodLength.set(periodLength*1e-3)
if 'X' in undulator._symmetries :
	undu_prog_paras.create_z_sym.set(1)
if 'Y' in undulator._symmetries :
	undu_prog_paras.create_y_sym.set(1)

undu.set_magnet_objects(magn_objects=undulator)

# undu.run(add=add)
results = undu.get_results(add=add)

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
	dataFile=f'by{add}.txt'
	)
nfig=bz.plot_over(
	x_quant=trajx,
	nfig=nfig,
	file_name=f'byz{add}.png',
	plot=True,
	title=f'B$_y$ and B$_z$',
	dataFile=	f'bz{add}.txt'
	)

intBy.plot_over(
	x_quant=trajx,
	nfig=nfig,
	nosave=True,
	clear=True,
	dataFile=f'intBy{add}.txt'
	)
nfig=intBz.plot_over(
	x_quant=trajx,
	nfig=nfig,
	file_name=f'intByz{add}.png',
	plot=False,
	title=f'1. Integral of \n B$_y$ and B$_z$',
	dataFile=	f'intBz{add}.txt'
	)

intBy2.plot_over(
	x_quant=trajx,
	nfig=nfig,
	nosave=True,
	clear=True,
	dataFile=f'int2By{add}.txt'
	)
nfig=intBz2.plot_over(
	x_quant=trajx,
	nfig=nfig,
	file_name=f'int2Byz{add}.png',
	plot=False,
	title=f'2. Integral of \n B$_y$ and B$_z$',
	dataFile=	f'int2Bz{add}.txt'
	)

ByProfile.plot_over(
	x_quant=zProfile,
	nfig=nfig,
	nosave=True,
	clear=True,
	dataFile=f'profileBy{add}.txt'
	)
nfig=BzProfile.plot_over(
	x_quant=zProfile,
	nfig=nfig,
	file_name=f'profileBz{add}.png',
	plot=False,
	title=f'Profile of B$_y$ and B$_z$',
	dataFile=	f'profileBz{add}.txt'
	)

pdb.set_trace()
