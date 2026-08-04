import unduwave as uw 
from unduwave.unduwave_incl import *
from unduwave import undu_blocks
from unduwave import undulatorComponents

try :
	# works when calling script with python3 script_file
	dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
except:
	# works when calling script with exec from python console
	dir_path = Path(os.getcwd())

yamlFile=dir_path/'myUndu.yaml'
res_folder='res'
res_folder_full=dir_path/f'{res_folder}/'

undu = uw.undu(undu_mode='from_undu_magns')
undu_prog_paras = undu._prog_paras
undu_prog_paras.res_folder.set(res_folder_full)
undu_prog_paras.plotGeometry.set(1)
undu_prog_paras.create_z_sym.set(0)
undu_prog_paras.create_y_sym.set(0)
undu_prog_paras.bmap_nz.set(10)
undu_prog_paras.bmap_nx.set(10)
undu_prog_paras.center_magnet_struct.set(0)
undu_prog_paras.nthreads.set(20)
undu_prog_paras.periodLength.set(0.024)

undu_center = np.array([0.0,0.0,0.0])
magn_unit_vec=undu_blocks.create_magnetization_unit_vec(magn_string='+x')

"""
Some general undu parameters
"""
gap = 5.5 # 5.5 
nperiods=3 #73
shift = 0.25

period_length=24
glue_slit= 0.102
keeper_slit_fm= 0.258
row_slit=0.795

magnetization=1.58 

slitLenPerPeriod=2*glue_slit+2*keeper_slit_fm
magnLen=(period_length-slitLenPerPeriod)/4

"""
Create the periodic poles, end poles, periodic magnets and endmagnet 1 and 2
"""

magnParas=undu_blocks.magParameters(
	len_x_main=magnLen, 
	len_y_main=23.5,
	len_z_main=23.5, 
	segm_x=2,
	segm_y=3,
	segm_z=3,
	material_id="pm_rec", 
)
llMagnPos=np.array([
	0.0,
	-magnParas._len_y_main()/2.0,
	-magnParas._len_z_main()/2.0-row_slit/2.0,
	])
magnet = undu_blocks.undumagBlockObject(
	center=llMagnPos,
	magnParas=magnParas,
	pnts=None,
	name="mag",
	parentName='',
	api=None,
	)

"""
Create a period
Period starts at beginning of 1st mag
"""

mag1=magnet.get_copy(name='m1')
mag1.move_it(vec=np.array([
	-period_length/2.0+keeper_slit_fm/2.0+magnParas._len_x_main()/2.0,
	0.0,
	0.0,
	]))

mag2=magnet.get_copy(name='m2')
mag2.move_it(vec=np.array([
	mag1._center[0]+magnParas._len_x_main()+glue_slit,
	0.0,
	0.0,
	]))

mag3=magnet.get_copy(name='m3')
mag3.move_it(vec=np.array([
	mag2._center[0]+magnParas._len_x_main()+keeper_slit_fm,
	0.0,
	0.0,
	]))

mag4=magnet.get_copy(name='m4')
mag4.move_it(vec=np.array([
	mag3._center[0]+magnParas._len_x_main()+glue_slit,
	0.0,
	0.0,
	]))

period_objects=[mag1,mag2,mag3,mag4]

period=undulatorComponents.period(
	period_length=0.024,
	objects=period_objects, # list of objects making up one period
	center=np.array([0.0,0.0,0.0]),
	name='period1'
)

len_upStrm=magnParas._len_x_main()*2.5

endMagn1Paras=copy.deepcopy(magnParas)
endMagn1Paras._len_x_main.set(endMagn1Paras._len_x_main()*0.25)

endMagn2Paras=copy.deepcopy(magnParas)
endMagn2Paras._len_x_main.set(endMagn2Paras._len_x_main()*0.5)

endMagn3Paras=copy.deepcopy(magnParas)
endMagn3Paras._len_x_main.set(endMagn3Paras._len_x_main()*0.75)

endmag1Pos_upStrm=copy.deepcopy(llMagnPos)
endmag1Pos_upStrm[0]=-len_upStrm/2.0+0.5*endMagn1Paras._len_x_main()
endMagn1_upStrm = undu_blocks.undumagBlockObject(
	center=endmag1Pos_upStrm,
	magnParas=endMagn1Paras,
	name="em1UpS",
	)

endmag2Pos_upStrm=copy.deepcopy(llMagnPos)
endmag2Pos_upStrm[0]=endmag1Pos_upStrm[0]+0.5*endMagn1Paras._len_x_main()+\
		magnParas._len_x_main()*0.5+endMagn2Paras._len_x_main()*0.5
endMagn2_upStrm = undu_blocks.undumagBlockObject(
	center=endmag2Pos_upStrm,
	magnParas=endMagn2Paras,
	name="em2UpS",
	)

endmag3Pos_upStrm=copy.deepcopy(llMagnPos)
endmag3Pos_upStrm[0]=endmag2Pos_upStrm[0]+0.5*endMagn2Paras._len_x_main()+\
		magnParas._len_x_main()*0.5+endMagn3Paras._len_x_main()*0.5
endMagn3_upStrm = undu_blocks.undumagBlockObject(
	center=endmag3Pos_upStrm,
	magnParas=endMagn3Paras,
	name="em3UpS",
	)

"""
dwn stream
"""

endmag3Pos_dwnStrm=copy.deepcopy(llMagnPos)
endmag3Pos_dwnStrm[0]=-len_upStrm/2.0+endMagn3Paras._len_x_main()*0.5
endMagn3_dwnStrm = undu_blocks.undumagBlockObject(
	center=endmag3Pos_dwnStrm,
	magnParas=endMagn3Paras,
	name="em3DwnS",
	)

endmag2Pos_dwnStrm=copy.deepcopy(llMagnPos)
endmag2Pos_dwnStrm[0]=endmag3Pos_dwnStrm[0]+0.5*endMagn3Paras._len_x_main()+\
		magnParas._len_x_main()*0.5+endMagn2Paras._len_x_main()*0.5
endMagn2_dwnStrm = undu_blocks.undumagBlockObject(
	center=endmag2Pos_dwnStrm,
	magnParas=endMagn2Paras,
	name="em2DwnS",
	)

endmag1Pos_dwnStrm=copy.deepcopy(llMagnPos)
endmag1Pos_dwnStrm[0]=endmag2Pos_dwnStrm[0]+0.5*endMagn2Paras._len_x_main()+\
		magnParas._len_x_main()*0.5+endMagn1Paras._len_x_main()*0.5
endMagn1_dwnStrm = undu_blocks.undumagBlockObject(
	center=endmag1Pos_dwnStrm,
	magnParas=endMagn1Paras,
	name="em1DwnS",
	)

downstream_end_objects=[
	endMagn3_dwnStrm,endMagn2_dwnStrm,endMagn1_dwnStrm
	]
endDownStream=undulatorComponents.endConfig(
	dist_period=0.0,
	objects=downstream_end_objects, # list of objects making up one period
	center=np.array([0.0,0.0,0.0]),
	name=''
)

upstream_end_objects=[
	endMagn1_upStrm,endMagn2_upStrm,endMagn3_upStrm
	]
endUpStream=undulatorComponents.endConfig(
	dist_period=0.0,
	objects=upstream_end_objects, # list of objects making up one period
	center=np.array([0.0,0.0,0.0]),
	name=''
)

periodicMagnetizationSequence=[
	np.array([ 0.0, magnetization,0.0 ]),
	np.array([ -magnetization, 0.0,0.0 ]),
	np.array([ 0.0, -magnetization, 0.0 ]),
	np.array([ magnetization, 0.0,0.0 ]),
	]
row=undulatorComponents.row(
	pos='ll',
	period=period,
	downstream_end=endDownStream,
	upstream_end=endUpStream,
	period_length=24,
	center=np.array([0.0,0.0,0.0]),
	name='row1',
	nperiods=nperiods,
	periodicMagnetizationSequence=periodicMagnetizationSequence,
	)

undulator=undulatorComponents.undulator(
	name='myUndu',
	period_length=period_length,
	accelerator='bessyII',
	rows=[row],
	center=np.array([0.0,0.0,0.0]),
	gap=6.0,
	gap_range=[6.0,80.0],
	shift=0.25,
	shift_range=[0.0,0.5],
	construct_quadrants=['ll','lr','ul','ur'],
	api=undu,
	)

undulator.writeToParameterFile(file=yamlFile)
periodLength=undulator.get_period_length(name_hints=['ll'])
print(f"periodLength: {periodLength}")

undulatorLoad=undulator.loadFromDict(
	name='myUndu',
	file=yamlFile,
	gap=gap,
	nperiods=nperiods,
	shift=shift,
	)
undu.set_magnet_objects(magn_objects=undulatorLoad)

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

add=''

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
bmap.plot_fld_map(
	bWhat='By', # "Bx", "By" or "Bz"
	xPos=None,
	yPos=0.0,
	zPos=None,
	nfig=nfig,
	filename=res_folder_full/"bymap.png",
	title="By Map",
	)

pdb.set_trace()
