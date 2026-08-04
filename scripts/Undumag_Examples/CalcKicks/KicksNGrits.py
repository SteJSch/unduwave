import os
import pdb
import sys
import unittest
import numpy as np
import getpass
import matplotlib.pyplot as pyplot
import time
import unduwave as uw
from unduwave.unduwave_incl import *
from unduwave import undulatorComponents
from unduwave import undu_blocks
from unduwave import magneticObjectGeometries
import unduwave.analytical_module.ana_undu.grids as grid

from matplotlib.pylab import matplotlib
from matplotlib.ticker import ScalarFormatter

try :
	# works when calling script with python3 script_file
	dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
except:
	# works when calling script with exec from python console
	Path(dir_path = os.getcwd())

if __name__ == "__main__" :

	wave=uw.wave()

	ebeam_paras = wave._ebeam_paras
	ebeam_paras = ebeam_paras.get_std_bessy_III_paras()
	beamEnGeV=ebeam_paras.beam_en()

	res_folder='kicks'
	res_folder_full=dir_path/res_folder
	loadFile=dir_path/"field_map.map"

	processes=6
	period_length=42.5

	bmap=uw.bfield.bfield(
		unitsXB=[1.0,1.0] # setting the units
		)

	bmap.load_field_from_file(
				file=loadFile,
				fieldMap=True,
				unduFile = True, 
			)

	bmap=bmap.center()
	unduRepr=bmap.create_grid_interpolation()
	unduRepr.processes=processes

	for plt in [[0,'Bx'],[1,'By'],[2,'Bz']]:
		unduRepr.plot_grid_map(
			indPlot=plt[0],
			zlab='',
			yPos=0.0,
			save=True,
			filename=res_folder_full/f"map_{plt[1]}.png",
			title=f"{plt[1]}-Map",
			)

	firstIntGrid, scndIntGrid=bmap.calc_integrals_grid(
		intrgnt_limit=400, 
		epsrel=1e-3,
		epsabs=1e-5,
		processes=processes,
		)

	for plt in [[firstIntGrid,'first'],[scndIntGrid,'second']]:
		for indw, what in enumerate(['x','y','z']) :
			plt[0].plot_grid_map(
				indPlot=indw,
				zlab='',
				xPos=plt[0]._g_xvals[-1],
				save=True,
				filename=res_folder_full/f"{plt[1]}_Integral_B{what}.png",
				title=f"{plt[1]} integral",
				)
			plt[0].write_grid_data(
				file=res_folder_full/f"{plt[1]}_Integral_B{what}.dat",
				cols=['x','y','z',f'{plt[1]}_intgrlBx',f'{plt[1]}_intgrlBy',f'{plt[1]}_intgrlBz'],
				)

	beffs, ys, zs, longIntrvl=bmap.calc_beff_grid(
		longIntrvl=[-period_length/2.0,period_length/2.0],
		prd_lngth=period_length,
		intrgnt_limit=None, 
		n_max = 10, 
		epsrel=1e-3,
		epsabs=1e-5,
		processes=processes,
		)

	for plt in [[0,'Bx'],[1,'By'],[2,'Bz']]:
		beffs.plot_grid_map(
			indPlot=plt[0],
			zlab='',
			xPos=period_length,
			save=True,
			filename=res_folder_full/f"{plt[1]}_eff.png",
			title=f"Effective {plt[1]}-Map",
			)

	kicksy_h, kicksz_h, pot_map =bmap.calc_kicks_grid(
		longIntrvl=[-period_length/2.0,period_length/2.0],
		prd_lngth=period_length,
		n_max = 10, 
		epsrel=1e-3,
		epsabs=1e-5,
		processes=processes,
		method='harmonic',
		beamEnGeV=beamEnGeV, # Beam Energy in GeV, beamEnGeV
		)

	for plt in [[kicksy_h,'Ky'],[kicksz_h,'Kz']]:
		plt[0].plot_grid_map(
			indPlot=0,
			zlab='',
			xPos=plt[0]._g_xvals[-1],
			save=True,
			filename=res_folder_full/f"{plt[1]}_harm.png",
			title=f"Harmonic {plt[1]}-Map",
			)

	kicksy, kicksz, pot_map_f, gY,gZ =bmap.calc_kicks_grid(
		prd_lngth=period_length,
		intrgnt_limit=400, 
		n_max = 10, 
		epsrel=1e-3,
		epsabs=1e-5,
		processes=processes,
		method='full',
		beamEnGeV=beamEnGeV, # Beam Energy in GeV, beamEnGeV
		)

	for plt in [[kicksy,'Ky'],[kicksz,'Kz']]:
		plt[0].plot_grid_map(
			indPlot=0,
			zlab='',
			xPos=plt[0]._g_xvals[-1],
			save=True,
			filename=res_folder_full/f"{plt[1]}_full.png",
			title=f"Full {plt[1]}-Map",
			)

	pdb.set_trace()