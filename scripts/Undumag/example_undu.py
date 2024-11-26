
import pdb
import sys
import os
sys.path.insert(0, '../../')

try :
	# works when calling script with python3 script_file
	dir_path = os.path.dirname(os.path.realpath(__file__))
except:
	# works when calling script with exec from python console
	dir_path = os.getcwd()

import unduwave as uw
from unduwave.api.undu_api_root import create_magnetization_unit_vec
from unduwave.api.undu_api_root import point_coords
from unduwave.api.undu_api_root import undu_magnet_block_coords

undu = uw.undu(res_folder=dir_path+'/res/')
undu_prog_paras = undu._undu_prog_paras

undu_center = point_coords()
magn_unit_vec=create_magnetization_unit_vec(magn_string='+x')

main_block = undu_magnet_block_coords(
	p_center=undu_center,
	len_x=10,
	len_y=20,
	len_z=20,
	magnetization=1.12,
	magn_unit_vec=magn_unit_vec,
	name="mag",
	mother="mom",
	segm_x = 2, 
	segm_y = 2, 
	segm_z = 2, 
	frac_y=1,
	frac_z=1,
	material = 'mag',
	chamf = 0.3,
	)
undu.add_element(element=main_block)

undu.run()

pdb.set_trace()

