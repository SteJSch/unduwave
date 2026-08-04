import unittest

import unduwave as uw
from unduwave.unduwave_incl import *
import unduwave.helpers.file_folder_helpers as f_h
from unduwave import undu_list_loading

try :
	# works when calling script with python3 script_file
	dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
except:
	# works when calling script with exec from python console
	dir_path = Path(os.getcwd())

class unduList_test(unittest.TestCase) :

	def test_all_list_loads(self) :

		self.res_folders = []
		self.in_folder = dir_path/'yamls/'

		files=f_h.find_files_exptn(folder=self.in_folder, hints=['.yaml'], exptns=['_geometries'])

		for file in files:

			res_folder=file.split('.yaml')[0]
			res_folder=res_folder+'_res'
			res_folder_full=dir_path/f'{res_folder}/'
			self.res_folders.append(res_folder_full)

			official_yaml_folder=self.in_folder
			official_yaml_file=file

			ule=undu_list_loading.undu_list_element(
				listFile=official_yaml_file,
				listFolder=official_yaml_folder,
				)

			undulator=ule.constructUndulator(
				gap=None,
				shift=0.0,
				center=None,#np.array([-200.0,0.0,0.0]),
				nperiods=3,
				)

			undu = uw.undu(undu_mode='from_undu_magns')
			undu_prog_paras = undu._prog_paras
			undu_prog_paras.res_folder.set(res_folder_full)
			undu_prog_paras.plotGeometry.set(1)
			if 'X' in undulator._symmetrie :
				undu_prog_paras.create_z_sym.set(1)
			if 'Y' in undulator._symmetrie :
				undu_prog_paras.create_y_sym.set(1)
			undu_prog_paras.bmap_nz.set(1)
			undu_prog_paras.bmap_ny.set(1)
			undu_prog_paras.bmap_nx.set(10)
			undu_prog_paras.center_magnet_struct.set(0)
			undu_prog_paras.nthreads.set(6)
			undulator.add_to_clc(api=undu)

			undu.run(add='')


		pdb.set_trace()
		self.tearDown()

	def setUp(self):
		self.resource = "Resource allocated"
		print("Setting up test resources...")

	def tearDown(self):
		print("Cleaning up resource...")

		ROOT_DIR_TEST = Path(__file__).resolve().parent
		del_dirs=[
			ROOT_DIR_TEST/'__pycache__',
			]
		del_dirs=del_dirs+self.res_folders

		for del_dir in del_dirs :
			if os.path.exists(del_dir) :		
				shutil.rmtree(del_dir)

if __name__ == '__main__':
	unittest.main()
