from unduwave.unduwave_incl import *
from unduwave.wave_modules.wave_parameters import *
import unduwave.helpers.file_folder_helpers as f_h
import unduwave.analytical_module.ana_undu.bfield as bfield

class undu_postprocess:
	"""
	A class for postprocessing Undumag files.
	
	Args:
		undu_api 
	"""
	
	def __init__(self, undu_api):
		"""
		Takes external api and create postprocess class
		"""
		self._undu_api = undu_api

	def copy_results(self,add=''):
		"""
		Cleans the undumag-stage folder and copies the desired files to their location,
		deletes non-desired files, and zips the results based on the undu_res_copy_behaviour setting.
		"""
		#'copy_all', 'copy_none' - only writes res_summary, 'copy_essentials' 
		undu_folder    = self._undu_api._prog_paras.undumag_curr_folder.get()
		res_folder     = self._undu_api._prog_paras.res_folder.get()
		res_undu       = res_folder / self._undu_api._prog_paras.undu_data_res_folder.get()
		os.makedirs(res_folder, exist_ok=True)
		os.makedirs(res_undu, exist_ok=True)
		copy_behav     = self._undu_api._prog_paras.undu_res_copy_behaviour.get()
		zip_res_folder = self._undu_api._prog_paras.zip_res_folder.get()
		files_dont_del = []

		if copy_behav == 'copy_all' :
			undu_res_extract   = self._undu_api._prog_paras.undu_ending_extract.get()
			undu_res_copy      = self._undu_api._prog_paras.undu_ending_copy.get()
			undu_files_no_copy = self._undu_api._prog_paras.no_copy.get()
			files_del = []
		elif copy_behav == 'copy_del_none' :
			undu_res_extract   = []
			undu_res_copy      = []
			undu_files_no_copy = []
			files_del          = []
		elif copy_behav == 'copy_essentials' :
			undu_res_extract   = self._undu_api._prog_paras.undu_files_essentials.get()
			undu_res_copy      = self._undu_api._prog_paras.undu_ending_copy.get()
			undu_files_no_copy = []
			files_del          = self._undu_api._prog_paras.undu_ending_extract.get()
			files_dont_del     = self._undu_api._prog_paras.no_copy.get()

		if not os.path.exists(res_undu):
			os.makedirs(res_undu)

		files_moved = f_h.mv_cp_files(hints = undu_res_extract, exptns = undu_files_no_copy\
				, folder_in = undu_folder / 'stage/', folder_out = res_undu, move = True, add_string = add)
		files_copied = f_h.mv_cp_files(hints = undu_res_copy, exptns = undu_files_no_copy\
				, folder_in = undu_folder / 'stage/', folder_out = res_undu, move = False, add_string = add)
		f_h.del_files(hints = files_del, exptns = files_dont_del, folder = undu_folder / 'stage/')

		for file_m in files_moved :
			if file_m.find('undumag_msh_radia') >= 0 :
				self.process_radia_out(radia_file=res_undu/file_m,add=add)
		# if zip_res_folder : 
		#     f_h.zip_files_in_folder(folder_to_pack = res_undu)
		return files_moved + files_copied

	def process_radia_out(self,radia_file,add='') :

		prog_paras=self._undu_api._prog_paras
		with open(radia_file, "r") as file:
			content = file.readlines()

		periodLength=prog_paras.periodLength()*1e3
		zsym=prog_paras.create_z_sym()
		xsym=prog_paras.create_x_sym()
		ysym=prog_paras.create_y_sym()

		dx=prog_paras.bmap_dx()
		zmin=prog_paras.bmap_z_min()
		zmax=prog_paras.bmap_z_max()
		nz=prog_paras.bmap_nz()
		xmin=prog_paras.bmap_x_min()
		xmax=prog_paras.bmap_x_max()
		nx=prog_paras.bmap_nx()
		ymin=prog_paras.bmap_y_min()
		ymax=prog_paras.bmap_y_max()
		ny=prog_paras.bmap_ny()

		with open( radia_file, 'w') as o_f:
			for ind, line in enumerate(content) :

				if line.find('open("unduradia.map",')>=0:

					o_f.write(f"  dir_path = os.path.dirname(os.path.realpath(__file__))\n")
					o_f.write(f"  radiaResPath=dir_path+'/radia_results/'\n")
					o_f.write(f"  os.makedirs(radiaResPath, exist_ok=True)\n")
					o_f.write(f"  kickMap = rad.FldFocKickPer(UnduSetUp, [0,0,0],[0,1,0],{periodLength:.4f},1,[1,0,0],{(zmax-zmin):.4f},{nz},{(ymax-ymin):.4f},{ny},'fileheader',[11,50,0,0],'rad')\n")
					o_f.write(f"  kickMapFile=kickMap[4]\n")
					o_f.write(f"  with open(radiaResPath+f'radiaKick{add}.map','w') as fo :\n")
					o_f.write(f"      for line in kickMapFile :\n")
					o_f.write(f"        fo.write(line)\n")
					line=f'  FILEMAP = open(radiaResPath+"radia_map{add}.map","w")\n'

				if line.find("nUnduNoPolMap = ")>=0:
					line=f"nUnduNoPolMap =            1 \n"
				elif line.find("nUnduNoMagMap = ")>=0:
					line=f"nUnduNoMagMap =            1 \n"
				elif line.find("iUnduXsym =")>=0:
					line=f"iUnduXsym ={xsym}\n"
				elif line.find("iUnduYsym =")>=0:
					line=f"iUnduYsym ={ysym}\n"
				elif line.find("iUnduZsym =")>=0:
					line=f"iUnduZsym ={zsym}\n"
				elif line.find("msh.") >= 0:
					line=line.replace("msh.","")
				elif line.find('rad.ObjDivMag(') >= 0 :
					endL=",'kxkykz->Numb')"
					begin="rad.ObjDivMag("
					spltLine=line.split(endL)[0]
					spltLine=spltLine.split(begin)[-1]
					objName=spltLine.split(',')[0]
					segmList=eval('[['+spltLine.split('[[')[-1])
					segmList[0][1]=1/segmList[0][1]
					segmList[2][1]=1/segmList[2][1]
					newLine=begin+objName+','+str(segmList)+endL+'\n'
					line=newLine
				o_f.write(line)
				if line.find('from __future__') >= 0 :

					o_f.write("import sys\n")
					o_f.write("import getpass\n")
					o_f.write("user=getpass.getuser()\n")
					o_f.write("sys.path.insert(0,f'/home/{user}/Documents/Very_Unimportant/Arbeit/HZB/Undus/Programming/Python_Code/apy/unduwave/External-Software/UNDUMAG/python')\n")
					o_f.write("sys.path.insert(0,f'/home/{user}/Documents/Very_Unimportant/Arbeit/HZB/Undus/Software/Radia/env/radia_python')\n")
	
	def cleanup(self,touchedFiles=[],add=''):
		"""
		Cleans up the Undumag run

		"""
		self.clear_map_file(touchedFiles=touchedFiles,add=add)

	def clear_map_file(self,touchedFiles=[],add='') :

		mapFile=None
		for file in touchedFiles :
			if file.find('.map') >= 0 :
				mapFile=file
				break
		if mapFile is None :
			return

		prog_paras=self._undu_api._prog_paras
		res_folder     = prog_paras.res_folder.get()
		res_undu       = res_folder / prog_paras.undu_data_res_folder.get()

		mapFile=res_undu/mapFile
		if not (mapFile is None) :
			if os.path.isfile(mapFile) :

				horInt=np.linspace(prog_paras.bmap_z_min(),prog_paras.bmap_z_max(),prog_paras.bmap_nz())
				vertInt=np.linspace(prog_paras.bmap_y_min(),prog_paras.bmap_y_max(),prog_paras.bmap_ny())

				bmap=bfield.bfield(
					unitsXB=[1e-3,1.0] # setting the units
					)
				mapOutFile=str(mapFile).split('.map')[0]+'_cleared.map'
				bmap.load_clear_undu_map(
					mapFile=mapFile,
					vHor=horInt,
					vVert=vertInt, 
					mapOutName=mapOutFile,
					)

				shutil.move(f_h.convert_path_to_win(mapOutFile), f_h.convert_path_to_win(mapFile))
