import unittest

import unduwave as uw
from unduwave.unduwave_incl import *
import unduwave.helpers.file_folder_helpers as f_h

try :
	# works when calling script with python3 script_file
	dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
except:
	# works when calling script with exec from python console
	dir_path = Path(os.getcwd())

class wave_settings_test(unittest.TestCase) :

	def test_wave_settings(self) :

		self.res_folder = dir_path/'res/'

		"""
		Getting wave
		"""
		wave = uw.wave(wave_mode='undu_ellip')
		wave_prog_paras = wave._prog_paras

		"""
		Setting Program Parameter
		"""
		wave_prog_paras.res_folder.set(self.res_folder)
		wave_prog_paras.calc_spectrum.set(True)
		wave_prog_paras.nthreads.set(6)

		"""
		Setting Spectrometer Parameter
		"""

		spectrometer_paras = wave._spectrometer_paras
		spectrometer_paras.spectrum_n_energies.set(10)
		spectrometer_paras.spectrum_min_energy.set(400)
		spectrometer_paras.spectrum_max_energy.set(500)
		spectrometer_paras.spectrum_undu_mode.set(1)

		"""
		Setting Screen Parameter
		"""
		screen_paras = wave._screen_paras
		screen_paras.screen_segm_hor.set(11) 
		screen_paras.screen_segm_vert.set(11)
		screen_paras.screen_extent_hor.set(40) # pinhole width mm
		screen_paras.screen_extent_vert.set(40) # pinhole height mm

		"""
		Setting Undulator Parameter
		"""

		bEffY=round(random.uniform(0,2),4)
		bEffZ=round(random.uniform(0,2),4)
		numPeriods=int(random.uniform(1,40))
		periodLength=round(random.uniform(0.01,0.5),4)
		elliptUnduPerShift=[0.0,0.25,0.5,-0.5,0.75][int(random.uniform(0,5))]

		undu_paras = wave._undu_paras # getting parameter object
		undu_paras.bEffY.set(bEffY)
		undu_paras.bEffZ.set(bEffZ)
		undu_paras.numPeriods.set(numPeriods)
		undu_paras.periodLength.set(periodLength)
		undu_paras.elliptUnduPerShift.set(elliptUnduPerShift)

		"""
		Setting Beam Parameter
		"""
		ebeam_paras = wave._ebeam_paras
		ebeam_paras.get_std_bessy_II_paras()

		"""
		Run
		"""

		wave.run()

		with open(self.res_folder/"wave_res/wave.out", 'r') as o_f:
			wave_out = o_f.readlines()

		res=[]
		for ind,line in enumerate(wave_out) :
			if line.find("number of periods, period and device length [m]:") >= 0 :
				splits=wave_out[ind+1].split(' ')
				for split in splits :
					split=split.strip()
					if len(split) == 0 :
						continue
					res.append(round(float(split),4))
				continue
			if line.find("shift parameter [periods], B_h/B_v:") >= 0 :
				line=line.split("shift parameter [periods], B_h/B_v:")[-1]
				line=line.strip().split(" ")
				shift=float(line[0])
				break
		if elliptUnduPerShift==-0.5:
			self.assertEqual(shift,0.0)
		else :
			self.assertEqual(elliptUnduPerShift,shift)
		self.assertEqual(res[0:2],[numPeriods,periodLength])

		with open(self.res_folder/"wave_res/wave.in", 'r') as o_f:
			wave_in = o_f.readlines()

		self.tearDown()

	def test_wave_easy_ellipt_undu(self) :

		self.res_folder = dir_path/'res/'

		"""
		Getting wave
		"""
		wave = uw.wave(wave_mode='undu_ellip')
		wave_prog_paras = wave._prog_paras

		"""
		Setting Program Parameter
		"""
		wave_prog_paras.res_folder.set(self.res_folder)
		wave_prog_paras.calc_spectrum.set(True)
		wave_prog_paras.nthreads.set(6)

		"""
		Setting Spectrometer Parameter
		"""

		spectrometer_paras = wave._spectrometer_paras
		spectrometer_paras.spectrum_n_energies.set(10)
		spectrometer_paras.spectrum_min_energy.set(400)
		spectrometer_paras.spectrum_max_energy.set(500)
		spectrometer_paras.spectrum_undu_mode.set(1)

		"""
		Setting Screen Parameter
		"""
		screen_paras = wave._screen_paras
		screen_paras.screen_segm_hor.set(11) 
		screen_paras.screen_segm_vert.set(11)
		screen_paras.screen_extent_hor.set(40) # pinhole width mm
		screen_paras.screen_extent_vert.set(40) # pinhole height mm

		"""
		Setting Undulator Parameter
		"""

		beffY=random.uniform(0.0,1.5)
		beffZ=random.uniform(0.0,1.5)
		numPeriods=int(100*random.uniform(0.1,1.5))
		periodLength=random.uniform(0.01,0.03)
		elliptUnduPerShift=[0.0,0.25,0.5,-0.5,0.75][int(random.uniform(0,5))]

		undu_paras = wave._undu_paras # getting parameter object
		undu_paras.bEffY.set(beffY)
		undu_paras.bEffZ.set(beffZ)
		undu_paras.numPeriods.set(numPeriods)
		undu_paras.periodLength.set(periodLength)
		undu_paras.elliptUnduPerShift.set(elliptUnduPerShift)

		"""
		Setting Beam Parameter
		"""
		ebeam_paras = wave._ebeam_paras
		ebeam_paras.get_std_bessy_II_paras()

		"""
		Run
		"""

		wave.run()

		with open(self.res_folder/"wave_res/wave.in", 'r') as o_f:
			wave_in = o_f.readlines()

		res=[]
		for ind,line in enumerate(wave_in) :
			if line.find("KELLIP=") >= 0 :
				line=line.split('KELLIP=')[-1]
				line=line.split("! helical undulator, parameter namelist ELLIPN")[0]
				kellip=float(line)
			elif line.find("B0ELLIPV=") >= 0 :
				line=line.split('B0ELLIPV=')[-1]
				line=line.split("! vertical peak field   [T]")[0]
				beffYW=float(line)
			elif line.find("B0ELLIPH=") >= 0 :
				line=line.split('B0ELLIPH=')[-1]
				line=line.split("! horizontal peak field   [T]")[0]
				beffZW=float(line)
			elif line.find("PERELLIP=") >= 0 :
				line=line.split('PERELLIP=')[-1]
				line=line.split("! number of periods")[0]
				numPeriodsW=int(line)
			elif line.find("XLELLIP=") >= 0 :
				line=line.split('XLELLIP=')[-1]
				line=line.split("! period length [m]")[0]
				periodLengthW=float(line)
			elif line.find("ELLSHFT=") >= 0 :
				line=line.split('ELLSHFT=')[-1]
				line=line.split("! shift of horiz. and vert. fields [periods]")[0]
				elliptUnduPerShiftW=float(line)

		self.assertEqual(kellip,1)

		if elliptUnduPerShift == 0.0 :
			self.assertEqual(np.isclose(beffYW,beffY, rtol=1e-02, atol=1e-02),True)
			self.assertEqual(np.isclose(beffZW,0.0, rtol=1e-02, atol=1e-02),True)
		elif elliptUnduPerShift==0.25 :
			self.assertEqual(np.isclose(beffYW,beffY, rtol=1e-02, atol=1e-02),True)
			self.assertEqual(np.isclose(beffZW,beffZ, rtol=1e-02, atol=1e-02),True)
		elif elliptUnduPerShift==-0.5 : # antiparallel in this crazy world
			self.assertEqual(np.isclose(beffYW,beffY, rtol=1e-02, atol=1e-02),True)
			self.assertEqual(np.isclose(beffZW,beffY, rtol=1e-02, atol=1e-02),True)
		elif elliptUnduPerShift==0.5 :
			self.assertEqual(np.isclose(beffYW,0.0, rtol=1e-02, atol=1e-02),True)
			self.assertEqual(np.isclose(beffZW,beffZ, rtol=1e-02, atol=1e-02),True)
		elif elliptUnduPerShift==0.75:
			self.assertEqual(np.isclose(beffYW,beffY, rtol=1e-02, atol=1e-02),True)
			self.assertEqual(np.isclose(beffZW,beffZ, rtol=1e-02, atol=1e-02),True)

		self.assertEqual(numPeriodsW,numPeriods)
		self.assertEqual(periodLengthW,periodLength)
		if elliptUnduPerShift==-0.5 :
			self.assertEqual(elliptUnduPerShiftW,0.0)
		else :
			self.assertEqual(elliptUnduPerShiftW,elliptUnduPerShift)
		self.tearDown()

	def test_wave_easy_easy_undu(self) :

		self.res_folder = dir_path/'res/'

		wave = uw.wave(wave_mode='undu_easy') # Simple undulator model with endpoles
		wave_prog_paras = wave._prog_paras

		"""
		Setting Undulator Parameter

		We give the K-parameters and the B-Amplitudes are calculated from there
		"""

		k0=0.0
		beffY=random.uniform(0.0,1.5)
		numPeriods=int(100*random.uniform(0.1,1.5))
		periodLength=random.uniform(0.01,0.03)

		undu_paras = wave._undu_paras # getting parameter object
		undu_paras.bEffY.set(beffY)
		undu_paras.numPeriods.set(numPeriods) # we count the number of B-field peaks here - one extra for the end-fields (odd)
		undu_paras.periodLength.set(periodLength)

		"""
		Setting Program Parameter
		"""
		wave_prog_paras.res_folder.set(self.res_folder)
		wave_prog_paras.calc_spectrum.set(True)
		wave_prog_paras.nthreads.set(6)

		"""
		Setting Spectrometer Parameter
		"""

		spectrometer_paras = wave._spectrometer_paras
		spectrometer_paras.spectrum_n_energies.set(10)
		spectrometer_paras.spectrum_min_energy.set(450)
		spectrometer_paras.spectrum_max_energy.set(500)
		spectrometer_paras.spectrum_undu_mode.set(1)

		"""
		Setting Screen Parameter
		"""
		screen_paras = wave._screen_paras
		screen_paras.screen_segm_hor.set(15) 
		screen_paras.screen_segm_vert.set(15)
		screen_paras.screen_extent_hor.set(40) # pinhole width mm
		screen_paras.screen_extent_vert.set(40) # pinhole height mm

		"""
		Setting Beam Parameter
		"""
		ebeam_paras = wave._ebeam_paras
		ebeam_paras = ebeam_paras.get_std_bessy_III_paras()

		"""
		Run
		"""

		wave.run()

		with open(self.res_folder/"wave_res/wave.in", 'r') as o_f:
			wave_in = o_f.readlines()

		res=[]
		for ind,line in enumerate(wave_in) :
			if (line.find(" KHALBA=") >= 0):
				line=line.split('KHALBA=')[-1]
				line=line.split("! insertion device described by HALBACH's formulas")[0]
				KHALBA=float(line)
			elif line.find(" B0HALBA=") >= 0 :
				line=line.split('B0HALBA=')[-1]
				line=line.split("! peak field [T]")[0]
				beffYW=float(line)
			elif line.find("PERHAL=") >= 0 :
				line=line.split('PERHAL=')[-1]
				line=line.split("! number of periods (NOT NECESSARYLY INTEGER)")[0]
				numPeriodsW=float(line)
			elif line.find("      ZLHALBA=") >= 0 :
				line=line.split('ZLHALBA=')[-1]
				line=line.split("! 2*pi/kz [m]")[0]
				periodLengthW=float(line)

		self.assertEqual(KHALBA,1)
		self.assertEqual(np.isclose(beffYW,beffY, rtol=1e-05, atol=1e-05),True)
		self.assertEqual(numPeriodsW,numPeriods)
		self.assertEqual(periodLengthW,periodLength)

		self.tearDown()

	def test_wave_easy_endp_undu(self) :

		self.res_folder = dir_path/'res/'

		beffY=random.uniform(0.0,1.5)
		numPeriods=int(100*random.uniform(0.1,1.5))
		periodLength=random.uniform(0.01,0.03)

		wave = uw.wave(wave_mode='undu_endp') # Simple undulator model with endpoles
		wave_prog_paras = wave._prog_paras

		"""
		Setting Undulator Parameter

		We give the K-parameters and the B-Amplitudes are calculated from there
		"""

		undu_paras = wave._undu_paras # getting parameter object
		undu_paras.bEffY.set(beffY)
		undu_paras.numPeriods.set(numPeriods) # number full periods
		undu_paras.periodLength.set(periodLength)

		"""
		Setting Program Parameter
		"""
		wave_prog_paras.res_folder.set(self.res_folder)
		wave_prog_paras.calc_spectrum.set(True)
		wave_prog_paras.nthreads.set(6)

		"""
		Setting Spectrometer Parameter
		"""

		spectrometer_paras = wave._spectrometer_paras
		spectrometer_paras.spectrum_n_energies.set(10)
		spectrometer_paras.spectrum_min_energy.set(450)
		spectrometer_paras.spectrum_max_energy.set(500)
		spectrometer_paras.spectrum_undu_mode.set(1)

		"""
		Setting Screen Parameter
		"""
		screen_paras = wave._screen_paras
		screen_paras.screen_segm_hor.set(15) 
		screen_paras.screen_segm_vert.set(15)
		screen_paras.screen_extent_hor.set(40) # pinhole width mm
		screen_paras.screen_extent_vert.set(40) # pinhole height mm

		"""
		Setting Beam Parameter
		"""
		ebeam_paras = wave._ebeam_paras
		ebeam_paras = ebeam_paras.get_std_bessy_III_paras()

		"""
		Run
		"""

		wave.run()

		with open(self.res_folder/"wave_res/wave.in", 'r') as o_f:
			wave_in = o_f.readlines()

		res=[]
		for ind,line in enumerate(wave_in) :
			if (line.find("      KHALBASY=") >= 0):
				line=line.split('KHALBASY=')[-1]
				line=line.split("! simple WLS or undulator/wiggler model")[0]
				KHALBA=float(line)
			elif line.find(" B0HALBASY=") >= 0 :
				line=line.split('B0HALBASY=')[-1]
				line=line.split("! peak field [T]")[0]
				beffYW=float(line)
			elif line.find("AHWPOL=") >= 0 :
				line=line.split('AHWPOL=')[-1]
				line=line.split("! number of main poles, should be an odd number")[0]
				numPeriodsW=(float(line)-1)/2
			elif line.find("      ZLHALBASY=") >= 0 :
				line=line.split('ZLHALBASY=')[-1]
				line=line.split("! 2*pi/kz of main poles [m]")[0]
				periodLengthW=float(line)

		self.assertEqual(KHALBA,1)
		self.assertEqual(np.isclose(beffYW,beffY, rtol=1e-05, atol=1e-05),True)
		self.assertEqual(numPeriodsW,numPeriods)
		self.assertEqual(periodLengthW,periodLength)

		self.tearDown()

	def setUp(self):
		self.resource = "Resource allocated"
		print("Setting up test resources...")

	def tearDown(self):
		print("Cleaning up resource...")

		ROOT_DIR_TEST = Path(__file__).resolve().parent
		del_dirs=[
			ROOT_DIR_TEST/'__pycache__',
			ROOT_DIR_TEST/self.res_folder,
			]

		for del_dir in del_dirs :
			if os.path.exists(del_dir) :		
				shutil.rmtree(del_dir)

if __name__ == '__main__':
	unittest.main()
