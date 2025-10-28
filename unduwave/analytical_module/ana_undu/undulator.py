"""
Defining Undulator
"""

import unduwave as unduwave
from unduwave.unduwave_incl import *
from unduwave.attribute_classes.attributes import _attribute
from unduwave.attribute_classes.attributes import _attribute_collection
import unduwave.helpers.file_folder_helpers as f_h
import unduwave.constants as uc
import unduwave.analytical_module.ana_undu.bfield as bfield

def getUnduObjFromBField(bfield) :
	"""
	returns undulatorObj created from the field
	"""
	pass

class undulatorObj(_attribute_collection) :

	def __init__(self,ebeam=None):
		self._undulatorGeoMatParas=None
		self._undulatorCharacter=None
		self._undulatorGeoMatRepres=None
		self._ebeam=ebeam
		self._anaCharacter=undulatorCharacterization()
		self._numCharacter=undulatorCharacterization()

	def setEbeam(self,ebeam) :
		self._ebeam=ebeam

	def fromBField(self,bfield) :
		"""
		tries to get paras from the bfield? bmax bmin beff lambda
		"""
		pass

	def fromRepres(self,
			undulatorRepres,
			periodLength=None,
			nPeriods=None,
			unduAPI=None,
			setSym=True,
			simulate=True,
			beffs=None,
			resFolder=None,
			) :
		"""
		tries to get data from the clc-representation for undumag?
		"""
		self._undulatorGeoMatRepres=undulatorRepres
		if periodLength is None :
			periodLength=self._undulatorGeoMatRepres.get_period_length()
		center, maxs, mins=self._undulatorGeoMatRepres.get_max_extent(maxs=None,mins=None)
		lenX=maxs[0]-mins[0]
		if nPeriods is None :
			nPeriods=int(lenX/periodLength)

		if (beffs is None) and simulate :

			if unduAPI is None :
				unduAPI = unduwave.undu(undu_mode='from_undu_magns')
			if resFolder is None :
				print("ERROR: fromRepres: In order to simulate, you have to set a result folder")
				return
			self.setUndumag(unduAPI=unduAPI,setSym=setSym)
			undu_prog_paras = unduAPI._prog_paras
			undu_prog_paras.res_folder.set(resFolder)
			undu_prog_paras.plotGeometry.set(1)
			self._undulatorGeoMatRepres.add_to_clc(api=unduAPI)
			# unduAPI.run()
			results = unduAPI.get_results()
			trajx = results.get_result(which='trajx')
			by = results.get_result(which='by')
			bz = results.get_result(which='bz')
			bfield=unduwave.bfield.bfield(
				unitsXB=[0.001,1.0] # setting the units
				)
			bfield.fromQuants(
				xQuant=trajx,
				byQuant=by,
				bzQuant=bz,
				)
			# bfield.plot_fld()
			beffFull, beffs = bfield.calc_beff(prd_lngth=periodLength, n_max = 20,colx = 'x')

		self._undulatorCharacter=undulatorCharacterization(
			bEffY=beffs[0],
			periodLength=periodLength*1e-3,
			numPeriods=nPeriods,
			bEffZ=beffs[1],
			lengthEndPeriodsRelative=0.0,
			ebeam=self._ebeam,
			thetaObservation=0.0,
			unduKY=0.0,
			unduKZ=0.0
			)

	def fromCharacter(self,undulatorRepres) :
		"""
		tries to get data from the characterization
		"""
		pass

	def fromUnduList(self,unduName) :
		"""
		loading undu from bessyII or III list, check if there is
		a model to load, or create one simple model
		"""
		pass

	def getGeoMatRepres(self) :
		"""
		Creates some easy representation in clc form of undulator from charact.
		for this, need some scan b_eff(g,lamb) - nope, lamb set and beff is not(!!)
		guaranteed. depends also on magnetizations.
		anyway, get something close to beff with some std values?!
		"""
		pass

	def getAnaSpectrum(self) :
		"""
		creates aspectrum object with self-paras, returns that
		"""
		pass

	def setUndumag(self,unduAPI,setSym=True) :
		"""
		sets undumag paras, bfield range, symmetries?
		since we want to simulate with undumag, we need the undulatorRepres
		"""
		if self._undulatorGeoMatRepres is None :
			print("No undulator representation set")
			return
		periodLength=self._undulatorGeoMatRepres.get_period_length()
		center, maxs, mins=self._undulatorGeoMatRepres.get_max_extent(maxs=None,mins=None)
		lenX=maxs[0]-mins[0]
		undu_prog_paras = unduAPI._prog_paras
		undu_prog_paras.periodLength.set(periodLength*1e-3)
		undu_prog_paras.bmap_x_min.set(mins[0]-5*periodLength)
		undu_prog_paras.bmap_x_max.set(maxs[0]+5*periodLength)
		if setSym :
			if mins[1] <= 0 :
				if maxs[1] <= 0 :
					undu_prog_paras.create_y_sym.set(1)
			elif mins[1] >= 0 :
				if maxs[1] >= 0 :
					undu_prog_paras.create_y_sym.set(1)
			if mins[2] <= 0 :
				if maxs[2] <= 0 :
					undu_prog_paras.create_z_sym.set(1)
			elif mins[2] >= 0 :
				if maxs[2] >= 0 :
					undu_prog_paras.create_z_sym.set(1)

	def setWAVE(self,wave) :
		"""
		make the wave input for some generic simu
		setting energies around 1st harm
		"""
		pass

class undulatorCharacterization(_attribute_collection) :

	def __init__(self,
			bEffY=0.0,
			periodLength=0.0,
			numPeriods=0.0,
			bEffZ=0.0,
			lengthEndPeriodsRelative=0.0,
			ebeam=None,
			thetaObservation=0.0,
			unduKY=0.0,
			unduKZ=0.0
			) :
		"""
		Definition
		"""
		self.bEffY = _attribute(1)
		self.bEffZ = _attribute(0)
		self.bEff = _attribute(1)
		self.periodLength = _attribute(0.02)
		self.numPeriods = _attribute(78)
		self.lengthEndPeriodsRelative = _attribute(2*1.5)

		self.undulatorLength = _attribute()
		self.undulatorWaveNum = _attribute()
		self.undulatorParameterKY = _attribute()
		self.undulatorParameterKZ = _attribute()
		self.undulatorParameterK = _attribute()
		self.reducedBeta = _attribute()

		self.undulatorFrequency = _attribute()
		self.undulatorPeriodTime = _attribute()
		self.reducedGamma = _attribute()
		self.reducedUndulatorParameterK = _attribute()
		self.maxDeflection = _attribute()
		self.maxAngle = _attribute()
		self.naturalOpeningAngle = _attribute()
		self.properFrequencyWeak = _attribute()
		self.properFrequencyStrong = _attribute()
		self.firstHarmEnergyWeak = _attribute()
		self.firstHarmEnergyStrong = _attribute()
		self.thetaObservation = _attribute()
		self.ebeam = _attribute()

		"""
		Setting
		"""
		self.bEffY.set(bEffY)
		self.bEffZ.set(bEffZ)
		self.periodLength.set(periodLength)
		self.numPeriods.set(numPeriods)
		self.lengthEndPeriodsRelative.set(lengthEndPeriodsRelative)
		self.thetaObservation.set(thetaObservation)
		self.undulatorParameterKY.set(unduKY)
		self.undulatorParameterKZ.set(unduKZ)
		self.ebeam.set(ebeam)
		self.update_non_beam_values()
		if not (self.ebeam.get() is None) :
			self.update_values(
				ebeam=self.ebeam.get(),
				thetaObservation=self.thetaObservation.get()
				)
		super().__init__()

	def update_non_beam_values(self) :
		if self.periodLength.get() == 0.0 :
			return
		self.undulatorLength.set( self.periodLength.get()*(self.numPeriods.get()+self.lengthEndPeriodsRelative.get()) )
		self.undulatorWaveNum.set( 2 * math.pi / ( self.periodLength.get() ) )
		if not(self.undulatorParameterKY.get() == 0.0) :
			self.bEffY.set(self.undulatorParameterKY.get()*uc.m_el * uc.v_c * self.undulatorWaveNum.get()/uc.q_el)
		elif not (self.bEffY.get() == 0.0 ) :
			ky=uc.q_el * self.bEffY.get() / ( uc.m_el * uc.v_c * self.undulatorWaveNum.get() )
			self.undulatorParameterKY.set(ky)
		if not(self.undulatorParameterKZ.get() == 0.0) :
			self.bEffZ.set(self.undulatorParameterKZ.get()*uc.m_el * uc.v_c * self.undulatorWaveNum.get()/uc.q_el)
		elif not (self.bEffZ.get() == 0.0 ) :
			kz=uc.q_el * self.bEffZ.get() / ( uc.m_el * uc.v_c * self.undulatorWaveNum.get() )
			self.undulatorParameterKZ.set(kz)
		unduk=math.sqrt(self.undulatorParameterKY.get()**2+self.undulatorParameterKZ.get()**2)
		self.undulatorParameterK.set(unduk)
		beff=math.sqrt(self.bEffY.get()**2+self.bEffZ.get()**2)
		self.bEff.set(beff)

	def update_values(self,ebeam,thetaObservation=0.0) :
		self.update_non_beam_values()
		if ebeam is None :
			ebeam=self.ebeam.get()
		if (ebeam is None) or (self.undulatorWaveNum.get() is None) :
			return
		gamma=ebeam.gammaFactor.get()
		self.reducedGamma.set(gamma/math.sqrt(1+self.undulatorParameterK.get()**2/2))
		rgamma=self.reducedGamma.get()
		self.reducedBeta.set( ebeam.betaFactor.get() * ( 1 - self.undulatorParameterK.get()**2/(4*ebeam.betaFactor.get()**2 * gamma**2) ))
		self.undulatorFrequency.set( self.undulatorWaveNum.get()*uc.v_c*self.reducedBeta.get())
		self.undulatorPeriodTime.set( 2*math.pi / self.undulatorFrequency.get())
		self.reducedUndulatorParameterK.set( self.undulatorParameterK.get() / math.sqrt( 1 + self.undulatorParameterK.get()**2/2 ))
		self.maxDeflection.set( self.undulatorParameterK.get()/(ebeam.betaFactor.get()*gamma*self.undulatorWaveNum.get()))
		self.maxAngle.set( self.undulatorParameterK.get()/gamma)
		if not (self.reducedUndulatorParameterK.get() == 0.0) :
			self.naturalOpeningAngle.set( self.maxAngle.get() / self.reducedUndulatorParameterK.get())
		else:
			self.naturalOpeningAngle.set( 0.0 )
		self.properFrequencyWeak.set( 2 * gamma**2/(1+gamma**2*thetaObservation**2)*self.undulatorFrequency.get())
		self.properFrequencyStrong.set( 2 * rgamma**2/(1+rgamma**2*thetaObservation**2)*self.undulatorFrequency.get())
		self.firstHarmEnergyWeak.set( self.properFrequencyWeak.get()*uc.hbar/uc.q_el )
		self.firstHarmEnergyStrong.set( self.properFrequencyStrong.get()*uc.hbar/uc.q_el )

	def set_paras_from_bessyII_undu_list(self, unduName) :
		fullList=f_h.loadBessyIIundulatorList()
		indFnd=None
		for indEl, el in enumerate(fullList) : 
			if el['name'] == unduName :
				indFnd=indEl
				break
		if not (indFnd is None) :
			unduFnd=fullList[indFnd]
			self.bEffY.set(unduFnd['beffy [T]'])
			if unduFnd['type'] == 'Apple II' : # elliptical undu
				self.bEffZ.set(unduFnd['beffz [T]'])
			self.periodLength.set(unduFnd['period length [mm]']*1e-3)
			self.numPeriods.set(unduFnd['periods'])
			self.update_values(
				thetaObservation=self.thetaObservation.get(),
				ebeam=self.ebeam.get(),
				)
		else:
			print(f"ana_undulator: set_paras_from_undu_list: Undulator {unduName} not found in list.")

	def getHarmonicBField(self,unitsXB=[1.0,1.0],wEnds=True) :
		bfieldHarm=bfield.bfield(
			unitsXB=unitsXB # setting the units
			)
		if wEnds :
			bfieldHarm=bfieldHarm.create_harm_field_with_ends(
				periodLength=self.periodLength.get(), 
				amplitude=self.bEffY.get(), 
				nperiods=self.numPeriods.get(),
				num_pnts_per_period = 200, 
				colx = 'x', 
				coly = 'By',
				shift=0,
				)
		else :
			bfieldHarm=bfieldHarm.create_harm_field(
				periodLength=self.periodLength.get(), 
				amplitude=self.bEffY.get(), 
				nperiods=self.numPeriods.get(),
				phase_shift = 0, 
				num_pnts_per_period = 200, 
				colx = 'x', 
				coly = 'By',
				deltaX=None,
				)
		return bfieldHarm

	def createEasyRepresentation(self) :
		"""
		creates a magnets list object representing the structure
		"""
		pass