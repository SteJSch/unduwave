"""
Contains functionality for handling files.

Module to handle various file operations including finding, moving, copying, and deleting files.
"""

from unduwave.unduwave_incl import *
from unduwave.attribute_classes.attributes import _attribute
from unduwave.attribute_classes.attributes import _attribute_collection
import unduwave.helpers.file_folder_helpers as f_h
import unduwave.constants as uc

class ana_undulator(_attribute_collection) :

	bEff = _attribute(1)
	bEff2 = _attribute(0)
	periodLength = _attribute(0.02)
	numPeriods = _attribute(78)
	lengthEndPeriodsRelative = _attribute(2*1.5)

	undulatorLength = _attribute()
	undulatorWaveNum = _attribute()
	undulatorParameterK = _attribute()
	reducedBeta = _attribute()

	undulatorFrequency = _attribute()
	undulatorPeriodTime = _attribute()
	reducedGamma = _attribute()
	reducedUndulatorParameterK = _attribute()
	maxDeflection = _attribute()
	maxAngle = _attribute()
	naturalOpeningAngle = _attribute()
	properFrequencyWeak = _attribute()
	properFrequencyStrong = _attribute()
	firstHarmEnergyWeak = _attribute()
	firstHarmEnergyStrong = _attribute()
	thetaObservation = _attribute()

	def __init__(self,
			bEff,
			periodLength,
			numPeriods,
			bEff2=0.0,
			lengthEndPeriodsRelative=0.0,
			ebeam=None,
			thetaObservation=0.0,
			unduK=None
			) :
		self.bEff.set(bEff)
		self.periodLength.set(periodLength)
		self.numPeriods.set(numPeriods)
		self.lengthEndPeriodsRelative.set(lengthEndPeriodsRelative)
		self.thetaObservation.set(thetaObservation)
		self.undulatorParameterK.set(unduK)
		self.update_non_beam_values()
		if not (ebeam is None) :
			self.update_values(ebeam=ebeam,thetaObservation=self.thetaObservation.get())

	def update_non_beam_values(self) :
		self.undulatorLength.set( self.periodLength.get()*(self.numPeriods.get()+self.lengthEndPeriodsRelative.get()) )
		self.undulatorWaveNum.set( 2 * math.pi / ( self.periodLength.get() ) )
		if not(self.undulatorParameterK.get() is None) :
			self.bEff.set(self.undulatorParameterK.get()*uc.m_el * uc.v_c * self.undulatorWaveNum.get()/uc.q_el)
		elif (self.bEff.get() is None) :
			print("Warning: ana_undulator:update_values - No B and no K given")
			return
		else :
			k1=uc.q_el * self.bEff.get() / ( uc.m_el * uc.v_c * self.undulatorWaveNum.get() )
			if not (self.bEff2.get() is None) :
				k2=uc.q_el * self.bEff2.get() / ( uc.m_el * uc.v_c * self.undulatorWaveNum.get() )
				k1=math.sqrt(k1**2+k2**2)
			self.undulatorParameterK.set(k1)

	def update_values(self,ebeam,thetaObservation=0.0) :
		self.update_non_beam_values()
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
