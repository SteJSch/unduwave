
import os
import pdb
import sys

def get_std_bessy_III_paras(ebeam_paras) :

	ebeam_paras.beam_en.set(2.5) # Beam energy in [GeV]
	ebeam_paras.current.set(0.3) # current in [A]
	ebeam_paras.bsigz.set(275e-6) # horizontal beam size [m]
	ebeam_paras.bsigzp.set(28.1e-6) # Horizontal beam divergence [rad]
	ebeam_paras.bsigy.set(22.5e-6) # vertical beam size [m]
	ebeam_paras.bsigyp.set(6.8e-6) # vertical beam divergence [rad]
	ebeam_paras.espread.set(0.963e-3) # energy spread [%]
	ebeam_paras.emitt_h.set(9.804e-11) #  horizontal  emittance [mrad]
	ebeam_paras.emitt_v.set(1.961e-12) #  vertical emittance [mrad]
	ebeam_paras.betfunh.set(0) #  horizontal beta functions [m]
	ebeam_paras.betfunv.set(0) #  vertical beta functions [m]
	ebeam_paras.circumference.set(350) # Ring circumference in [m]
	ebeam_paras.rdipol.set(2.78) # Bending radius of dipoles [m]
	return ebeam_paras
