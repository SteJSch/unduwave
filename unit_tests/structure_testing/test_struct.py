"""
Contains the import statements for wavepy modules
"""

import unittest
import math
import pathlib
import sys 
import numpy as np
import pdb
import os
import copy 
import random
import h5py

sys.path.insert(0,'../../../')
import unduwave as uw

class tA(uw.attribute_classes.attributes._attribute_collection) :
	varA=4
	varC= uw.attribute_classes.attributes._attribute(value=1,name='varC')
	def __init__(self) :
		self.varB=5

class my_lib_test(unittest.TestCase) :

	def test_class_members(self) :
		"""
		test test
		"""
		self.tA=tA()
		self.tA.varA=5
		self.tA.varC.set(4)
		self.tA2=tA()
		self.assertEqual( self.tA.varA,5 )
		self.assertEqual( self.tA2.varA,4 )
		# self.assertEqual( self.tA.varC.get(),4 )
		# self.assertEqual( self.tA2.varC.get(),1 )


	def test_class_members2(self) :
		eBeamPara=uw.ebeam_parameters()
		attributes = vars(eBeamPara.__class__)
		# for attr_name, attr_value in dict(eBeamPara).items():
		# for attr_name, attr_value in attributes.items():
			# if isinstance(attr_value, uw.attribute_classes.attributes._attribute):
			# 	print("Found")
			# 	pdb.set_trace()

	def test_tensors(self) :
		numX=3
		minX=-3
		maxX=3
		numY=3
		minY=2
		maxY=3
		numZ=3
		minZ=-4
		maxZ=-5

		xCoords=np.linspace(minX,maxX,numX)
		yCoords=np.linspace(minY,maxY,numY)
		zCoords=np.linspace(minZ,maxZ,numZ)

		gridX, gridY, gridZ = np.meshgrid(xCoords,yCoords,zCoords,indexing='ij')

		pdb.set_trace()

if __name__ == '__main__':
	unittest.main()
