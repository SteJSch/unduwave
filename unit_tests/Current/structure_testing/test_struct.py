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

class my_lib_test(unittest.TestCase) :

	def test_tensor_hdf5(self) :
		"""
		test test
		"""
		self.assertEqual( 1,1 )

if __name__ == '__main__':
	unittest.main()
