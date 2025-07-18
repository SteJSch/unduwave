"""
Contains the import statements for undupy modules
"""

def to_scn(number : float,norm: bool=True):
	"""
	Converts number to string using scientific notation

	:param number: The number to convert
	:param norm: If true, std scientific notation, False: With leading 0.
	:return: The resulting string.
	"""
	if norm : 
		return '{:.5E}'.format(number)
	else :
		a, b = '{:.4E}'.format(number).split('E')
		return '{:.5f}E{:+03d}'.format(float(a)/10, int(b)+1)

import os
import pandas as pd
import h5py
import pdb
import numbers
import copy
import ast
import math
import numpy as np
import scipy as scipy
import scipy.special as ss
import matplotlib.pyplot as plt
# import apy.src.helpers.physical_constants as pc
# from apy.src.helpers import plot_helpers as plt_h
import sys
import subprocess
import shutil
import random
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
from datetime import datetime
from scipy import interpolate
from scipy import integrate
from scipy.integrate import nquad
from scipy.interpolate import RectBivariateSpline
from scipy.signal import find_peaks
from matplotlib import cm
from pathlib import Path
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
import time
random.seed(datetime.now().timestamp())
cm_inch = 1/2.54
"""
Conversion factor cm to inch
"""
