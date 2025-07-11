"""
Contains physical constants for unduwave
"""

import math

fein_const=1/137 # Fein Structure Constant
hbar=1.054e-34 # Planck's constant devided by 2 pi
q_el = 1.602e-19 # Charge of electron, [C]
m_el = 9.109e-31 # mass pf electron [kg]
v_c = 299792458 # velocity of light [m/s]
mu0=8.854e-12 # Vacuum permittivity [F/m]
R_el = q_el**2/(4*math.pi*mu0*m_el*v_c**2) # classical electron radius [m]

