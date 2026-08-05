# unduwave

## What is it Good For?

unduwave is a python library that allows for magnetic field calculations of very general magnet-arrangements and the calculation of the properties of the radiation produced by relativistic electron beams traversing those magnetic fields. This allows, among other things, for design and analysis of synchrotrons and other particle accelerator functional elements - undulators, wave-length shifters, dipoles, sextupoles, etc. The magnetic geometries can be constructed from different materials and elements like magnetblocks (general convex geometries). Helper functions are offered for easy construction of different kinds of undulators - planar, elliptical, hybrid - from basic building blocks. 

Synchrotron radiation spectra can be calculated on screens positioned downstream from radiation-producing sources (magnetic fields through which relativistic energy beams traverse). Fluxes, flux densities, stokes parameters, powr distributions and more are easily calculable. The calculations can take the beam characteristics like energy-spread and emittance into account. 

Under the hood unduwave uses the fortran programs WAVE and Undumag developed by Michael Scheer. See below for a short explanation.

## How to install unduwave

To install this package run:

`python -m pip install unduwave`

## Get started using unduwave

To get started, see the documentation on [github pages](https://stejsch.github.io/unduwave), or on [readthedocs](https://unduwave.readthedocs.io/en/latest/).

## How to cite unduwave

[1] S. Schäfer. M. Scheer. (2026). unduwave. [Online]. Available: https://github.com/SteJSch/unduwave

## Version

- Current Version: 1.0.0

## License

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
