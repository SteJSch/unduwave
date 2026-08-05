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

- Current Version: 0.9.2
