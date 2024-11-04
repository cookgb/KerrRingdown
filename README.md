# KerrRingdown

This repository contains a *Mathematica paclet* called __KerrRingdown\`__ for fitting a Numerical Relativity ring-down signal as a linear combination of Quasinormal modes.
The various fitting methods implemented in this package are based on the methods outlined in [Cook (2020)](https://link.aps.org/doi/10.1103/PhysRevD.102.024027) and [Gao and Cook (2024)]().  Please see these reference for details and cite these references when publishing results based on these packages.

---
## Paclet

#### KerrRingdown\`
The actual *Mathematica paclet* is found in the KerrRingdown directory of the top level of the repository.  This paclet supports the following main operations:
1. Reading in Numerical Relativity gravitational waveforms.
2. Loading in tables of Quasinormal Modes.
3. Performing linear and non-linear fitting using:
  + The Eigenvalue Methods (both standard and simulation subspace limited).
  + Linear Least-Squares based on SVD fitting on the Design Matrix or the Normal Equation.
4. Plotting the signal and fit functions.
5. Plotting the Quasinormal Mode expansion coefficients.

---
## Build and Installing the Paclet

The file __CreatePaclets.nb__ is included in the top level of the repository.  It can be used to build the paclet from source.  Public releases will also be available containing pre-built versions of the KerrRingdown\` paclet which can be installed directly into your *Mathematica* installation.

Paclet installation is accomplished by using the *Mathematica* function __PacletInstall__.

___
## Tables of Quasinormal Modes.

The __KerrRingdown\`__ paclet requires a set of tabulated quasinormal modes stored in a set of HDF5 files.  These files are not included in the Github repository, but can be found on [Zenodo](https://doi.org/10.5281/zenodo.14024959 "Quasinormal Mode Repo").  Note that the current version of __KerrRingdown\`__ requires the files from __*Version V3*__.  Prior versions use a data format that will not work correctly.  The tabulated quasinormal modes are stores in files named __KerrQNM___*nn*__.h5__, where *nn* is a 2-digit number indicating the value of the overtones for the modes stored in each file.

---
## Mathematica-Style Documentation

Full Mathematica-style documentation is included with the paclet.  This includes detailed documentation for each publicly callable function, a *Function Guide*, and a *Tech Note*.  The *Tech Note* includes details on how fitting works with the paclet.  To view the documentation, first use __PacletInstall__ to intall the __KerrRingdown\`__ paclet, then execute the command __Needs["KerrRingdown\`"]__ to load the functions and their documentation.  At this point, you can search the *Mathematica* documentation for *Kerr Ringdown*, or you could enter the name of one of its functions to access its documentation.  For example, the function __HDF5QNMDir__ is used to specify the directory where the tabulated quasinormal mode data files are found, and the documentation for this function will also lead you to the *Function Guide* and a *Tech Note*.
