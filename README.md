# GeoTech2D

This repository contains the source code of `GeoTech2D` which is a 2D unstructured finite element code for geomechanical applications with visco-elasto-viscoplastic material properties and mode-I and mode-II regularized plasticity. It employs LBB-stable conforming Crouzeix-Raviart triangular elements using a mixed displacement and pressure increment formulation.

The code and the numerical formulation to solve problems with combined mode-I/mode-II plasticity models in a robust manner is described in detail in:

- A.A. Popov, N. Berlie and B.J.P. Kaus (2025, submitted) A dilatant visco-elasto-viscoplasticity model with globally continuous tensile cap: stable two-field mixed formulation. *Geoscientific Model Development*

***
### Instructions to run GeoTech2D code

To start simulation do the following:

- create directory mesh
- place binary .npz file in directory mesh as specified in setup script (e.g. crust.npz)
- invoke the setup script (e.g. python crust.py)

Input parameters must be defined in a separate calling script

Parameter definition should be followed by a call to runGeoTech2D function

See description in GeoTech2D.py module

To facilitate input preparation use helper functions form utils.py module

#### Requirements

conda install -c anaconda numpy scipy

conda install -c conda-forge pyevtk

***
### Instructions to generate mesh with MeshPy package

MeshPy package implements Python API for Triangle quality mesh generator

To generate mesh simply invoke corresponding setup script (e.g. python crust.py)

Binary .npz file will be placed in directory mesh (will be created if necessary)

Input parameters must be defined in a separate calling script

Parameter definition should be followed by a call to runMeshPy function

See description in meshpy_triangle_api.py module

To facilitate input preparation use helper functions form utils.py module

#### Requirements

conda install -c anaconda numpy scipy matplotlib

conda install -c conda-forge meshpy distinctipy
