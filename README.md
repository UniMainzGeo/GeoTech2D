# GeoTech2D

This repository contains the source code of `GeoTech2D` which is a 2D unstructured finite element code for geomechanical applications with visco-elasto-viscoplastic material properties and mode-I and mode-II regularized plasticity. It employs LBB-stable conforming Crouzeix-Raviart triangular elements using a mixed displacement and pressure increment formulation.

The code and the numerical formulation to solve problems with combined mode-I/mode-II plasticity models in a robust manner is described in detail in:

- A.A. Popov, N. Berlie and B.J.P. Kaus (2025, submitted) A dilatant visco-elasto-viscoplasticity model with globally continuous tensile cap: stable two-field mixed formulation. *Geoscientific Model Development*


![Mode-I propagation](/VIDEO/dyke_Vx_50.gif)
![Mode-I propagation](/VIDEO/dyke_Pf_50.gif)
![Crustal scale extension wioth mode-I & mode-II plasticity](/VIDEO/ductile_EII_50.gif)


### Instructions to run the GeoTech2D code

To start a simulation do the following:

- create a directory called `mesh`
- place a binary `.npz` file in directory `mesh` as specified in setup script (e.g. `crust.npz`)
- invoke the setup script from python (e.g. `python crust.py`)

Input parameters must be defined in a separate calling script

Parameter definition should be followed by a call to the `runGeoTech2D` function.
See description in `CODE/src/GeoTech2D.py` module

To facilitate input preparation use the helper functions from the `CODE/src/utils.py` module

#### Requirements

You need python with the following packages:
```
conda install -c anaconda numpy scipy

conda install -c conda-forge pyevtk
```

### Instructions to generate mesh with MeshPy package

`The` MeshPy package implements the Python API for the `Triangle` quality mesh generator

To generate a mesh simply invoke the corresponding setup script (e.g.` python crust.py`)

Binary .npz file will be placed in directory `mesh` (will be created if necessary)

Input parameters must be defined in a separate calling script

Parameter definition should be followed by a call to `runMeshPy` function

See description in `MESH/src/meshpy_triangle_api.py` module

To facilitate input preparation use helper functions from the `MESH/src/utils.py` module

#### Requirements
```
conda install -c anaconda numpy scipy matplotlib

conda install -c conda-forge meshpy distinctipy
```
