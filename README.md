# GeoTech2D
[![DOI](https://zenodo.org/badge/984744668.svg)](https://doi.org/10.5281/zenodo.15496842)

This repository contains the source code of `GeoTech2D` which is a 2D unstructured finite element code for geomechanical applications with visco-elasto-viscoplastic material properties and mode-I and mode-II regularized plasticity. It employs LBB-stable conforming Crouzeix-Raviart triangular elements using a mixed displacement and pressure increment formulation.

The code and the numerical formulation to solve problems with combined mode-I/mode-II plasticity models in a robust manner is described in detail in:

- [A.A. Popov, N. Berlie and B.J.P. Kaus (2025, submitted) A dilatant visco-elasto-viscoplasticity model with globally continuous tensile cap: stable two-field mixed formulation. *Geoscientific Model Development*](https://egusphere.copernicus.org/preprints/2025/egusphere-2025-2469/)

![Mode-I propagation](/VIDEO/dyke_Vx_50.gif)
![Mode-I propagation](/VIDEO/dyke_Pf_50.gif)
![Crustal scale extension wioth mode-I & mode-II plasticity](/VIDEO/ductile_EII_50.gif)

### 1. Instructions to run GeoTech2D

To start a simulation do the following:

- Create a directory called `mesh`
- Create a mesh as described below
- Place a binary `.npz` file in directory `mesh` as specified in the setup script (e.g. `crust.npz`)
- Invoke the setup script from python (e.g. `python crust.py`)

#### 1.1 Prepare new simulation

Define all input parameters in a separate calling script (use the supplied scripts for reference).

Parameter definition should be followed by a call to the `runGeoTech2D` function.

See description in `CODE/src/GeoTech2D.py` module.

To facilitate input preparation, you can use helper functions from the `CODE/src/utils.py` module.

#### 1.2 Requirements

```
conda install -c anaconda numpy scipy

conda install -c conda-forge pyevtk
```

### 2. Instructions to generate mesh with MeshPy package

The `MeshPy` package implements the Python API for the `Triangle` quality mesh generator.

To generate a mesh, simply invoke the corresponding setup script (e.g.` python crust.py`).

Binary `.npz` file will be placed in directory `mesh` (will be created if necessary).

#### 2.1 Generate new grid

Define all input parameters in a separate calling script (use the supplied scripts for reference).

Parameter definition should be followed by a call to `runMeshPy` function.

See description in `MESH/src/meshpy_triangle_api.py` module.

To facilitate input preparation use helper functions from the `MESH/src/utils.py` module.

#### 2.2 Requirements
```
conda install -c anaconda numpy scipy matplotlib

conda install -c conda-forge meshpy distinctipy
```
