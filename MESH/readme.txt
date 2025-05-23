**********************************
Calling scripts for MeshPy package
**********************************

Python API for Triangle quality mesh generator

To generate mesh simply invoke corresponding setup script, e.g.:

python crust.py

Binary .npz file will be placed in directory mesh (will be created if necessary)

Input parameters must be defined in a separate calling script

Parameter definition should be followed by a call to runMeshPy function

See description in meshpy_triangle_api.py module

To facilitate input preparation use helper functions form utils.py module

Requirements:

conda install -c anaconda numpy scipy matplotlib
conda install -c conda-forge meshpy distinctipy