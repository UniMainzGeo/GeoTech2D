**********************************
Calling scripts for GeoTech2D code
**********************************

To start simulation do the following:

- create directory mesh
- place binary .npz file in directory mesh as specified in setup script (e.g. crust.npz)
  mesh files can be generated with scrips provided in this repository
- invoke the setup script, e.g.:

python crust.py

Input parameters must be defined in a separate calling script

Parameter definition should be followed by a call to runGeoTech2D function

See description in GeoTech2D.py module

To facilitate input preparation use helper functions form utils.py module

Requirements:

conda install -c anaconda numpy scipy
conda install -c conda-forge pyevtk