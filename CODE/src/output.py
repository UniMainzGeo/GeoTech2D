import os
import numpy as np
from pyevtk.hl import unstructuredGridToVTK
from pyevtk.vtk import VtkTriangle, VtkGroup
import src.scale as scl
from src.history import outputSolVar

#------------------------------------------------------------------------------
def export2Paraview(vtkgrp, outfile, step, nodes, elems, nc, nipe, attr, Ut, U, pf, T, time, dt, svar):

    # get number of elements, coordinates, connectivity and displacements
    ne   = elems.shape[0]
    x    = np.copy(nodes[0:nc, 0])
    y    = np.copy(nodes[0:nc, 1])
    z    = np.zeros((nc))
    conn = np.copy(elems[:,  0:3]).flatten().astype(int)
    Ux   = np.copy(Ut[range(0, 2*nc, 2)])
    Uy   = np.copy(Ut[range(1, 2*nc, 2)])
    Vx   = np.copy(U [range(0, 2*nc, 2)]/scl.mm)/(dt/scl.yr)
    Vy   = np.copy(U [range(1, 2*nc, 2)]/scl.mm)/(dt/scl.yr)

    # define offset of last vertex of each element
    offset = np.cumsum(3*np.ones(ne, dtype=int))

    # define cell types
    ctype    = np.zeros(ne, dtype=int)
    ctype[:] = VtkTriangle.tid

    cellData  = {}
    pointData = {}

    # store vectors
    cellData ["Units[     ]"] = np.copy(attr)
    pointData["Disp [m    ]"] = (Ux,  Uy,  z)
    pointData["Vel  [mm/yr]"] = (Vx,  Vy,  z)
    pointData["Pf   [MPa  ]"] = np.copy(pf/scl.MPa)
    pointData["Temp [C    ]"] = np.copy(T)

    # output solutin variables
    outputSolVar(elems, nc, nipe, pointData, svar)

    # get time step output file name
    fname = outfile + str(step).zfill(5)

    # write data to file
    unstructuredGridToVTK(fname, x+Ux, y+Uy, z, connectivity = conn.flatten(), offsets = offset, cell_types = ctype, cellData = cellData, pointData = pointData)

    vtkgrp.addFile(filepath = fname + '.vtu', sim_time = time/scl.yr)

#------------------------------------------------------------------------------
def savePVDtmp(outdir, pfile, step, time):

    # save pvd of current state of simulation
    vtkgrp = VtkGroup(os.path.join(outdir, 'TEMP-' + pfile))

    for istep in range(step):

        # get time step output file name
        fname = os.path.join(outdir, pfile) + str(istep+1).zfill(5)

        vtkgrp.addFile(filepath = fname + '.vtu', sim_time = time[istep]/scl.yr)

    vtkgrp.save()

#------------------------------------------------------------------------------
