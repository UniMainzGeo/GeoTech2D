import os
import sys
import gc
import datetime
import numpy as np
from scipy.sparse.linalg import spsolve
from pyevtk.vtk import VtkGroup
import src.scale as scl
from src.functions import getShapeQuadLine, getShapeQuadTria
from src.bc import setBCVals
from src.assemble import getElemMatVecMech, getFaceTracVecMech, assembleMech
from src.output import export2Paraview, savePVDtmp
from src.history import allocHistVar, storeHistVar
from src.setup import initT, initAPS, initPf, updatePf

#------------------------------------------------------------------------------
def runGeoTech2D(data):

#==============================================================================
# INPUT PARAMETERS
#
# data      - dictionary containing all input parameters (see below)
#
# setupname - output sub-directory name (created automatically in output directory)
# meshfile  - mesh file name expected to locate in mesh directory (should exist)
# phases    - material parameters (dictionary of dictionaries see utils.py) 
#             friction and dilatation angles must be specified in degrees
# params    - time step parameters and solution controls (dictionary, see utils.py)
#             time must be specified in years
# bc_nums   - boundary marker indices (start from 2) (list)
# bc_type   - boudary condition types "disp" or "track" (list)
#             "disp" - velocity in mm/yr, "track" - traction vector
# bc_vals   - boundary condition values (list of numpy arrays) 
# dyke      - optional dyke setup parameters (dictionary, see dyke.py)
# aps       - optional APS perturbation parameters (dictionary, see crust.py)
# temp      - optional temperature initialization parameters (dictionary, see ductile.py)
#==============================================================================

    # start time tracking   
    tb = datetime.datetime.now()
       
    # access input parameters
    setupname = data['setupname']
    meshfile  = data['meshfile']
    params    = data['params']
    phases    = data['phases']
    bc_nums   = data['bc_nums']
    bc_type   = data['bc_type']
    bc_vals   = data['bc_vals']
    dyke      = data['dyke']
    aps       = data['aps']
    temp      = data['temp']

    # get output directory and file name 
    outdir   = os.path.join("output", setupname)
    filename = setupname

    # create output directory
    if not os.path.isdir('output'):
        os.makedirs('output') 

    if not os.path.isdir(outdir):
        os.makedirs(outdir) 
        
    # get output file name
    outfile = os.path.join(outdir, filename)

    # create group output file
    vtkgrp = VtkGroup(outfile)

    # read mesh
    mesh  = np.load(os.path.join("mesh", meshfile) + '.npz')
    nc    = mesh['nc']
    nodes = mesh['nodes']
    elems = mesh['elems']
    faces = mesh['faces']
    marks = mesh['marks']
    attr  = mesh['units']
    mesh.close()

    # initialize parameters
    dtref  = params['dt']
    dtmin  = params['dtmin']
    tmax   = params['tmax']
    tinit  = params['tinit']
    nrstep = params['nrstep']
    actnl  = params['actnl']
    rtol   = params['rtol']
    dtol   = params['dtol']
    maxit  = params['maxit']

    # set constants
    ndim  = 2  # spatial dimensions
    nnel  = 7  # nodes per element
    ncel  = 3  # corner nodes per element
    nipe  = 6  # integration points per element
    nnfc  = 3  # nodes per face
    ncfc  = 2  # corner nodes per face
    nipf  = 3  # integration points per face
    nedof = 17 # degrees of freedom per element

    # assign mesh parameters
    nc   = int(nc)        # number of corner nodes
    nn   = nodes.shape[0] # number of nodes
    ne   = elems.shape[0] # number of elements
    udof = ndim*nn        # displacement dof
    pdof = ncel*ne        # pressure dof
    ndof = udof + pdof    # total dof

    # get shape functions and derivatives in quadrature points
    Nu, dNu, w  = getShapeQuadTria(nnel, nipe)
    Np, dNp, _  = getShapeQuadTria(ncel, nipe)
    NL, _,  wL  = getShapeQuadLine(nnfc, nipf)
    Nc,  dNc, _ = getShapeQuadLine(ncfc, nipf)

    # setup matrices, vectors
    sz      = ne*nedof**2
    rows    = np.zeros((sz), dtype = np.int32)
    cols    = np.zeros((sz), dtype = np.int32)
    vals    = np.zeros((sz), dtype = np.float64)
    F       = np.zeros((ndof))
    Fint    = np.zeros((ndof))
    Fext    = np.zeros((ndof))
    U       = np.zeros((ndof))
    Ut      = np.zeros((ndof))
    T       = np.zeros((nc))
    pf      = np.zeros((nc))

    dbcflags = np.zeros((ndof), dtype='bool')
    
    # allocate history variables
    hist, svar = allocHistVar(ne, nipe)

    # get boundary displacements, tractions
    idisp, itrac, vdisp, vtrac = setBCVals(faces, marks, bc_nums, bc_type, bc_vals)

    # store flags of constrained dof
    dbcflags[idisp] = True

    # initialize APS perturbation
    if aps:
        initAPS(aps['xsig'], aps['xmu'], aps['ysig'], aps['ymu'], aps['daps'], aps['size'], elems, nodes, Np, nipe, hist)

    # initialize dyke pressure source
    if dyke:
        initPf(dyke['xsig'], dyke['xmu'], dyke['ysig'], dyke['ymu'], dyke['psource'], nodes, nc, pf)

    # initialize temperature
    if temp:
        initT(temp, nodes, nc, T)

    # TIME STEP LOOP
    dt    = dtref
    time  = 0.
    times = []
    step  = 0
    init  = tinit > 0.0
    if not init:
        actnl = True
    restart = 0

    print('===============================================')
    print('Running simulation ...')
    print('Number of elements : ', ne)
    print('Number DOF         : ', ndof)

    while time < tmax:

        # increase time step after two successfull restart steps
        if dt < dtref and restart == int(nrstep):
            dt      *= 2.0
            restart  = 0

        # zero out displacement increment vector
        U[:] = 0.

        # apply displacment increments
        if not init :
            U[idisp] = vdisp*dt
        
        print('===============================================')
        print('STEP:', step+1, 'dt:', dt/scl.yr, '[yr]', 'Time:', (time+dt)/scl.yr, '[yr]')
        print('   ITER: FRES:         VRES:         ASSEMBLY:       SOLVE:')
        
        # NEWTON-RAPHSON ITERATION
        conv = False

        for it in range(int(maxit)) :

            t1 = datetime.datetime.now()

              # ELEMENT COMPUTATION & SPARSE ASSEMBLY
            vres, ndiv = getElemMatVecMech(nodes, elems, params, phases, attr, hist, svar, dbcflags, U, pf, T, Np, dNp, Nu, dNu, w, rows, cols, vals, Fint, Fext, actnl, it, dt)

            getFaceTracVecMech(nodes, itrac, vtrac, params, NL, dNc, wL, Fext)

            K, fres = assembleMech(ndof, udof, rows, cols, vals, Fint, Fext, F, dbcflags, params)

            t2 = datetime.datetime.now()

            # LINEAR SOLVE
            dU = spsolve(K, -F)

            t3 = datetime.datetime.now()

            print('   %-4d  %e  %e  %s  %s'  % (it+1, fres, vres, str(t2-t1), str(t3-t2)))

            # update displacment increment
            U += dU

            # check local iteration divergence
            if ndiv:
                print('   *** WARNING! time step diverged due to diverged local iterations')
                break

            # check displacement divergence
            if np.abs(dU[0:udof]).max() > dtol:
                print('   *** WARNING! time step diverged due to large displacement')
                print(str(round((np.abs(dU[0:udof]).max()/ dtol) *100.,2)),"% of dtol")
                break

            if fres < rtol :

                print('   TIME STEP CONVERGED!')
                conv = True
                break

            # check iteration overflow
            if it+1 == maxit:
                print('*** WARNING! time step diverged due to iteration overflow')
                break

        # check covergence
        if not conv:

            print('   *** RESTARTING TIME STEP')
            dt      /= 2.0
            restart  = 0

            if dt < dtmin:
                vtkgrp.save()
                sys.exit('   *** ERROR MINIMUM TIME STEP REACHED. STOP')
            else:
                continue

        # increment time & time step index
        time += dt
        step += 1

        # store time stamp for output
        times.append(time)

        # update restart counter
        if dt < dtref:
            restart += 1

        # update total displacements
        Ut += U

        # store history variables
        storeHistVar(hist, svar)

        # zero out displacement and APS after initialization stage
        if init and time >= tinit:

            print('===============================================')
            print('INITIALIZATION PHASE COMPLETE ...')
            print('RESETING DISPLACEMENTS ...')

            Ut[0:udof] = 0.
            init       = False
            actnl      = True

        # store time step
        export2Paraview(vtkgrp, outfile, step, nodes, elems, nc, nipe, attr, Ut, U, pf, T, time, dt, svar)

        # update pressure in the propagating dyke
        if dyke:
            updatePf(dyke['xsig'], dyke['xmu'], dyke['avsref'], elems, nodes, nipe, hist, nc, pf)

        # store temporary pvd file
        savePVDtmp(outdir, filename, step, times)

        # clear memory
        gc.collect()

    # store group output file
    vtkgrp.save()

    te = datetime.datetime.now()

    print('#===========================================')
    print('Total solution time         :', str(te-tb))
    print('Total number DOF            :', ndof)
    print('Total time reached          :', time/scl.yr, '[yr]')
    print('#===========================================')

#------------------------------------------------------------------------------
