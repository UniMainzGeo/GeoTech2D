import numpy as np
from src.utils import initParams, initPhase, scaleParams
from src.GeoTech2D import runGeoTech2D

#------------------------------------------------------------------------------
def main():
    
    # set setup name
    setupname = 'ductile'
    meshfile  = 'ductile'

    # solution parameters
    params           = initParams()
    params['dt']     = 1e3 
    params['tmax']   = 35e3
    params['dtmin']  = 1e-2
    params['tinit']  = 1e3
    params['nrstep'] = 15
    params['ltol']   = 1e-8
    params['lmaxit'] = 100
    params['lstol']  = 0.9
    params['amin']   = 0.5
    params['rtol']   = 1e-4
    params['dtol']   = 1e6
    params['maxit']  = 15

    # phases
    phases          =  { }
    phase           =  initPhase()
    phase['K']      =  1.1e+11
    phase['G']      =  5e10
    phase['rhob']   =  3000.0
    phase['Bd']     =  5.000000e-24 
    phase['Bn']     =  8.897100e-25
    phase['n']      =  3.3
    phase['En']     =  1.900000e+05
    phase['phi']    =  30.0
    phase['psi']    =  3.0 
    phase['cmc']    =  2e7
    phase['Hcmc']   = -5e8
    phase['cmcmin'] =  5e6
    phase['ps']     = -1e6
    phase['eta']    =  1e19

    phases[1] = phase

    # boundary conditions 
    bc_nums = [2, 3, 4]
    bc_type = ["disp", "disp", "disp"]
    bc_vals = [
        np.array([-1.577880,  np.nan]),
        np.array([ 1.577880,  np.nan]),
        np.array([np.nan, 0.0])]

    #  APS perturbations
    aps         = { }
    aps['xsig'] = 7000.0
    aps['xmu']  = 50000.0
    aps['ysig'] = 1000.0
    aps['ymu']  = 25000.0 
    aps['daps'] = 0.005
    aps['size'] = 0.01

    # initial temperature
    temp         = { }
    temp['top']  = 5.0
    temp['grad'] = 20.0
    
    # scale parameters 
    scaleParams(params, phases, bc_type, bc_vals)

    # store input data
    data              = { }
    data['setupname'] = setupname
    data['meshfile']  = meshfile
    data['params']    = params 
    data['phases']    = phases 
    data['bc_nums']   = bc_nums
    data['bc_type']   = bc_type
    data['bc_vals']   = bc_vals
    data['aps']       = aps    
    data['dyke']      = { }   
    data['temp']      = temp
    
    runGeoTech2D(data)

#------------------------------------------------------------------------------
if __name__ == "__main__":

    main()

#------------------------------------------------------------------------------
