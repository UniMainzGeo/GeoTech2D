import numpy as np
from src.utils import initParams, initPhase, scaleParams
from src.GeoTech2D import runGeoTech2D

#------------------------------------------------------------------------------
def main():
    
    # set setup name
    setupname = 'tensile'
    meshfile  = 'tensile'

    # solution parameters
    params           = initParams()
    params['dt']     = 100.0 
    params['tmax']   = 100000.0
    params['dtmin']  = 0.05
    params['tinit']  = 100.0
    params['nrstep'] = 20
    params['ltol']   = 1e-8
    params['lmaxit'] = 100
    params['lstol']  = 0.9
    params['amin']   = 0.5
    params['rtol']   = 1e-5
    params['dtol']   = 1e6
    params['maxit']  = 10

    # phases
    phases          =  { }
    phase           =  initPhase()
    phase['K']      =  6.4e10
    phase['G']      =  4e10
    phase['rhob']   =  3000.0 
    phase['phi']    =  30.0 
    phase['cmc']    =  2e7
    phase['Hcmc']   = -1e8
    phase['cmcmin'] =  5e6
    phase['ps']     = -1e6
    phase['eta']    =  1e19
    
    phases[1] = phase
    
    phase           =  initPhase()
    phase['K']      =  6.4e10
    phase['G']      =  4e9
    phase['rhob']   =  3000.0 
    phase['phi']    =  30.0 
    phase['cmc']    =  2e7
    phase['Hcmc']   = -1e8
    phase['cmcmin'] =  5e6
    phase['ps']     = -1e6
    phase['eta']    =  1e19
    
    phases[2] = phase

    # boundary conditions 
    bc_nums = [2, 3, 4, 5]
    bc_type = ["disp", "disp", "disp", "disp"]
    bc_vals = [
        np.array([-1e-4,    np.nan]),
        np.array([ 1e-4,    np.nan]),
        np.array([ np.nan,  0.0]),
        np.array([ np.nan,  0.0])]

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
    data['aps']       = { }    
    data['dyke']      = { }   
    data['temp']      = { }
    
    runGeoTech2D(data)

#------------------------------------------------------------------------------
if __name__ == "__main__":

    main()

#------------------------------------------------------------------------------
