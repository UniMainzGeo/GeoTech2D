import numpy as np
from src.utils import initParams, initPhase, scaleParams
from src.GeoTech2D import runGeoTech2D

#------------------------------------------------------------------------------
def main():
    
    # set setup name
    setupname = 'crust'
    meshfile  = 'crust'

    # solution parameters
    params           = initParams()
    params['dt']     = 50.0 
    params['tmax']   = 3000.0
    params['dtmin']  = 0.05
    params['tinit']  = 50.0
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
    phase['eta']    =  1e18

    phases[1] = phase

    # boundary conditions 
    bc_nums = [2, 3, 4]
    bc_type = ["disp", "disp", "disp"]
    bc_vals = [
        np.array([-10.0,  np.nan]),
        np.array([ 10.0,  np.nan]),
        np.array([np.nan, 0.0])]

    #  APS perturbations
    aps         = { }
    aps['xsig'] = 7000.0
    aps['xmu']  = 20000.0
    aps['ysig'] = 1000.0
    aps['ymu']  = 7000.0
    aps['daps'] = 0.005
    aps['size'] = 0.0

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
    data['temp']      = { }
    
    runGeoTech2D(data)

#------------------------------------------------------------------------------
if __name__ == "__main__":

    main()

#------------------------------------------------------------------------------
