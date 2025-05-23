import numpy as np
import src.scale as scl

#------------------------------------------------------------------------------
def initParams() :
    
    params           = {}
    params['g']      = 9.81    # [m/s^2]   gravity acceleration
    params['dt']     = 1.0     # [yr]      time step
    params['dtmin']  = 1.0     # [yr]      min. time step
    params['tmax']   = 1.0     # [yr]      simulation time
    params['tinit']  = 0.0     # [yr]      initialization time
    params['nrstep'] = 5       # [ ]       number of constant steps after restart
    params['ltol']   = 1e-9    # [ ]       local iteration tolerance
    params['lmaxit'] = 20      # [ ]       max. num. local iterations
    params['lstol']  = 0.95    # [ ]       line search tolerance
    params['amin']   = 0.05    # [ ]       line search minimum step length
    params['rtol']   = 1e-6    # [ ]       relative tolerance
    params['dtol']   = 0.5     # [ ]       displacement tolerance
    params['maxit']  = 20      # [ ]       max. num. iterations
    params['Rugc']   = 8.31446 # [J/mol/K] universal gas constant
    params['Tshift'] = 273.15  # [K]       temperature shift
    params['actnl']  = False   # [ ]       nonlinearity at initialization stage activation flag

    return params

#------------------------------------------------------------------------------
def initPhase() :
    
    phase           = {}
    phase['K']      = 0.0 # [Pa]       bulk modulus
    phase['G']      = 0.0 # [Pa]       shear modulus
    phase['rhob']   = 0.0 # [kg/m^3]   bulk density
    phase['phi']    = 0.0 # [deg]      friction angle
    phase['psi']    = 0.0 # [deg]      dilatation angle
    phase['cmc']    = 0.0 # [Pa]       Mohr-Coulomb cohesion
    phase['Hphi']   = 0.0 # [deg]      friction hardening modulus
    phase['Hcmc']   = 0.0 # [Pa]       cohesion hardening modulus
    phase['phimin'] = 0.0 # [deg]      min. friction angle
    phase['cmcmin'] = 0.0 # [Pa]       min cohesion
    phase['ps']     = 0.0 # [Pa]       tensile strength (pressure, must be negative)
    phase['eta']    = 0.0 # [Pa*s]     stabilization viscosity
    phase['Bd']     = 0.0 # [1/Pa/s]   linear creep pre-factor
    phase['Ed']     = 0.0 # [J/mol]    linear creep activation energy
    phase['Bn']     = 0.0 # [1/Pa^n/s] power-law creep pre-factor
    phase['En']     = 0.0 # [J/mol]    power-law creep activation energy
    phase['n']      = 1.0 # [ ]        power-law creep exponent

    return phase

#------------------------------------------------------------------------------
def scaleParams(params, phases, bc_type, bc_vals) :

    # switch from years to seconds
    params['dt']    = scl.yr*params['dt']
    params['tmax']  = scl.yr*params['tmax']
    params['tinit'] = scl.yr*params['tinit']
    params['dtmin'] = scl.yr*params['dtmin']

    for phase in phases.values():
    
        # switch from degrees to radians
        phase['phi']    = np.radians(phase['phi'])
        phase['psi']    = np.radians(phase['psi'])
        phase['Hphi']   = np.radians(phase['Hphi'])
        phase['phimin'] = np.radians(phase['phimin'])

    for i, tp in enumerate(bc_type):
        
        # scale boundary velocity 
        if tp == 'disp':
            bc_vals[i] *= scl.mmperyr

#------------------------------------------------------------------------------