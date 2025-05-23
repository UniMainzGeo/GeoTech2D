import numpy as np
from src.local import viscoElastic, returnMap, getf, getII

#------------------------------------------------------------------------------
def stressUpdate(deps, dp, pfl, T, hist, svar, cnt, params, phase, ID, m, actnl, dt) :

    # get solution parameters
    ltol   = params['ltol']
    lmaxit = params['lmaxit']
    lstol  = params['lstol']
    amin   = params['amin']
    rhob   = phase ['rhob']

    # get material parameters
    rhob = phase['rhob']

    # get history variables
    sn  = hist['svtau'][cnt, :]
    pn  = hist['svp'  ][cnt]
    aps = hist['svaps'][cnt]
    avs = hist['svavs'][cnt]

    # compute elastic, creep and plastic parameters
    K, G, AL, AN, n, k, kf, c, ps, a, b, pd, py, pf, R, eta = getMatParam(params, phase, actnl, T, aps, dt)

    # get strain norm
    snrm = np.linalg.norm(deps, np.inf)

    # compute trial pressure (positive in compression)
    pstar = pn - dp

    # perform volumetric-deviatoric strain increment decomposition
    de  = np.dot(ID, deps)
    dth = np.dot(deps, m)

    # compute effective deviatoric strain increment
    dd = de + sn/2./G

    # get normalized deviatoric tensor
    ddII = getII(dd)

    if ddII:
        N = dd/ddII
    else:
        N = np.zeros((4))

    # compute volumetric, effective and total deviatoric strain rates
    th  = dth/dt
    dII = ddII/dt
    eII = getII(de)/dt

    # compute visco-elastic stress predictor
    sII, eta_eff, beta_1, pdiv = viscoElastic(dII, AL, AN, n, ltol, lmaxit)

    # compute yield function
    f = getf(k, c, a, pd, py, R, pstar, sII)

    # check for plasticity
    if f > 0. and actnl:

        # compute plastic stress corrector
        sII, pcor, dIIpl, thpl, eta_eff, beta_1, beta_2, beta_3, beta_4, cdiv = returnMap(dII, pstar, sII, K, AL, AN, n, k, kf, c, ps, a, b, pd, py, pf, R, eta, dt, ltol, lmaxit, lstol, amin, snrm)

        # update APS and AVS
        aps += dIIpl*dt
        avs += thpl *dt

    else:
        pcor   =  pstar
        beta_2 =  0.0
        beta_3 =  0.0
        beta_4 = -1.0
        cdiv   =  0

    # compute deviatoric stress
    tau = sII*N

    # compute effective stress (positive in extension)
    sig = tau - (pcor + pfl)*m

    # compute volumetric residual
    rth = dth - dp/K

    # compute normalized volumetric residual
    nrth = np.abs(rth)

    # compute tangent matrices
    Duu = (2.*eta_eff*ID + beta_1*np.outer(N, N) + beta_2*np.outer(m, N))/dt
    Dup = -beta_3*N - beta_4*m
    Dpu = m
    Dpp = -1./K

    # pack solution variables
    svar['svtau'][cnt, :] = tau[:]
    svar['svp'  ][cnt]    = pcor
    svar['svaps'][cnt]    = aps
    svar['svavs'][cnt]    = avs
    svar['sveII'][cnt]    = eII
    svar['svsII'][cnt]    = sII
    svar['svth' ][cnt]    = th
    
    return sig, rth, nrth, Duu, Dup, Dpu, Dpp, rhob, pdiv + cdiv

#------------------------------------------------------------------------------
def getMatParam(params, phase, actnl, T, aps, dt):

    # get material parameters  
    K      = phase['K']
    G      = phase['G']
    Bd     = phase['Bd']
    Ed     = phase['Ed']
    Bn     = phase['Bn']
    En     = phase['En']
    n      = phase['n']
    phi    = phase['phi']
    psi    = phase['psi']
    cmc    = phase['cmc']
    Hphi   = phase['Hphi']
    Hcmc   = phase['Hcmc']
    ps     = phase['ps']
    phimin = phase['phimin']
    cmcmin = phase['cmcmin']
    eta    = phase['eta']
    
    # get solution parameters
    Rugc   = params['Rugc']
    Tshift = params['Tshift']

    # shift temperature from C to K
    T = T + Tshift

    # switch off power-law for initialization
    if not actnl:
        Bn = 0.0
        En = 0.0
        n  = 1.0

    # get visco-elastic constants
    RT = Rugc*T
    AL = 1./2./G/dt + Bd*np.exp(-Ed/RT)
    AN = Bn*np.exp(-En/RT)

    # apply hardening /softening
    phi += Hphi*aps
    cmc += Hcmc*aps

    if phi < phimin:
        phi  = phimin
        # hardening modulii must be set consistently for a coupled model
        # Hphi = 0.

    if cmc < cmcmin:
        cmc  = cmcmin
        # hardening modulii must be set consistently for a coupled model
        # Hcmc = 0.

    # DP and flow potential parameters
    k   = np.sin(phi)
    kf  = np.sin(psi)
    c   = np.cos(phi)*cmc

    # cap center, radius and scaling
    a    = np.sqrt(1 + k*k)
    cosa = 1/a
    sina = k/a
    py   = (ps + c*cosa)/(1 - sina)
    R    = py - ps

    # cap-DP delimiter
    pd = py - R*sina

    # flow potential center and scaling
    pf = pd + kf*(c + k*pd)
    b  = np.sqrt(1 + kf*kf)

    # return parameters
    return K, G, AL, AN, n, k, kf, c, ps, a, b, pd, py, pf, R, eta

#------------------------------------------------------------------------------