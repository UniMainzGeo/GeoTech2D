import numpy as np

#------------------------------------------------------------------------------
# LOCAL STRESS UPDATE FUNCTIONS
#------------------------------------------------------------------------------
def viscoElastic(dII, AL, AN, n, ltol, lmaxit) :

    #=========================================
    # NONLINEAR VISCO-ELASTIC STRESS PREDICTOR
    #=========================================

    # set divergence flag
    div = 0

    if dII != 0. and AN != 0.:

        #============================
        # non-linear non-trivial case
        #============================

        # get initial guess
        sL  =  dII/AL
        sN  = (dII/AN)**(1./n)
        sII = 1./(1./sL + 1./sN)

        # stress iteration
        for i in range(lmaxit) :

            # get residual
            r = dII - AL*sII - AN*sII**n

            # get Jacobian
            J = -AL - n*AN*sII**(n-1.)

            # update stress
            sII -= r/J

            # check convergence
            rnorm = np.abs(r/dII)

            if rnorm <= ltol:
                break

        # check diverged iteration
        if rnorm > ltol :
            div = 1

        # get effective stiffness constants
        A       = -1./J
        eta_eff =  sII/dII/2.
        beta_1  =  A/2. - eta_eff

    else:

        #====================
        # linear/trivial case
        #====================

        # get linear solution
        sII     = dII/AL
        eta_eff = 1./AL/2.
        beta_1  = 0.

    return sII, eta_eff, beta_1, div

#------------------------------------------------------------------------------
def lineSearch(x, dx, r, rnorm, rm, dII, pstar, AL, AN, n, K, k, kf, c, a, b, pd, py, pf, R, eta, dt, lstol, amin, snrm):

    # initialize step length
    alpha = 1.0

    # iterate unless step length becomes too small
    while alpha > amin:

        # apply scaled update
        xm = x + alpha*dx

        # get updated residual
        rmnorm = getRes(xm, rm, dII, pstar, AL, AN, n, K, k, kf, c, a, b, pd, py, pf, R, eta, dt, snrm)

        # check whether residual is sufficiently reduced
        if rmnorm <= lstol*rnorm:
            break

        # bisect step length
        alpha /= 2.0

    return alpha

#------------------------------------------------------------------------------
def returnMap(dII, pstar, sII, K, AL, AN, n, k, kf, c, ps, a, b, pd, py, pf, R, eta, dt, ltol, lmaxit, lstol, amin, snrm) :

    #====================================================
    # PLASTIC STRESS CORRECTOR (RETURN MAPPING ALGORITHM)
    #====================================================

    # WARNING! this function must be called for plastic cases only (check beforehand)

    x  = np.zeros((3))
    r  = np.zeros((3))
    rm = np.zeros((3))
    J  = np.zeros((3, 3))

    # initial guess
    x[0] = sII
    x[1] = pstar
    x[2] = 0.0

    # set divergence flag
    div = 0

    # stress iteration
    for i in range(lmaxit) :

        # get residual
        rnorm = getRes(x, r, dII, pstar, AL, AN, n, K, k, kf, c, a, b, pd, py, pf, R, eta, dt, snrm)

        # check convergence
        if rnorm <= ltol:
            break

        # get Jacobian
        getJac(x, J, AL, AN, n, K, k, kf, c, a, b, pd, py, pf, eta, dt)

        # get update vector
        dx = -np.linalg.solve(J, r)

        # get update length
        alpha = lineSearch(x, dx, r, rnorm, rm, dII, pstar, AL, AN, n, K, k, kf, c, a, b, pd, py, pf, R, eta, dt, lstol, amin, snrm)

        # update variables
        x += alpha*dx

    # get Jacobian
    getJac(x, J, AL, AN, n, K, k, kf, c, a, b, pd, py, pf, eta, dt)

    # compute plastic strain rates and effective stiffness constants
    sII, p, dIIpl, thpl, eta_eff, beta_1, beta_2, beta_3, beta_4 = getResults(x, J, dII, K, dt, k, kf, c, b, pd, pf)

    # check diverged iteration
    if rnorm > ltol :
        div = 1

    return sII, p, dIIpl, thpl, eta_eff, beta_1, beta_2, beta_3, beta_4, div

#------------------------------------------------------------------------------
def getII(s) :

    # get square root of second invariant of a deviatoric tensor

    sxx = s[0]
    syy = s[1]
    szz = s[2]
    sxy = s[3]
    J2  = 1./2.*(sxx**2 + syy**2 + szz**2) + sxy**2

    return np.sqrt(J2)

#------------------------------------------------------------------------------
def getf(k, c, a, pd, py, R, p, s) :

    # delimiter point
    sd = c + k*pd

    # yield function
    if s*(py - pd) >= (py - p)*sd:

        f =  s - k*p - c

    else:

        Ry = np.sqrt(s**2 + (p - py)**2)
        f  = a*(Ry - R)

    return f

#------------------------------------------------------------------------------
def getRes(x, r, dII, pstar, AL, AN, n, K, k, kf, c, a, b, pd, py, pf, R, eta, dt, snrm):

    # get deviatoric stress, pressure and plastic multiplyer
    s   = x[0]
    p   = x[1]
    lam = x[2]

    # cutoff negative deviatoric stresses
    if s < 0.0:
        s = 0.0

    # delimiter point
    sd = c + k*pd

    # yield function
    if s*(py - pd) >= (py - p)*sd:

        f = s - k*p - c

    else:

        Ry = np.sqrt(s**2 + (p - py)**2)
        f  = a*(Ry - R)

    # flow potential
    if s*(pf - pd) >= (pf - p)*sd:

        gs = 1./2.
        gp = kf

    else:

        Rf =  np.sqrt(s**2 + (p - pf)**2)
        gs =  1./2.*b/Rf*s
        gp = -b/Rf*(p-pf)

    # get residual
    r[0] = dII - AL*s - AN*s**n - lam*gs
    r[1] = (p - pstar)/K/dt     - lam*gp
    r[2] = f                    - lam*eta

    # get normalized residual
    rnorm = np.abs(r[0])/snrm + np.abs(r[1])/snrm + np.abs(r[2])/c

    return rnorm

#------------------------------------------------------------------------------
def getJac(x, J, AL, AN, n, K, k, kf, c, a, b, pd, py, pf, eta, dt) :

    # get deviatoric stress, pressure and plastic multiplyer
    s   = x[0]
    p   = x[1]
    lam = x[2]

    # cutoff negative deviatoric stresses
    if s < 0.0:
        s = 0.0

    # delimiter point
    sd = c + k*pd

    # yield function
    if s*(py - pd) >= (py - p)*sd:

        dfds =  1.
        dfdp = -k

    else:

        Ry   = np.sqrt(s**2 + (p - py)**2)
        dfds = a/Ry*s
        dfdp = a/Ry*(p-py)

    # flow potential
    if s*(pf - pd) >= (pf - p)*sd:

        gs    = 1./2.
        gp    = kf
        dgsds = 0.
        dgsdp = 0.
        dgpds = 0.
        dgpdp = 0.

    else:

        Rf    =  np.sqrt(s**2 + (p - pf)**2)
        gs    =  1./2.*b/Rf*s
        gp    = -b/Rf*(p-pf)
        dgsds =  1./2.*b*(Rf**2-s**2)/Rf**3
        dgsdp = -b*s*(p - pf)/2./Rf**3
        dgpds =  b*s*(p - pf)/Rf**3
        dgpdp = -b * (Rf**2 - (p - pf)**2)/Rf**3

    # get Jacobian
    J[0, 0] = -lam*dgsds -AL - n*AN*s**(n-1.)
    J[0, 1] = -lam*dgsdp
    J[0, 2] = -gs

    J[1, 0] = -lam*dgpds
    J[1, 1] = -lam*dgpdp + 1./K/dt
    J[1, 2] = -gp

    J[2, 0] =  dfds
    J[2, 1] =  dfdp
    J[2, 2] = -eta

#------------------------------------------------------------------------------
def getResults(x, J, dII, K, dt, k, kf, c, b, pd, pf):

    # get deviatoric stress, pressure and plastic multiplyer
    s   = x[0]
    p   = x[1]
    lam = x[2]

    # delimiter point
    sd = c + k*pd

    # flow potential
    if s*(pf - pd) >= (pf - p)*sd:

        gs = 1./2.
        gp = kf

    else:

        Rf =  np.sqrt(s**2 + (p - pf)**2)
        gs =  1./2.*b/Rf*s
        gp = -b/Rf*(p-pf)

    # get plastic strain rates
    dIIpl = lam*gs
    thpl  = lam*gp

    # linearize local iteration
    H =  np.array(((1., 0.), (0., -1./K/dt), (0., 0.)))
    A = -np.linalg.solve(J, H)

    # compute effective stiffness constants
    eta_eff =  s/2./dII
    beta_1  =  A[0, 0]/2. - eta_eff
    beta_2  = -A[1, 0]/2.
    beta_3  =  A[0, 1]
    beta_4  = -A[1, 1]

    return s, p, dIIpl, thpl, eta_eff, beta_1, beta_2, beta_3, beta_4

#------------------------------------------------------------------------------
