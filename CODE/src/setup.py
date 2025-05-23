import numpy as np

#------------------------------------------------------------------------------
def initT(temp, nodes, nc, T):

    km    =  1000.0
    grad  =  temp['grad']/km
    top   =  nodes[0:nc, 1].max()
    d     =  top - nodes[0:nc, 1] 
    T [:] =  temp['top'] + d*grad

#------------------------------------------------------------------------------
def initAPS(xsig, xmu, ysig, ymu, daps, size, elems, nodes, Np, nipe, hist):

    # initialize random generator for reproducibility
    np.random.seed(7)

    # initialize number of elements and counter
    ne  = elems.shape[0]
    cnt = 0
    nip = ne*nipe
    npp = int(nip*size)

    pert = np.zeros(nip, dtype=bool)
    pert[np.random.randint(0, nip, npp)] = True

    # ELEMENT LOOP
    for i in range(ne) :

        # get corner indices and coordinates
        cidx  = elems[i, 0:3]
        coord = nodes[cidx, :]

        # INTEGRATION POINT LOOP
        for j in range(nipe) :

            # get coordinate shape functions
            Npj = Np[j]

            # get integration point coordinate
            crd = np.dot(Npj, coord)

            # generate APS perturbation
            x  = crd[0]
            y  = crd[1]
            xg = np.exp(-((x-xmu)**2.0)/(xsig**2.)/2.)
            yg = np.exp(-((y-ymu)**2.0)/(ysig**2.)/2.)

            aps = daps*xg*yg

            if pert[cnt]:
                aps *= 2.0

            hist['svaps'][cnt] = aps

            cnt += 1

#------------------------------------------------------------------------------
def initPf(xsig, xmu, ysig, ymu, psource, nodes, nc, pf):

    # NODE LOOP
    for i in range(nc) :

        # nodal coordinates
        x = nodes[i, 0]
        y = nodes[i, 1]

        xg = np.exp(-((x-xmu)**2.0)/(xsig**2.)/2.)
        yg = np.exp(-((y-ymu)**2.0)/(ysig**2.)/2.)

        pf[i] = psource*xg*yg

#------------------------------------------------------------------------------
def updatePf(xsig, xmu, avsref, elems, nodes, nipe, hist, nc, pf):

    # FLUID PRESSURE QUASI-DIFFUSION IN THE HIGH VOLUMETRIC PLASTIC STRAIN ZONE

    # initialize number of elements and counter
    ne  = elems.shape[0]
    cnt = 0

    # ELEMENT LOOP
    for i in range(ne) :

        # get average avs
        avgavs = 0.0

        # INTEGRATION POINT LOOP
        for j in range(nipe) :

            avgavs += hist['svavs'][cnt]

            cnt += 1

        avgavs /= float(nipe)

        # compute diffusion function
        fdif = avgavs/avsref

        if fdif > 1.0:
            fdif = 1.0

        # get corner indices
        cidx = elems[i, 0:3]

        # get nodal fluid pressure
        p1 = pf[cidx[0]]
        p2 = pf[cidx[1]]
        p3 = pf[cidx[2]]

        # get maximum pressure
        pmax = p1

        if p2 > pmax:
            pmax = p2

        if p3 > pmax:
            pmax = p3

        # get minimum pressure
        pmin = p1

        if p2 < pmin:
            pmin = p2

        if p3 < pmin:
            pmin = p3

        # diffuse pressure
        p1 += fdif*(pmax-pmin)
        p2 += fdif*(pmax-pmin)
        p2 += fdif*(pmax-pmin)

        # limit by maximum
        if p1 > pmax:
            p1 = pmax
        if p2 > pmax:
            p2 = pmax
        if p3 > pmax:
            p3 = pmax

        if  pf[cidx[0]] < p1:
            pf[cidx[0]] = p1

        if  pf[cidx[1]] < p2:
            pf[cidx[1]] = p2

        if  pf[cidx[2]] < p3:
            pf[cidx[2]] = p3

    # NODE LOOP
    for i in range(nc) :

        # nodal coordinates
        x = nodes[i, 0]

        xg = np.exp(-((x-xmu)**2.0)/(xsig**2.)/2.)

        pf[i] *= xg

#------------------------------------------------------------------------------
