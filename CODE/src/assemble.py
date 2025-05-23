import numpy as np
from src.update import stressUpdate
from scipy.sparse import coo_matrix

#------------------------------------------------------------------------------
def getElemMatVecMech(nodes, elems, params, phases, attr, hist, svar, dbcflags, U, pf, T, Np, dNp, Nu, dNu, w, rows, cols, vals, Fint, Fext, actnl, it, dt) :

    nedof = 17
    ne    = elems.shape[0]
    nn    = nodes.shape[0]
    nipe  = w.shape[0]
    Ke    = np.empty((nedof, nedof))
    fint  = np.empty((nedof))
    fext  = np.empty((nedof))
    idof  = np.empty((nedof), dtype=np.int32)
    Nm    = np.zeros((2, 14))
    Bm    = np.zeros((4, 14))

    # get gravity acceleration vector
    g = params['g']
    g = np.array((0., -g))

    # matrix constants
    m  = np.array((1., 1., 1., 0.))
    ID = 1./2.*np.diag(np.array((2., 2., 2., 1.))) - 1./3.*np.outer(m, m)

    # zero out residual
    Fint[:] = 0.0
    Fext[:] = 0.0

    cnt  = 0
    ndiv = 0
    vres = 0.0

    # deactivate nonlinearity for the first iteration
    if it == 0:
        actnl = 0

    # ELEMENT LOOP
    for i in range(ne) :

        # get corner indices and coordinates
        cidx  = elems[i, 0:3]
        coord = nodes[cidx, :]

        # get dof indices
        uidx = elems[i, :]
        idof[np.arange(0,  14, 2)] = 2*uidx
        idof[np.arange(1,  14, 2)] = 2*uidx + 1
        idof[np.arange(14, 17, 1)] = np.arange(2*nn + 3*i, 2*nn + 3*(i+1), 1)

        # get element displacement and pressure increments, fluid pressure and temperature
        Ue   = U [idof[0:14]]
        pe   = U [idof[14:17]]
        pfe  = pf[cidx]
        Te   = T [cidx]

        # setup element stiffness matrix and RHS vector
        Ke  [:] = 0.0
        fint[:] = 0.0
        fext[:] = 0.0

        # INTEGRATION POINT LOOP
        for j in range(nipe) :

            # get shape functions, derivatives & integration weight
            Npj  =  Np[j]
            dNpj = dNp[j]
            Nuj  =  Nu[j]
            dNuj = dNu[j]
            wj   =   w[j]

            # get Jacobian & determinant
            J    = np.dot(dNpj, coord)
            detJ = np.linalg.det(J)

            # get shape functions derivatives in global coordinates
            dNuj = np.linalg.solve(J, dNuj)

            # shape function and derivative matrix
            for k in range(7) :

                # shape function matrix
                Nm[0, 2*k  ] = Nuj[k]
                Nm[1, 2*k+1] = Nuj[k]

                # strain matrix
                Bm[0, 2*k  ] = dNuj[0, k]
                Bm[1, 2*k+1] = dNuj[1, k]
                Bm[3, 2*k  ] = dNuj[1, k]
                Bm[3, 2*k+1] = dNuj[0, k]

            Nt = Nm.T
            Bt = Bm.T

            # compute strain and pressure increments
            deps = np.dot(Bm,  Ue)
            dp   = np.dot(Npj, pe)

            # interpolate fluid pressure and temperature
            pfp = np.dot(Npj, pfe)
            Tp  = np.dot(Npj, Te)

            # compute stress, tangent stiffness matrix, and solution variales
            sig, rth, nrth, Duu, Dup, Dpu, Dpp, rhob, div = stressUpdate(deps, dp, pfp, Tp, hist, svar, cnt, params, phases[attr[i]], ID, m, actnl, dt)

            # get element tangent stiffness matrix and residual vector
            Ke[ 0:14,  0:14] += np.dot  (Bt,              np.dot(Duu, Bm)) * detJ * wj
            Ke[ 0:14, 14:17] += np.outer(np.dot(Bt, Dup), Npj)             * detJ * wj
            Ke[14:17,  0:14] += np.outer(Npj,             np.dot(Dpu, Bm)) * detJ * wj
            Ke[14:17, 14:17] += np.outer(Npj,             Npj)*Dpp         * detJ * wj

            fint[ 0:14] += np.dot(Bt, sig)    * detJ * wj
            fint[14:17] += Npj*rth            * detJ * wj
            fext[ 0:14] += np.dot(Nt, rhob*g) * detJ * wj

            # update iterator & diverged point counter
            cnt  += 1
            ndiv += div

            # update volumetric strain rate residual
            if nrth > vres:
                vres = nrth

        # constrain matrix
        for k in range(nedof) :

            if dbcflags[idof[k]] :
                Ke[k, :] = 0.0
                Ke[:, k] = 0.0
                Ke[k, k] = 1.0

        # store local matrix & assemble local vector
        sidx        = np.arange(i*nedof**2, (i+1)*nedof**2)
        irow        = np.repeat(idof, nedof).reshape(nedof, nedof)
        icol        = irow.T
        rows[sidx]  = irow.flatten()
        cols[sidx]  = icol.flatten()
        vals[sidx]  = Ke.flatten()
        Fint[idof] += fint
        Fext[idof] += fext
        
    return vres, ndiv

#------------------------------------------------------------------------------
def getFaceTracVecMech(nodes, itrac, vtrac, params, N, dNc, w, Fext) :

    nedof = 6
    ne    = itrac.shape[0]
    nipe  = w.shape[0]
    fext  = np.empty((nedof))
    idof  = np.empty((nedof), dtype=np.int32)
    Nm    = np.zeros((2, nedof))

    # FACE LOOP
    for i in range(ne) :

        # get nodal indices & coordinates, & surface traction
        uidx  = itrac[i,    :    ]
        cidx  = itrac[i,    0:3:2]
        coord = nodes[cidx, :    ]
        trac  = vtrac[i,    :    ]

        # get dof indices
        idof[np.arange(0, nedof, 2)] = 2*uidx
        idof[np.arange(1, nedof, 2)] = 2*uidx+1

        # setup local vector
        fext[:] = 0.0

        # INTEGRATION POINT LOOP
        for j in range(nipe) :

            # get shape functions, derivatives & integration weight
            Nj   =  N [j]
            dNcj = dNc[j]
            wj   =   w[j]

            # shape function matrix
            for k in range(3) :

                Nm[0, 2*k]   = Nj[k]
                Nm[1, 2*k+1] = Nj[k]

            Nt = Nm.T

            # compute tangent vector
            t = np.dot(dNcj, coord)

            # get outward normal vector
            n = np.array([t[1], -t[0]])

            # get area element
            dA = np.linalg.norm(n)

            # integrate
            fext += np.dot(Nt, trac) * dA * wj

        # assemble
        Fext[idof] += fext

#------------------------------------------------------------------------------
def assembleMech(ndof, udof, rows, cols, vals, Fint, Fext, F, dbcflags, params) :

    # assemble global sparse matrix in compressed sparse column format
    K = coo_matrix((vals, (rows, cols)), shape=(ndof, ndof)).tocsc()

    # assemble residual
    F[:] = Fint - Fext

    # constrain residual
    F[dbcflags] = 0.0

    # add reaction forces to the external force vector
    Fext[dbcflags] -= Fint[dbcflags]

    # compute normalized force residual norm
    fres = np.abs(F[0:udof]).max()/np.abs(Fext[0:udof]).max()

    return K, fres

#------------------------------------------------------------------------------
