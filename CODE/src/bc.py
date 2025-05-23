import numpy as np

#------------------------------------------------------------------------------
def setBCVals(faces, marks, bc_nums, bc_type, bc_vals) :

    idisp = np.zeros([0],   dtype = np.int32)
    itrac = np.zeros([0,0], dtype = np.int32)
    vdisp = np.zeros([0])
    vtrac = np.zeros([0,0])

    for i, ID in enumerate(bc_nums):

        # access bc
        tp = bc_type[i]

        # get subset of boundary faces
        ifaces = faces[marks == ID, :]

        # get unique node indices
        inodes = np.unique(ifaces.flatten())

        if tp == "disp" :

            # get number of nodes
            nn = inodes.shape[0]

            # get displacement components
            dx = bc_vals[i][0]
            dy = bc_vals[i][1]

            # x-displacement component
            if not np.isnan(dx):

                # dof indices
                idof = 2*inodes

                # generate boundary values
                vals = dx*np.ones(nn)

                # store boundary values
                idisp, vdisp = appendData(idisp, vdisp, idof, vals)

            # y-displacement component
            if not np.isnan(dy) :

                # dof indices
                idof = 2*inodes + 1

                # generate boundary values
                vals = dy*np.ones(nn)

                # store boundary values
                idisp, vdisp = appendData(idisp, vdisp, idof, vals)

        elif tp == "trac" :

            # get number of cells
            nc = ifaces.shape[0]

            # get traction components
            tx = bc_vals[i][0]
            ty = bc_vals[i][1]

            # generate boundary values
            vals       = np.zeros((nc, 2))
            vals[:, 0] = tx
            vals[:, 1] = ty

            # store boundary values
            itrac, vtrac = appendData(itrac, vtrac, ifaces, vals)

    # get unique set of boundary constraints
    idisp, idx = np.unique(idisp, return_index=True)
    vdisp      = vdisp[idx]

    return idisp, itrac, vdisp, vtrac

#------------------------------------------------------------------------------
def appendData(idx, val, uidx, uval) :

    if idx.size :

        idx = np.append(idx, uidx, axis=0)
        val = np.append(val, uval, axis=0)

    else:

        idx = uidx
        val = uval

    return idx, val

#------------------------------------------------------------------------------
