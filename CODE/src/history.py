import numpy as np
from dataclasses import dataclass
import src.scale as scl

#------------------------------------------------------------------------------
# history variables nomenculature including output labels and scales
# change it here, if necessary, and the rest will be parametrized 
#------------------------------------------------------------------------------
@dataclass
class var_dtype:
    
    vname  : str   # variable name
    ncomp  : int   # number of components
    label  : list  # output lables
    ishist : bool  # history variable flag
    isout  : bool  # output variable flag
    pscale : float # output scaling factor

#------------------------------------------------------------------------------
def getHistVar():
  
    stress_lables = ['TXX  [MPa  ]', 'TYY  [MPa  ]', 'TZZ  [MPa  ]', 'TXY  [MPa  ]']
    
    hist_var = []
    
    hist_var.append(var_dtype('svtau', 4, stress_lables,    True,  True, -scl.MPa))
    hist_var.append(var_dtype('svp',   1, ['P    [MPa  ]'], True,  True,  scl.MPa))
    hist_var.append(var_dtype('svaps', 1, ['APS  [     ]'], True,  True,  1.0))
    hist_var.append(var_dtype('svavs', 1, ['AVS  [     ]'], True,  True,  1.0))
    hist_var.append(var_dtype('sveII', 1, ['EII  [1/s  ]'], False, True,  1.0))
    hist_var.append(var_dtype('svsII', 1, ['SII  [MPa  ]'], False, True,  scl.MPa))
    hist_var.append(var_dtype('svth',  1, ['TH   [1/s  ]'], False, True,  1.0))

    return hist_var

#------------------------------------------------------------------------------
def allocHistVar(ne, nipe):

    sz    = ne*nipe
    hist  = {}
    svar  = {}
   
    hist_var = getHistVar()
 
    for var in hist_var:

        if var.ncomp > 1:

            svar[var.vname] = np.zeros((sz, var.ncomp))
        
        else:
            
            svar[var.vname] = np.zeros((sz))
     
        if var.ishist:

            if var.ncomp > 1:

                hist[var.vname] = np.zeros((sz, var.ncomp))
        
            else:
            
                hist[var.vname] = np.zeros((sz))      

    return hist, svar

#------------------------------------------------------------------------------
def storeHistVar(hist, svar):
        
    hist_var = getHistVar()
 
    for var in hist_var:

        if var.ishist:
            
            hist[var.vname][:] = svar[var.vname][:]

#------------------------------------------------------------------------------
def outputSolVar(elems, nc, nipe, pointData, svar):

    nelnod  = np.zeros((nc))
    nodcomp = np.zeros((nc))

    hist_var = getHistVar()
    
    for var in hist_var:

        if var.isout:

            if var.ncomp > 1:
                
                for i in range(var.ncomp):
        
                    interpSolVarComp(elems, nipe, svar[var.vname][:, i], nodcomp, nelnod)
                    
                    pointData[var.label[i]] = nodcomp/nelnod/var.pscale

            else:
                
                interpSolVarComp(elems, nipe, svar[var.vname], nodcomp, nelnod)
                
                pointData[var.label[0]] = nodcomp/nelnod/var.pscale


#------------------------------------------------------------------------------
def interpSolVarComp(elems, nipe, svcomp, nodcomp, nelnod) :

    # sum up solution variable from integration points to nodes
    # count number of contributions for each node for normalization

    # setup integration point index for each corner node
    if nipe == 3:

        # N1 -> p2
        # N2 -> p3
        # N3 -> p1

        pidx = np.array((2, 3, 1), dtype=np.int32)

    elif nipe == 6:

        # N1 -> p4
        # N2 -> p5
        # N3 -> p6

        pidx = np.array((4, 5, 6), dtype=np.int32)

    nodcomp[:] = 0.0
    nelnod [:] = 0.0

    ne  = elems.shape[0]
    cnt = 0

    for i in range(ne) :

        idx = elems[i, 0:3]

        nelnod[idx] += 1.0

        for j in range(3) :

            nodcomp[idx[j]] += svcomp[cnt + pidx[j] - 1]

        cnt += nipe

#------------------------------------------------------------------------------
