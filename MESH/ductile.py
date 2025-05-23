import numpy as np
from src.utils import add_rectangle
from src.meshpy_triangle_api import runMeshPy

#------------------------------------------------------------------------------
def refinement(vertices, area):

    # define refined domain and mesh sizes
    xl      = 30000.0
    xr      = 70000.0
    yb      = 18000.0
    refined = 60.0
    average = 500.0
    maxdist = 20000.0
    
    # get node coordinates
    x = np.array(vertices)

    # get element center coordinates
    cent = np.sum(x, axis=0)/3.0
    xe   = cent[0]
    ye   = cent[1]
    
    # compute distance to refined domain
    if xe < xl:
        if ye < yb:
            d = np.sqrt((xe-xl)**2.0 + (ye-yb)**2.0)
        else:
            d = xl-xe
    elif xe > xr:
        if ye < yb:
            d = np.sqrt((xe-xr)**2.0 + (ye-yb)**2.0)
        else:
            d = xe-xr
    elif ye < yb:
        d = yb-ye
    else:
        d = 0.0
    
    # get target element size by linear interpolation
    if d == 0.0:
        target = refined
    elif d > maxdist:
        target = average
    else:
        target = refined + (average - refined)*(d/maxdist)

    # set refinement flag
    if area > target**2:
 	    refine = True
    else:
        refine = False

    return refine

#------------------------------------------------------------------------------
def main():

    filename = 'ductile'
    cx       = 50000.0
    cy       = 12500.0
    L        = 100000.0
    H        = 25000.0
    area     = 1000000.0
    angle    = 30.0
    points   = []
    facets   = []
    holes    = []
    regions  = []

    #============================
    # initialize geometry objects
    #============================
 
    # create basic region (attribute 1)
    add_rectangle(cx, cy, L, H, points, facets)
    regions.append((cx, cy, 1, area))

    #============================
    # initialize boundary markers
    #============================
    
    markers    = [0]*len(facets)
    markers[3] = 2 # left facet
    markers[1] = 3 # right facet
    markers[0] = 4 # bottom facet

    runMeshPy(points, facets, markers, holes, regions, area, angle, filename, refinement)

#------------------------------------------------------------------------------
if __name__ == "__main__":

    main()

#------------------------------------------------------------------------------