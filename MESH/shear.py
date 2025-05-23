import numpy as np
from src.utils import add_arc, gen_facets
from src.meshpy_triangle_api import runMeshPy

#------------------------------------------------------------------------------
def refinement(vertices, area):

    return False

#------------------------------------------------------------------------------
def main():

    filename = 'shear'
    L        = 1.0
    H        = 0.7
    xc       = 0.0
    yc       = 0.0
    r        = 0.025
    n        = 5
    area     = 3e-4
    angle    = 30.0
    points   = []
    facets   = []
    holes    = []
    regions  = []


    #============================
    # initialize geometry objects
    #============================
    
    # generate outer contour with an arc
    arc        = add_arc(xc, yc, r, n, 0.0, np.pi/2) 
    add_points = arc + [(0.0, H), (L, H), (L, 0.0)]
    facets     = gen_facets(points, add_points)
    points     = points + add_points

    # close arc from the bottom left corner 
    m = len(points)
    points = points + [(xc, yc)]
    facets = facets + [(0, m), (m, n-1)]

    # create basic region (attribute 1)
    regions.append((L/2, H/2, 1, area))

    # create weak inclusion (attribute 2)
    regions.append((xc+r/2, yc+r/2, 2, area))

    #============================
    # initialize boundary markers
    #============================

    nf            = n-1
    markers       = [0]*len(facets)
    markers[nf]   = 2 # left facet
    markers[nf+5] = 2 #
    markers[nf+2] = 3 # right facet
    markers[nf+3] = 4 # bottom facet
    markers[nf+4] = 4
    markers[nf+1] = 5 # top facet

    runMeshPy(points, facets, markers, holes, regions, area, angle, filename, refinement)

#------------------------------------------------------------------------------
if __name__ == "__main__":

    main()

#------------------------------------------------------------------------------