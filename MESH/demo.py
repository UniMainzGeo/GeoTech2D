import numpy as np
from src.utils import add_rectangle, add_circle
from src.meshpy_triangle_api import runMeshPy

#------------------------------------------------------------------------------
def refinement(vertices, area):

    # define refinement domain (circle) and target element area
    xc     =  -1.0
    yc     =   0.5
    r      =   0.5
    target =   0.005
    
    # get node coordinates
    x = np.array(vertices)

    # get element center coordinates
    cent = np.sum(x, axis=0)/3.0
    xe   = cent[0]
    ye   = cent[1]

    # get distance to the center of the refinement domain
    d = np.sqrt((xe-xc)**2.0 + (ye-yc)**2.0)

    # refine large elements withing the refinement domain
    if(d < r) and (area > target):
 	    refine = True
    else:
        refine = False

    return refine

#------------------------------------------------------------------------------
# 
# use this function if no refinement is required
# 
# def refinement(vertices, area):
# 
#     return False
# 
#------------------------------------------------------------------------------
def main():
    
    filename = 'demo'
    area     = 0.3
    refined  = 0.1
    angle    = 30.0
    points   = []
    facets   = []
    holes    = []
    regions  = []

    #============================
    # initialize geometry objects
    #============================
    
    # every region or hole is specified by internal point
    # regions accept additionally:
    #   * attribute (unit or material) index
    #   * element volume (area)

    # create basic region (attribute 1)
    add_rectangle(0.0, 0.0, 6.0, 6.0, points, facets)
    regions.append((0.0, 0.0, 1, area))

    # add circular hole
    add_circle( 1.5,  1.5, 0.5, 30, points, facets)
    holes.append((1.5, 1.5))

    # add circular region (attribute 2)
    add_circle(-1.5, -1.5, 0.8, 30, points, facets)
    regions.append((-1.5, -1.5, 2, area))

    # add rectangular hole
    add_rectangle(-1.7, 2.0, 1.5, 0.7, points, facets)
    holes.append((-1.7, 2.0))

    # add rectangular region (attribute 3, refined elements)
    add_rectangle(1.5, -1.5, 1.0, 2.0, points, facets)
    regions.append((1.5, -1.5, 3, refined*area))

    #============================
    # initialize boundary markers
    #============================

    # set boundary indices starting from 2
    # 0 is reserved for unassigned inner faces
    # 1 is reserved for unassigned outer faces

    # use factes of the first object (basic region)
    # leave top facet unassigned
    
    markers    = [0]*len(facets)
    markers[0] = 2 # bottom facet
    markers[1] = 3 # right facet
    markers[3] = 4 # left facet

    runMeshPy(points, facets, markers, holes, regions, area, angle, filename, refinement)

#------------------------------------------------------------------------------
if __name__ == "__main__":

    main()

#------------------------------------------------------------------------------
