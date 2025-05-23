from src.utils import add_rectangle
from src.meshpy_triangle_api import runMeshPy

#------------------------------------------------------------------------------
def refinement(vertices, area):

    return False

#------------------------------------------------------------------------------
def main():

    filename = 'crust'
    cx       = 20000.0
    cy       = 3500.0
    L        = 40000.0
    H        = 7000.0
    area     = 2500.0
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
