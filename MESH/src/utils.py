import numpy as np

#------------------------------------------------------------------------------
def gen_facets(points, add_points):

    #================================================
    # connect additional points with factes in a loop
    # number additonal point after existing ones 
    #================================================

    start = len(points)
    end   = start + len(add_points) - 1

    add_facets = [(i, i+1) for i in range(start, end)] + [(end, start)]

    return add_facets

#------------------------------------------------------------------------------
def add_rectangle(x, y, dx, dy, points, facets):
    
    #======================================================
    # add a rectangular region (generate points and factes)
    #======================================================

    cx = dx/2.0
    cy = dy/2.0

    # add corner points
    add_points = [(x-cx, y-cy), (x+cx, y-cy), (x+cx, y+cy), (x-cx, y+cy)]

    # connect points by factes
    add_facets = gen_facets(points, add_points)

    # add rectangle
    points.extend(add_points)
    facets.extend(add_facets)

#------------------------------------------------------------------------------
def add_circle(x, y, r, n, points, facets):

    #===================================================
    # add a circular region (generate points and factes)
    #===================================================
    
    # generate angles with given step
    phi = np.linspace(0.0, 2.0*np.pi, n, endpoint=False)

    # generate points on the circle
    add_points = [(x+r*np.cos(a), y+r*np.sin(a)) for a in phi]

    # connect points by factes
    add_facets = gen_facets(points, add_points)

    # add circle
    points.extend(add_points)
    facets.extend(add_facets)
    
#------------------------------------------------------------------------------
def add_arc(x, y, r, n, b, e):

    #===================================================
    # add an arc path (generate points only)
    #===================================================
    
    # generate angles with given step
    phi = np.linspace(b, e, n, endpoint=True)

    # generate points on the circle
    arc = [(x+r*np.cos(a), y+r*np.sin(a)) for a in phi]

    return arc

#------------------------------------------------------------------------------

