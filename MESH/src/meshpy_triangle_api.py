import os
import shutil
import copy
import numpy as np
import distinctipy
import matplotlib.pyplot as plt
import meshpy.triangle as triangle

#------------------------------------------------------------------------------
def runMeshPy(points, facets, markers, holes, regions, area, angle, filename, refinement):

#==============================================================================
# INPUT PARAMETERS
# 
# points   - key point coordinates
# facets   - connectivity between the points (each facet connects two points)
# markers  - boundary identifier for each facet (set to 0 by default) 
#            unmarked faces will be assigned to 0 (internal), and to 1 (external)
#            start marker numbering with 2 to ensure unique identification 
# holes    - internal point coordinate for each hole
# regions  - internal point coordinate, material unit ID, and target element area for each region
# area     - default element area
# angle    - minimum element angle
# filename - name of the binary file to store the generated mesh
#
# NOTE: points and facets should ideally form a "clean" Planar Straight Line Graph (PSLG)
#       that means nodes must be inserted in all the facet intersections
#
# refinement is a function used to refine the mesh, it has the following API:
#
#   def refinement(vertices, area):
#
#      vertices - coordinates of the corners of a triangle
#      area     - area of a triangle
#
# return value:
#   
#   True  - if a triangle requires refinement
#   False - otherwise
#
# If no refinement is necessary just directly return False
#==============================================================================

    # create output directory
    if not os.path.isdir('mesh'):
        os.makedirs('mesh')    

    plt.close('all')

    # create input object
    info = triangle.MeshInfo()
    info.set_points(points)
    info.set_facets(facets, markers)
    if holes:
        info.set_holes(holes)
    info.regions.resize(len(regions))
    for i, region in enumerate(regions) :
        info.regions[i] = copy.copy(region)

    # call Triangle library
    mesh = triangle.build(info,
                          max_volume             = area,
                          min_angle              = angle,
                          volume_constraints     = True,
                          attributes             = True,
                          allow_boundary_steiner = True,
                          refinement_func        = refinement)

    # get nodes, elements, boundary faces, and face markers
    nodes = np.array(mesh.points)
    elems = np.array(mesh.elements)
    faces = np.array(mesh.facets)
    marks = np.array(mesh.facet_markers)
    units = np.array(mesh.element_attributes, dtype=np.int32)
    
    # get number of corner nodes
    nc = nodes.shape[0]

    # add nodes on the edges and in the centers
    nodes, elems, faces = addNodesEdgeCenter(nodes, elems, faces)

    # plot elements (use units as color)
    plt.figure()
    plt.tripcolor(nodes[:, 0], nodes[:, 1], elems[:, 0:3], facecolors=units, edgecolors='w')
    plt.title('Elements and units')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('equal')
    plt.draw()
    plt.pause(1e-1)
    
    # store basic output to file
    np.savez(os.path.join("mesh", filename + '.npz'),
        nc    = nc,
        nodes = nodes,
        elems = elems,
        faces = faces,
        marks = marks,
        units = units)

    # plot faces (use markers as colors)
    # black color indicates unassigned faces
    #    * dotted line - inner faces (marker 0)
    #    * dashed line - outer faces (marker 1)
    # solid color lines indicate assigned faces (marker 2 and higher)
    
    umarks = np.unique(marks)
    nmarks = umarks.shape[0]
    colors = []
    styles = []
    widths = []
    cnt    = 0

    if(umarks[cnt] == 0):
        colors.append((0,0,0))
        styles.append(':')
        widths.append(1)
        cnt     += 1

    if(umarks[cnt] == 1):
        colors.append((0,0,0))
        styles.append('--')
        widths.append(1)
        cnt     += 1

    ncolor = nmarks - cnt
    colors = colors + distinctipy.get_colors(ncolor)
    styles = styles + ['-']*(ncolor)
    widths = widths + [ 3 ]*(ncolor)

    plt.figure()
    for i in range(nmarks):

        mfaces = faces[marks == umarks[i], 0:3:2]

        # get start & end coordinates of the edges
        xs = nodes[mfaces[:, 0], 0]
        ys = nodes[mfaces[:, 0], 1]
        xe = nodes[mfaces[:, 1], 0]
        ye = nodes[mfaces[:, 1], 1]

        # plot edges
        x = np.vstack((xs, xe))
        y = np.vstack((ys, ye))

        plt.plot(x, y, linestyle=styles[i], color=colors[i], linewidth=widths[i])

    plt.title('Boundary faces and markers')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('equal')
    plt.draw()
    plt.pause(1e-1)

    plt.show()
    
    print('=========================================')
    print('Output file name  : ', filename + '.npz')
    print('Number of elements: ', elems.shape[0])
    print('Number of nodes   : ', nodes.shape[0])
    print('Number of corners : ', nc)
    print('=========================================')

    shutil.rmtree(os.path.join('src', '__pycache__'))

#------------------------------------------------------------------------------
def addNodesEdgeCenter(nodes, elems, faces) :

    #==================================================
    # add nodes to the edges, faces and element centers
    #==================================================
    
    # generate all edges in elements and faces
    n1    = np.concatenate((elems[:, 0], elems[:, 1], elems[:, 2], faces[:, 0]))
    n2    = np.concatenate((elems[:, 1], elems[:, 2], elems[:, 0], faces[:, 1]))
    edges = np.column_stack((n1, n2))

    # sort eges for unique identification
    edges = np.sort(edges, axis=1)

    # remove edge duplicates, get unique numbering
    # edgnum array contains unique global indices of all redundant eges
    # in the order of generation: element edges, face edges
    edges, edgnum = np.unique(edges, return_inverse=True, axis=0)

    # add edge node indices
    nn    = nodes.shape[0]
    ne    = elems.shape[0]
    nf    = faces.shape[0]
    edge1 = edgnum[   0 :   ne]
    edge2 = edgnum[  ne : 2*ne]
    edge3 = edgnum[2*ne : 3*ne]
    edgef = edgnum[3*ne : 3*ne+nf]
    elems = np.column_stack((elems, edge1+nn, edge2+nn, edge3+nn))
    faces = np.column_stack((faces[:, 0], edgef+nn, faces[:, 1]))

    # add edge node coordinates
    edgcrd = (nodes[edges[:, 0], :] + nodes[edges[:, 1], :])/2.0
    nodes  = np.vstack((nodes, edgcrd))

    # add center node indices
    nn      = nodes.shape[0]
    centers = np.arange(0, ne)
    elems   = np.column_stack((elems, centers+nn))

    # add center node coordinates
    cencrd = (nodes[elems[:, 0], :] + nodes[elems[:, 1], :] + nodes[elems[:, 2], :])/3.0
    nodes  = np.vstack((nodes, cencrd))

    return nodes, elems, faces

#------------------------------------------------------------------------------

