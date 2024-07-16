"""
Making cylinders.

By: Cléo Delêtre
Date: July 16th, 2024
"""
# Import modules:
import gmsh
import numpy as np
import sys

def cylinder(x0, radius, n_points, height):
    """
    Creates a gmsh of the rounded blade and saves it in ./gmsh
    :param radius: Radius at the bottom of the blade.
    :param n_points: Number of points used to create the radius on each sides of the blade
    :param angle_fraction: The last point on the radius will be at this value. Usually between 0.5 and 1.0
    :param dept: Thickness of the blade in the z direction.
    """

    # Initialize parameters
    x1 = x0 + height
    # Initialize gmsh:
    gmsh.initialize()

    # Delta theta used for the radius
    delta_theta = 2 * np.pi / (n_points - 1)

    # Initialize arrays
    points_id_0 = []
    points_id_1 = []
    lines_id_x0_2_x1 = []
    lines_id_diag = []
    lines_id_x0 = []
    lines_id_x1 = []
    line_loops = []
    plane_surfaces = []

    # Points at z0 and z1 
    for i in range(n_points):
        points_id_0.append(
            gmsh.model.geo.addPoint(x0,
                                    radius * (np.cos(i * delta_theta)), radius * (np.sin(i * delta_theta)),
                                    meshSize=1.))
        points_id_1.append(
            gmsh.model.geo.addPoint(x1,
                                    radius * (np.cos(i * delta_theta)), radius * (np.sin(i * delta_theta)),
                                    meshSize=1.))


    # Now we need to make triangles.
    # Z1[n]   <----------------- Z1[n]
    #  ^    \__                    ^
    #  |       \__                 |
    #  |          \__              |
    #  |             \__           |
    #  |                \__        |
    #  |                   \__>    |
    # Z1[n+1] <----------------- Z1[n+1]
    #
    for index, p0 in enumerate(points_id_0):
        # Lines going from x0 to x1
        lines_id_x0_2_x1.append(gmsh.model.geo.add_line(p0, points_id_1[index]))

    for j in range(len(points_id_0) - 1):
        # Lines going from x1 to x0[ + 1]
        lines_id_diag.append(gmsh.model.geo.add_line(points_id_1[j], points_id_0[j + 1]))

    for j in range(len(points_id_0) - 1):
        # Lines going from x0[+1] to x0[]
        lines_id_x0.append(gmsh.model.geo.add_line(points_id_0[j + 1], points_id_0[j]))
        lines_id_x1.append(gmsh.model.geo.add_line(points_id_1[j + 1], points_id_1[j] ))

    for j in range(len(points_id_1) - 1):
        line_loops.append(gmsh.model.geo.addCurveLoop([lines_id_x0_2_x1[j], lines_id_diag[j], lines_id_x0[j]]))
        line_loops.append(gmsh.model.geo.addCurveLoop([lines_id_diag[j], lines_id_x0_2_x1[j + 1], lines_id_x1[j]]))

    for j in range(len(line_loops)):
        plane_surfaces.append(gmsh.model.geo.addPlaneSurface([line_loops[j]]))

    gmsh.model.geo.synchronize()
    gmsh.model.addPhysicalGroup(2, plane_surfaces)
    gmsh.model.geo.synchronize()
    
    # Generate mesh:
    gmsh.model.mesh.generate(2)

    # Write mesh data:
    gmsh.write("cylinder.msh")

    # Creates  graphical user interface
    """
    if 'close' not in sys.argv:
        gmsh.fltk.run()
    """
    # Finalize the Gmsh API
    gmsh.finalize()

    print(f"cylinder.msh saved")
    return

def support(x0, radius, n_points, height):
    """
    Creates a gmsh of the rounded blade and saves it in ./gmsh
    :param radius: Radius at the bottom of the blade.
    :param n_points: Number of points used to create the radius on each sides of the blade
    :param angle_fraction: The last point on the radius will be at this value. Usually between 0.5 and 1.0
    :param dept: Thickness of the blade in the z direction.
    """

    # Initialize parameters
    x1 = x0 + height
    # Initialize gmsh:
    gmsh.initialize()

    # Delta theta used for the radius
    delta_theta = 2 * np.pi / (n_points - 1)

    # Initialize arrays
    points_id_0 = []
    points_id_1 = []
    lines_to_origin_x0 =[]
    lines_to_origin_x0_inverse =[]
    lines_to_origin_x1 =[]
    lines_to_origin_x1_inverse =[]
    lines_id_x0_2_x1 = []
    lines_id_diag = []
    lines_id_x0 = []
    lines_id_x1 = []
    line_loops = []
    plane_surfaces = []
    

    # Points at origin
    P0_x0 = gmsh.model.geo.addPoint(x0,0,0)
    P0_x1 = gmsh.model.geo.addPoint(x1,0,0)

    # Points at x0 and x1 
    for i in range(n_points):
        points_id_0.append(
            gmsh.model.geo.addPoint(x0,
                                    radius * (np.cos(i * delta_theta)), radius * (np.sin(i * delta_theta)),
                                    meshSize=1.))
        points_id_1.append(
            gmsh.model.geo.addPoint(x1,
                                    radius * (np.cos(i * delta_theta)), radius * (np.sin(i * delta_theta)),
                                    meshSize=1.))


    # Now we need to make triangles.
    # Z1[n]   <----------------- Z1[n]
    #  ^    \__                    ^
    #  |       \__                 |
    #  |          \__              |
    #  |             \__           |
    #  |                \__        |
    #  |                   \__>    |
    # Z1[n+1] <----------------- Z1[n+1]
    #
    for index, p0 in enumerate(points_id_0):
        # Lines going from x0 to x1
        lines_id_x0_2_x1.append(gmsh.model.geo.add_line(p0, points_id_1[index]))
        #Lines going from origin to circle
        lines_to_origin_x0.append(gmsh.model.geo.add_line(P0_x0,points_id_0[index]))
        lines_to_origin_x0_inverse.append(gmsh.model.geo.add_line(points_id_0[index], P0_x0))
        lines_to_origin_x1.append(gmsh.model.geo.add_line(P0_x1,points_id_1[index]))
        lines_to_origin_x1_inverse.append(gmsh.model.geo.add_line(points_id_1[index], P0_x1))

    for j in range(len(points_id_0) - 1):
        # Lines going from x1 to x0[ + 1]
        lines_id_diag.append(gmsh.model.geo.add_line(points_id_1[j], points_id_0[j + 1]))

    for j in range(len(points_id_0) - 1):
        # Lines going from x0[+1] to x0[]
        lines_id_x0.append(gmsh.model.geo.add_line(points_id_0[j + 1], points_id_0[j]))
        lines_id_x1.append(gmsh.model.geo.add_line(points_id_1[j + 1], points_id_1[j] ))
    
    # Line loops

    for j in range(len(points_id_1) - 1):
        line_loops.append(gmsh.model.geo.addCurveLoop([lines_id_x0_2_x1[j], lines_id_diag[j], lines_id_x0[j]]))
        line_loops.append(gmsh.model.geo.addCurveLoop([lines_id_diag[j], lines_id_x0_2_x1[j + 1], lines_id_x1[j]]))

        #line loops for base
        line_loops.append(gmsh.model.geo.addCurveLoop([lines_to_origin_x0[j+1], lines_id_x0[j], lines_to_origin_x0_inverse[j]]))
        line_loops.append(gmsh.model.geo.addCurveLoop([lines_to_origin_x1[j+1], lines_id_x1[j], lines_to_origin_x1_inverse[j]]))

    for j in range(len(line_loops)):
        plane_surfaces.append(gmsh.model.geo.addPlaneSurface([line_loops[j]]))

    gmsh.model.geo.synchronize()
    gmsh.model.addPhysicalGroup(2, plane_surfaces)
    gmsh.model.geo.synchronize()
    
    # Generate mesh:
    gmsh.model.mesh.generate(2)

    # Write mesh data:
    gmsh.write("support.msh")

    # Creates  graphical user interface
    """
    if 'close' not in sys.argv:
        gmsh.fltk.run()
    """
    # Finalize the Gmsh API
    gmsh.finalize()

    print(f"support.msh saved")
    return

# (x0, radius, n_points, height)
cylinder(-0.005, 0.005, 50, 0.02)
support(-0.01 ,0.005 ,50, 0.005)