===============================
Introduction on how to use GMSH
===============================

--------------------------
Installation
--------------------------

--------------------------
Geometry
--------------------------

""""""""""""""""""""""""""
Built-in kernel
""""""""""""""""""""""""""

""""""""""""""""""""""""""
OpenCASCADE kernel
""""""""""""""""""""""""""
OpenSource external CAD engined library.

""""""""""""""""""""""""""
Import CAD file
""""""""""""""""""""""""""
Importing CAD files (``.step`` or ``.stp`` format) can be particularly convenient for more complex fluid geometries (like pipes) or immersed solids (like an impeller):

.. hint::
  In the case of immersed solids, use a simplified CAD file of the outer shell of the solid, e.g. without any screws or bolts or threads.

1. ``File > New...``: create a new .geo file (can use OpenCASCADE or Built-in kernel)
2. ``Files > Merge...``: merge the CAD file (``.step`` or ``.stp`` format) with GMSH
3. ``Tools > Statistics``: check that the geometry is loaded (point, curves, surfaces, and if 3D volumes)

.. seealso::
  You can find a step-by-step video `here <https://www.youtube.com/watch?v=e7zA3joOWX8>`_, with very useful tools as how to inspect your mesh.

--------------------------
Physical group
--------------------------

.. warning::
  This step is essential for the mesh compatibility with Lethe

---------------------------
Mesh
---------------------------

""""""""""""""""""""""""""
Unstructured
""""""""""""""""""""""""""

Basic:

1. ``Left pannel: Modules > Mesh > 2D`` or ``3D``: create the mesh
2. ``Tools > Statistics``: check that the mesh is generated appropriately (by default, triangles for 2D and hexahedra for 3D)
3. (optional) ``Left pannel: Modules > Mesh > Refine by splitting``: refine the mesh (beware, it takes more and more time for each refinement)
4. ``Left pannel: Modules > Mesh > Save``: save the mesh in a ``.msh`` file, to be used in Lethe (see :doc:`../../parameters/cfd/mesh`)

""""""""""""""""""""""""""
Structured
""""""""""""""""""""""""""

Basic, from an unstructured mesh:

1. ``Left pannel: Modules > Mesh > 2D`` or ``3D``: create the mesh
2. ``Tools > Statistics``: check that the mesh is generated appropriately (by default, triangles for 2D and hexahedra for 3D)
3. (optional) ``Left pannel: Modules > Mesh > Refine by splitting``: refine the mesh (beware, it takes more and more time for each refinement)
4. (optional) ``Tools > Options > Mesh`` and ``General`` panel, check ``Recombine all triangular meshes``: generate a quad mesh. If necessary, the 


--------------------------
Hints
--------------------------


