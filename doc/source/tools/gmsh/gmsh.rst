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
In the gmsh geometry section of the gmsh GUI (see **Geometry** > **Elementary entities** > **Add**), you can add directly multiple 2D or 3D common geometries with a simple click thanks to OpenCASCADE kernel. After selecting the geometry that you need, gmsh will automatically open a window where you can easily set the characteristic lenghts of the geometry. By doing so, the .geo file will be updated.

.. note::
	At anytime, you can go back to your code to change manually some parameters of your geometry. In the **Elementary entities** section, you can **Edit script** to go back to the code and **Reload** script to visualize and update the .geo file on the gmsh GUI.
	
To represent a circle in a rectangle to simulate the flow around a sphere in 2D, we will first select the **Circle** geometry available with OpenCASCADE, set the radius to 1 and center it at :math:`(x,y)=(0,0)`. Then, we will select the **Rectangle** geometry, set the length (DX) to 20, the width (DY) to 10 and the left bottom corner to :math:`(x,y)=(-5,-5)`. The geometry should look like below.

.. image:: images/geo.png
    :alt: The geometry of a circle in a rectangle
    :align: center
    :name: geometry
    
.. note::
	If you click on **Edit script**, you will see that the OpenCASCADE kernel has been added to the code as ``SetFactory("OpenCASCADE");``. Also, the circle is set with ``Circle(1) = {0, 0, 0, 1, 0, 2*Pi};`` and the rectangle with ``Rectangle(1) = {-5, -5, 0, 20, 10, 0};``.

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


