===============================
Introduction on how to use GMSH
===============================

--------------------------
Installation
--------------------------

--------------------------
Geometry
--------------------------

.. image:: images/geo.png
    :alt: The geometry of a circle in a rectangle
    :align: center
    :name: geometry

""""""""""""""""""""""""""
Built-in kernel
""""""""""""""""""""""""""
It is quite easy to create a .geo file directly by coding line by line the geometry that you want to produce. You usally start by setting the points that delimits the shape being traced. The attribut ``Point`` is essential to set a point of the shape and to specify its coordinates. For the rectangle, we will have:

.. code-block::

	// Point(id) = {<X>, <Y>, <Z>, <Prescribed mesh size at point>}
	Point(1) = {-5, -5, 0, 1.0};
	Point(2) = {-5, 5, 0, 1.0};
	Point(3) = {15, -5, 0, 1.0};
	Point(4) = {15, 5, 0, 1.0};
	
For the circle, we need 3 points to trace a circle arc:

.. code-block::

	Point(5) = {0, 0, 0, 1.0};
	Point(6) = {-1, 0, 0, 1.0};
	Point(7) = {1, 0, 0, 1.0};
	
When the points are generated, we can close the shape with lines. For straight lines, the attribut ``Line`` is used. For the rectangle:

.. code-block::

	// Line(id) = {<id of start point>, <id of end point>}
	Line(1) = {1, 2};
	Line(2) = {2, 4};
	Line(3) = {4, 3};
	Line(4) = {3, 1};
	
And for the circle arcs:

.. code-block::

	// Circle(id) = {<id of start point>, <id of middle point>, <id of end point>}
	Circle(5) = {6, 5, 7};
	Circle(6) = {7, 5, 6};

.. warning::
	``Circle`` are also considered as ``Line``, so the ID needs to differ.
	
Shapes are now done, but gmsh does not recognize them as closed shapes or surfaces. The next steps consists to close the lines with ``Curve Loop`` and then create a surface out of it with ``Plane Surface``. Because we need to mesh inside the domain and all around the sphere, the final plane surface will be delimited by the curve loops of both the rectangle and the circle.

.. code-block::

	// Curve Loop(id) = {<id of line>, ...}
	Curve Loop(1) = {1, 2, 3, 4};
	Curve Loop(2) = {5, 6};
	// Plane Surface(id) = {<id of curve loop>, ...}
	Plane Surface(1) = {1, 2};
	
.. note::
	All the lines of code can be directly made with the GUI of gmsh with some clicks and keyboards shortcuts. In the ``Elementary entities`` section, you can ``Edit script`` to go back to the code and ``Reload script`` to visualize .geo file on the gmsh GUI.

""""""""""""""""""""""""""
OpenCASCADE kernel
""""""""""""""""""""""""""
In the gmsh geometry section of the gmsh GUI (see ``Geometry > Elementary entities > Add``), you can add directly multiple 2D or 3D common geometries with a simple click thanks to OpenCASCADE kernel. After selecting the geometry that you need, gmsh will automatically open a window where you can easily set the characteristic lenghts of the geometry. By doing so, the .geo file will be updated.
	
To represent a circle in a rectangle to simulate the flow around a sphere in 2D, we will first select the ``Disk`` geometry available with OpenCASCADE, set the radius to 1 (for X and Y) and center it at :math:`(x,y)=(0,0)`. Then, we will select the ``Rectangle`` geometry, set the length (DX) to 20, the width (DY) to 10 and the left bottom corner to :math:`(x,y)=(-5,-5)`.
    
.. note::
	If you click on ``Edit script``, you will see that the OpenCASCADE kernel has been added to the code as ``SetFactory("OpenCASCADE");``. Also, the rectangle is set with ``Rectangle(1) = {-5, -5, 0, 20, 10, 0};`` and the circle with ``Disk(2) = {0, 0, 0, 1, 1};``.
	
The ``Disk`` and ``Rectangle`` are already considered as surfaces in gmsh, so no need to pass from points, to curves and then surfaces. It is important to remark that the disk surface is overlapping the rectangle surface. The disk surface needs to be removed from the domain in order to simulate a flow around a sphere. OpenCASCADE offers boolean operations to do it. ``Intersection``, ``Union`` and ``Difference`` are available with the kernel and only works with OpenCASCADE geometries, such as ``Disk`` and ``Rectangle``. The code to use in order to do the difference between the rectangle and the disk is the following:

.. code-block::
	
	// BooleanDifference{ Surface{<id of surface to keep>}; Delete; }{ Surface{<id of surface to remove>}; Delete; }
	BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }
	
.. hint::
	If you do not remember what is the ID of a point, line, surface or volume, you can go to ``Tools > Options``, enable the labels boxes and choose the label type to make appear the ID on the GUI.

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
After generating your geometry, it is essential to set physical groups for boundary conditions identification that are compatible to the prm file.

.. warning::
  This step is essential for the mesh compatibility with Lethe.

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


