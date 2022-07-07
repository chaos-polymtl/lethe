===============================
Introduction on how to use GMSH
===============================

--------------------------
Installation
--------------------------

--------------------------
Geometry
--------------------------

This guide uses mainly this simple geometry:

.. image:: images/geo.png
    :alt: The geometry of a circle in a rectangle
    :align: center
    :name: geometry

""""""""""""""""""""""""""""""""""
Create and Edit the ``.geo`` file
""""""""""""""""""""""""""""""""""

The geometry is written in a ``.geo`` file:

1. Start ``GMSH``
2. ``File > New...``
3. Either use the folders navigator or directly specify the path and filename in ``Filename`` 
    for example: ``/home/<user>/Documents/Mesh/<filename>.geo``
4. Specify the Kernel to use: :ref:`built-in kernel` or :ref:`opencascade kernel`

You can then open the created ``.geo`` file in a simple text editor, either by:

* opening the file in your folder, or
* ``Left pannel: Modules > Geometry > Edit script``

.. warning::
	After each modification in the ``.geo`` file, save it and load the modifications in ``GMSH``:
	    ``Left pannel: Modules > Geometry > Reload script``

.. _built-in kernel:

""""""""""""""""""""""""""
Built-in kernel
""""""""""""""""""""""""""
It is quite easy to create a ``.geo`` file directly by coding line by line the geometry that you want to produce.

.. tip::
	Using a Built-in kernel gives you full control over the geometry creation, which is necessary to create a :ref:`structured mesh` (transfinite) mesh. However, it usually takes more time than using functions available in :ref:`opencascade kernel`.

1. Set the points that delimits the shape being traced. The attribut ``Point`` is essential to set a point of the shape and to specify its coordinates. 
	* For the rectangle, we will have:

	.. code-block::

		// Point(<id>) = {<X>, <Y>, <Z>, <Prescribed mesh size at point>}
		Point(1) = {-5, -5, 0, 1.0};
		Point(2) = {-5, 5, 0, 1.0};
		Point(3) = {15, -5, 0, 1.0};
		Point(4) = {15, 5, 0, 1.0};

	.. hint::
		See :ref:`mesh` for more information about the ``<Prescribed mesh size at point>``.
		
	* For the circle, we need 3 points to trace a circle arc:

	.. code-block::

		Point(5) = {0, 0, 0, 1.0};
		Point(6) = {-1, 0, 0, 1.0};
		Point(7) = {1, 0, 0, 1.0};
	
2. Close the shape with lines. 
	* For the rectangle (straight lines):

	.. code-block::

		// Line(id) = {<id of start point>, <id of end point>}
		Line(1) = {1, 2};
		Line(2) = {2, 4};
		Line(3) = {4, 3};
		Line(4) = {3, 1};
	
	* For the circle arcs:

	.. code-block::

		// Circle(id) = {<id of start point>, <id of middle point>, <id of end point>}
		Circle(5) = {6, 5, 7};
		Circle(6) = {7, 5, 6};

	.. warning::
		``Circle`` are also considered as ``Line``, so the ``<id>`` needs to differ from the ``Line`` ones.

3. Close the lines with ``Curve Loop`` and then create a surface out of it with ``Plane Surface``. The final plane surface will be delimited by the curve loops of both the rectangle and the circle.

.. code-block::

	// Curve Loop(id) = {<id of line>, ...}
	Curve Loop(1) = {1, 2, 3, 4};
	Curve Loop(2) = {5, 6};
	// Plane Surface(id) = {<id of curve loop>, ...}
	Plane Surface(1) = {1, 2};
	
.. tip::
	All the lines of code can be directly made with the GUI of gmsh with some clicks and keyboards shortcuts.

.. _opencascade kernel:

""""""""""""""""""""""""""
OpenCASCADE kernel
""""""""""""""""""""""""""
In the GMSH geometry section of the GMSH GUI (see ``Left pannel: Modules > Geometry > Elementary entities > Add``), you can add directly multiple 2D or 3D common geometries with a simple click thanks to OpenCASCADE kernel. GMSH will automatically open a window where you can easily set the characteristic lenghts of the geometry, and update the ``.geo`` file.

.. warning::
	Always save your ``.geo`` file in your text editor before modifying it through the GMSH GUI. If you modified the ``.geo`` file without saving it, GMSH will not update it. 

For our example (circle in a rectangle in 2D):

1. Select the ``Disk`` geometry available with OpenCASCADE, set the radius to 1 (for X and Y) and center it at :math:`(x,y)=(0,0)`. 
2. Select the ``Rectangle`` geometry, set the length (``DX``) to :math:`20`, the width (``DY``) to :math:`10` and the left bottom corner to :math:`(x,y)=(-5,-5)`.
    
.. note::
	If you click on ``Edit script``, you will see that the OpenCASCADE kernel has been added to the code as ``SetFactory("OpenCASCADE");``. 

	The rectangle is set with ``Rectangle(1) = {-5, -5, 0, 20, 10, 0};`` and the circle with ``Disk(2) = {0, 0, 0, 1, 1};``.
	
.. tip::
	The ``Disk`` and ``Rectangle`` are already considered as surfaces in gmsh, so no need to pass from points, to curves and then surfaces. 

3. Remove the disk surface from the rectangular domain, with OpenCASCADE boolean operation, either via the GUI (``Geometry > Elementary entities > Boolean``) or the code:

.. code-block::
	
	// BooleanDifference{ Surface{<id of surface to keep>}; Delete; }{ Surface{<id of surface to remove>}; Delete; }
	BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }
	
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

.. _mesh:

---------------------------
Mesh
---------------------------

""""""""""""""""""""""""""
Unstructured
""""""""""""""""""""""""""

Basic:

1. (optional) ``Tools > Options > Mesh`` and ``General`` panel, check ``Recombine all triangular meshes``: generate a quad mesh.
2. (optional) In the same panel, change ``Min/Max element size`` to have smaller/bigger elements, therefor a finer/coarser mesh.
3. ``Left pannel: Modules > Mesh > 2D`` or ``3D``: create the mesh
4. ``Tools > Statistics``: check that the mesh is generated appropriately (by default, triangles for 2D and hexahedra for 3D)
5. (optional) ``Left pannel: Modules > Mesh > Refine by splitting``: refine the mesh (beware, it takes more and more time for each refinement)
6. ``Left pannel: Modules > Mesh > Save``: save the mesh in a ``.msh`` file, to be used in Lethe (see :doc:`../../parameters/cfd/mesh`)

By following all the previous steps, the mesh generated looks like bellow.

.. image:: images/unstructured.png
    :alt: 2D mesh with quads
    :align: center
    :name: unstructured mesh
    
Attractors can also be used to refine the mesh towards specific edges or surfaces. In this example, attractors could be interesting if the mesh needs to be finer around the sphere. Attractors can only be added by code with the ``Field`` attribut.

.. code-block::

	// LcMax -                         /------------------
	//                               /
	//                             /
	//                           /
	// LcMin -o----------------/
	//        |                |       |
	//     Attractor       DistMin   DistMax

1. Set the attractor:

.. code-block::

	Field[1] = Attractor;
	
2. Specify where the refinement needs to be done (near the circle in this case):

.. code-block::

	Field[1].EdgesList = {5};
	
3. Set the minimum/maximum characteristic length and the minimum/maximum distance of the refinement:

.. code-block::

	Field[2] = Threshold;
	Field[2].IField = 1;
	Field[2].LcMin = 0.25;
	Field[2].LcMax = 1;
	Field[2].DistMin = 1;
	Field[2].DistMax = 2;
	
4. Apply the attractor:

.. code-block::

	Background Field = 2;
	
Here is the mesh generated with an attractor around the sphere:

.. image:: images/attractor.png
    :alt: The mesh generated with an attractor arround the sphere
    :align: center
    :name: attractor

.. _structured mesh:

""""""""""""""""""""""""""
Structured
""""""""""""""""""""""""""
.. warning::
	The ``.geo`` file must be built with the :ref:`built-in kernel`.

.. _tips:

--------------------------
Other tips
--------------------------

Use the ``Visibility`` options to get the ID of an element easily on the GUI: 
	* ``Tools > Options > Mesh > Tab: Visibility``
	* Check the adequate boxes (for example ``1D element labels`` for points, etc.) 
	* Choose the label type in the drop-down menu ``Label type`` (for example ``Elementary entity tag``).


