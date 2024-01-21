=============================
Floating Mesh (Solid Objects)
=============================

Floating meshes (solid objects) are finite (limited) auxiliary objects that can be stationary or moving. Rotating impellers, sliding surfaces, and finite stoppers are examples of floating meshes. The main differences between floating meshes and floating walls are:
1. floating wall is a surface, while a floating mesh can be a volume or a plane
2. floating wall is infinite, while floating mesh is finite
3. floating wall is stationary while floating mesh may be stationary or moving.

.. note:: 
    At the moment, solid objects (floating meshes) in Lethe have to be defined using triangular (simplex) meshes. Only triangular 2D meshes of 3D surfaces in 3D DEM simulations are presently supported. Quadrilateral 2D meshes of 3D surfaces and 1D mesh of 2D surfaces are not supported at the moment.

This subsection explains the solid objects (floating meshes) information. First of all, the ``number of solids`` is specified. Then, for each solid object, we need a ``mesh`` subsection. In these subsections, the ``type``, ``initial refinement``, and other properties of the objects are defined. Note that ``simplex`` must be ``true`` for all the objects since the solid objects are defined using simplex meshes. For more information on mesh subsection, visit `CFD mesh <https://lethe-cfd.github.io/lethe/documentation/parameters/cfd/mesh.html>`_

.. code-block:: text

 subsection solid objects
  # Total number of floating walls
    set number of solids  = 1
      subsection solid object 0
        subsection mesh
            set type                   = gmsh
            set file name              = none
            set simplex                = true
            set initial refinement     = 0
            set initial rotation axis  = 1, 0, 0
            set initial rotation angle = 0
            set initial translation    = 0, 0, 0

        end
    
        subsection translational velocity
            set Function expression = 0 ; 0 ; 0
        end

        subsection angular velocity
            set Function expression = 0 ; 0 ; 0
        end 

        subsection center of rotation
            # X COR
            set x = 0
            # Y COR
            set y = 0
            # Z COR
            set z = 0
        end

        set output solid object = true

    end
 end

* The ``number of solids`` parameter defines the total number of floating meshes we wish to use during the simulation.

* For each floating mesh, we need a separate subsection (for instance 	``subsection solid object 0``) in which the information of the floating mesh (``type``, ``file name``, ``initial refinement``, ``initial translation``, ``initial rotation axis``, ``initial rotation angle``  ) is defined. Each component of the ``initial translation`` and of the ``initial rotation axis`` parameters represent the ``x``, ``y`` and ``z`` axis. The rotation is apply before the translation.

* In the subsection ``translational velocity``, we define the translational velocity of the floating mesh.

* In the subsection ``angular velocity``, we define the angular velocity of the floating mesh.

* In the subsection ``center of rotation``, we define the center of rotation of the solid object (in the case of rotational motion).

* The ``output solid object`` defines if we want an output file to be generated for the solid object at every output time step.

