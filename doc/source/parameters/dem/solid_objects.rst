=============================
Solid Objects
=============================

Solid objects are finite auxiliary objects that can be stationary or in motion. Rotating impellers, sliding surfaces, and finite stoppers are examples of solid objects. The main differences between them and floating walls are:
1. floating wall is a plane, while a solid object can be either a volume or a surface
2. floating wall is infinite, while solid object is finite
3. floating wall is stationary while solid object may be stationary or moving.

.. note:: 
    At the moment, solid objects in Lethe have to be defined using triangular (simplex) meshes. Only triangular 2D meshes of 3D surfaces in simulations are currently supported, meaning that any other meshes will not work.

This subsection explains the solid objects information. First of all, the ``number of solids`` is specified. Then, for each solid object, we need a ``mesh`` subsection. In these subsections, the ``type``, ``initial refinement``, and other properties of the objects are defined. Note that currently, only ``solid surfaces``, meaning ``dim=2, spacedim=3`` meshes, are usable. The ``simplex`` parameter must be set to ``true`` for all ``solid surfaces``. For more information on mesh subsection, visit `CFD mesh <https://chaos-polymtl.github.io/lethe/documentation/parameters/cfd/mesh.html>`_

.. code-block:: text

 subsection solid objects
   subsection solid surfaces
     # Total number of solid surfaces
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

       set center of rotation    = 0., 0., 0.
       set output solid object = true

       # Choices are adiabatic|isothermal
       set thermal boundary type = adiabatic
       # If type = isothermal
       subsection temperature
         set Function expression = 0
       end
     end
   end
 end

* The ``number of solids`` parameter defines the total number of solid surfaces used during the simulation.

* For each solid object, we need a separate subsections (for instance 	``subsection solid object 0``) in which the information of the solid surface 0 (``type``, ``file name``, ``initial refinement``, ``initial translation``, ``initial rotation axis``, ``initial rotation angle``) is defined. Each component of the ``initial translation`` and of the ``initial rotation axis`` parameters represent the ``x``, ``y`` and ``z`` axes. The rotation is applied before the translation.

* In the subsection ``translational velocity``, we define the translational velocity function of the solid object.

* In the subsection ``angular velocity``, we define the angular velocity function of the solid object.

* The ``center of rotation`` defines the center of rotation of the solid object (in the case of rotational motion) in the ``x``, ``y`` and ``z`` directions. This center of rotation will move according to the ``translational velocity`` parameter.

* The ``output solid object`` defines if we want an output file to be generated for the solid object at every output time step.

* The ``thermal boundary type`` defines whether the solid object should be adiabatic or isothermal when launching a multiphysic DEM simulation.

* In the subsection ``temperature``, we define the temperature of the solid object as a function of time. It is used only if the solid object is not adiabatic.