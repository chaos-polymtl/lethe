======================================
Mortar method
======================================

This example showcases the use of the mortar method in Lethe to solve a 2D incompressible flow problem with a rotating impeller.

---------
Features
---------

- Solver: ``lethe-fluid`` (with Q1-Q1) 
- Transient problem using bdf1 time integrator
- Shows the usage of different parameters related to the mortar method
- Explains how to set boundary conditions in a gmsh file before using the mortar method 



----------------------------
Files Used in This Example
----------------------------

All files mentioned below are located in the example's folder (``examples/incompressible-flow/2d-mortar-method``).

- Geometry file: ``full-geometry.geo``
- Rotor mesh file: ``rotor.msh``
- Stator mesh file: ``stator.msh``
- Parameter file: ``mortar.prm``



-----------------------
Description of the Case
-----------------------
We simulate the flow around a rotating impeller in a cylindrical domain with the presence of a buffet. This domain is composed of the rotor part, which contains the impeller and is rotating, and the stationary stator part. These two geometries are attached with mortar elements, which maintain the continuity of the solution between the two domains.

The geometry of the problem is illustrated in the following figure:

.. image:: images/geometry-description.png
    :alt: The geometry
    :align: center
    :name: geometry_description



--------------
Parameter File
--------------

Simulation Control
~~~~~~~~~~~~~~~~~~

A 1st order backward differentiation is used for the time integration (``bdf1``), for a 0.01 second ``time step`` and a 10 seconds simulation. Outputs are saved every 5 steps, as indicated by the ``output frequency`` parameter.

.. code-block:: text

  subsection simulation control
    set method           = bdf1
    set time step        = 0.01
    set time end         = 10
    set output frequency = 10
    set output path      = ./output/
    set output name      = mortar-output
  end

FEM
~~~

For this example, Q1 elements are used for both velocity and pressure. The bubble enrichment function is also enabled for the velocity field to increase the accuracy and stability of the solution. 

.. code-block:: text

  subsection FEM
    set enable bubble function velocity = true
    set velocity order                  = 1
    set pressure order                  = 1
  end

Mesh
~~~~

When using the mortar method, only the static mesh needs to be included in the ``mesh`` subsection. In this example, ``statorPW.gmsh`` only contains the meshed stator shown in the above case description.

.. code-block:: text

  subsection mesh
    set type               = gmsh
    set file name          = statorPW.msh
    set initial refinement = 2
  end

Mortar
~~~~~~

The ``mortar`` subsection specifies all the parameters required to simulate a rotor-stator geometry using mortar elements to attach the two meshes. The ``mesh`` subsection embedded in the ``mortar`` subsection refers to the rotor domain and contains the same exact parameters as the above stator ``mesh`` subsection. Is is important that the nodes of the rotor and stator meshes initially match at the interface for the mortar method to work; uneven meshes may lead to nodes being ignored at the mortar interface.
After the inputs for the rotor geometry have been given, both ids at the rotor-stator interface need to be given. All the inputs related to the rotation of the rotor domain can now be entered :
- The ``center of rotation`` represents the coordinates around which the rotor domain will be rotating.
- The ``rotor rotation angle`` subsection contains the expression for the rotation angle, which can either be constant or time dependant.
- The ``rotor angular velocity`` subsection contains the expression for the angular velocity or the rotor and needs to correspond to the time derivative of the ``rotor rotation angle``.
We can now ``set verbosity`` to ``verbose`` to print information about the rotor rotation at every iteration. Lastly, we set the ``radius tolerance`` to 1e-6; every iteration, the radial distance between the center of rotation and every stator and rotor node is calculated to verify that the difference between the maximum and minimum value is lesser than the tolerance.

.. note::
    When creating the .msh files from the .geo files, the boundary ids need to be unique for the rotor and stator domains. For example, if the boundary id for the rotor-stator interface is 3 in the rotor mesh, it should be 4 in the stator mesh. This way, both boundaries can be specified in the ``mortar`` subsection without any conflict. As of today, a manual modification of the .msh files may be required to achieve this.

.. code-block:: text

  subsection mortar
    set enable = true
    subsection mesh
      set type               = gmsh
      set file name          = rotorPW.msh
    end
    set rotor boundary id  = 3
    set stator boundary id = 4
    set center of rotation = 0, 0
    subsection rotor rotation angle
      set Function expression = 1*t
    end
    subsection rotor angular velocity
      set Function expression = 1
    end
    set verbosity = verbose
    set radius tolerance = 1e-6
  end

Boundary Conditions
~~~~~~~~~~~~~~~~~~~

Both the outside wall (ID 5) and the buffet (ID 6) have a no-slip Dirichlet boundary condition. The mortar boundaries (ID 3 and 4) are already defined in the ``mortar`` subsection, so a ``none`` boundary condition is assigned to them in this subsection. Finally, the velocity of the fluid at the impeller boundary (ID 2) is defined as a function of the coordinates, since the impeller is rotating.

.. code-block:: text

  subsection boundary conditions
    set number = 3
    subsection bc 0
      set type               = noslip
      set id                 = 5,6
    end
    subsection bc 1
      set type               = none
      set id                 = 3,4
    end
    subsection bc 2
      set type = function
      set id   = 2
      subsection u
        set Function expression = -y
      end
      subsection v
        set Function expression = x
      end
    end
  end

Manifolds
~~~~~~~~~

Since the mortar boundary and the outside wall are both circular, we need to set-up manifolds to assure the circular shape is not lost when refining the mesh. For more details about manifolds parameters, 

.. code-block:: text

  subsection manifolds
    set number = 3
    subsection manifold 0
      set id = 3
      set type = spherical
      set point coordinates = 0, 0
    end
    subsection manifold 1
      set id = 4
      set type = spherical
      set point coordinates = 0, 0
    end
    subsection manifold 2
      set id = 5
      set type = spherical
      set point coordinates = 0, 0
    end
  end


Linear Solver
~~~~~~~~~~~~~

The only modification made in the linear solver is the use of the AMG preconditioner. Using the surface LIC visualization method, we obtain the following animation :

.. code-block:: text

subsection linear solver
  subsection fluid dynamics
    set verbosity                             = quiet
    set method                                = gmres
    set preconditioner                        = amg
  end
end



----------------------
Running the Simulation
----------------------

Assuming that the ``lethe-fluid`` executable is within your path, the simulation can be launched by typing

.. code-block:: text
  :class: copy-button

  lethe-fluid mortar.prm

Lethe will generate a number of files. The most important one bears the extension ``.pvd``. It can be read by visualization programs such as `Paraview <https://www.paraview.org/>`_.



----------------------
Results and Discussion
----------------------

Using Paraview, we can visualize the results of the simulation with the ``.pvd`` output file. Using the surface LIC visualization method, we obtain the following animation : 

.. raw:: html

    <iframe width="560" height="315" src="https://www.youtube.com/embed/rt6PAvgMkio" frameborder="0" allowfullscreen></iframe>

.. raw:: html

   <iframe width="560" height="315" src="https://www.youtube.com/embed/jvaT76qBBTs" frameborder="0" allowfullscreen></iframe>

As we can see, the presence of the buffet creates recirculation zones in the flow, which are visible in the animations. Continuity of the solution is maintained at the mortar interface, as expected. 


----------------------------
Possibilities for Extension
----------------------------
- Vary the Reynolds number to see the effect on the flow.
- Run the simulation for different impeller and/or buffet geometries.
- Make the case three-dimensional by using a cylindrical mortar interface.


----------
Reference
----------

.. [#larson2013] \M. G. Larson, F. Bengzon. *The Finite Element Method: Theory, Implementation, and Applications*, vol. 10. Springer Berlin, Heidelberg, 2013. 
