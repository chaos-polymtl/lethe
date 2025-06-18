=====================================
Turbulent Flow Around a Cylinder
=====================================

This example showcases the turbulent flow around a Cylinder at :math:`Re=3900`. . It features the matrix-free solver (``lethe-fluid-matrix-free``) which is more computationally efficient when solving problems using high-order elements and fine meshes. It also demonstrates the usage of static box refinement and boundary condition refinement to statically refine a mesh. 

---------
Features
---------

- Solvers: ``lethe-fluid-matrix-free`` (with Q2-Q2 or Q3-Q3)
- Transient problem using ``bdf2`` time integrator
- Static mesh refinement using the ``box refinement`` feature :doc:`../../../parameters/cfd/box_refinement`

---------------------------
Files Used in This Example
---------------------------

All files mentioned below are located in the example's folder (``examples/incompressible-flow/3d-turbulent-flow-around-a-cylinder``).

- Parameter file: ``turbulent-cylinder.prm``
- Postprocessing Python scripts: ``postprocess-cylinder.py`` and ``cylinder-functions.py``

------------------------
Description of the Case
------------------------

The flow around bluff bodies such as a cylinder is quite complicated and it often used as a benchmark problem for CFD. Such flow typically involved boundary-layer seperatoions, flow-regime transition, transition to turbulence, vortex shedding and coherent structures. If the body is symmetric, as is the case for a cylinder,  the wake usually exhibits self-induced periodicity from vortices being shed from alternate sides of the body, generating fluctuating forcs on the body. Taylor-Couette flow occurs in the annular space between two coaxial cylinders with different angular velocities. In this example, we study the flow around a cylinder at a Reynolds number of 3900 is considered to be in the subcritical turbulent regime.

This example is a canonical benchmark for LES, as explained in the book by Grinstein, Margolin and Rider [#wang2021]_. It also showcases the capabilities of Lethe to statically refine the mesh athe beggining of the simulation using user-defined box refinement. The mesh is refined in the vicinity of the cylinder to capture the boundary layer and the wake region more accurately.

The simulation set-up as well as the boundary ids are illustrated in the following figure:


.. image:: images/3d_cylinder_perspective_schematic.png
  :alt: The geometry and surface ID
  :align: center
  :name: geometry
  :height: 9cm

--------------
Parameter File
--------------

Mesh
~~~~

The ``mesh`` subsection specifies the computational grid. We use a custom mesh generated using the deal.II library's `GridGenerator <https://www.dealii.org/current/doxygen/deal.II/namespaceGridGenerator.html>`_ to create the flow using high-order elements.

.. code-block:: text
  
  subsection mesh
    set type               = dealii
    set grid type          = custom_channel_with_cylinder
    set grid arguments     = 25 : 8 : 52 : 4.71238898038 : 4 : 0.75 : 5 : 1:  false: true
    set initial boundary refinement = 1 
    set boundaries refined = 2
  end

.. note::

  This mesh generator is only present in the 9.7 version of the deal.II library.

Box refinement
~~~~~~~~~~~~~~~~
The ``box refinement`` subsection allows us to refine the mesh in a specific region of the domain. In this case, we refine the mesh around the cylinder to capture the boundary layer and wake region more accurately. The box refinement is defined by its center, size, and refinement level. The mesh is refined by the box refinement twice by setting ``ìnitial refinement = 2``.

.. code-block:: text

  subsection box refinement
    subsection mesh
      set type               = dealii
      set grid type          = subdivided_hyper_rectangle
      set grid arguments     = 1,1,1 : -2, -3, -1 : 52,3,5 : false
      set initial refinement = 0
    end
    set initial refinement = 2
  end

Boundary Conditions
~~~~~~~~~~~~~~~~~~~

The ``boundary conditions`` subsection establishes the constraints on different parts of the domain: 

.. code-block:: text

  subsection boundary conditions
    set number = 6
    subsection bc 0
      set type = function
      subsection u
        set Function expression = 1
      end
      subsection v
        set Function expression = 0
      end
      subsection w
        set Function expression = 0
      end
    end
    subsection bc 1
     set type = outlet
     set beta = 1
    end
    subsection bc 2
      set type = noslip
    end
    subsection bc 3
      set type = slip
    end
    subsection bc 4
      set type = slip
    end
    subsection bc 5
      set type = periodic
      set periodic_id = 6
      set periodic_direction = 2
    end
  end


Periodic boundary conditions are applied to the front (``id=5``) and the back (``id=6``) of the domain to mimic an infinite domain in the axial direction.

Physical Properties
~~~~~~~~~~~~~~~~~~~

In the present case, the Reynolds number is defined as: :math:`Re = \frac{UD}{\nu}`. Since we set the values of :math:`U` and :math:`D`, the Reynold number of 3900 can be set solely using the kinematic viscosity: 


.. code-block:: text

  subsection physical properties
    set number of fluids = 1
    subsection fluid 0
      set kinematic viscosity = 2.5641025e-04
    end
  end


FEM Interpolation
~~~~~~~~~~~~~~~~~

The results obtained for the turbulent flow around a cylinder are highly mesh and order dependent. The present eaxmples consider both :math:`Q_1Q_1` and :math:`Q_2Q_2` elements.

.. code-block:: text

    subsection FEM
      set velocity order = 1  #2 for Q3
      set pressure order = 1  #2 for Q3
    end

Forces
~~~~~~

The ``forces`` subsection controls the postprocessing of the torque and the forces acting on the boundaries of the domain: 

.. code-block:: text

  subsection forces
    set verbosity             = verbose
    set calculate force       = true
    set force name            = force
    set output precision      = 10
    set output frequency      = 10
  end

By setting ``calculate force = true``, the calculation of the force resulting from the fluid dynamics physics on every boundary of the domain is automatically calculated. 


Post-processing
~~~~~~~~~~~~~~~

.. code-block:: text

  subsection post-processing
    set calculate average velocities      = true
    set initial time for average velocity = 25
  end

To monitor the average velocity and pressure, we set ``calculate average velocities = true``. The average velocity is computed starting from the time step specified by ``initial time for average velocity = 25``. This allows us to focus on the statistically steady state of the flow. 

Simulation Control
~~~~~~~~~~~~~~~~~~

The ``simulation control`` subsection controls the flow of the simulation. To maximize the temporal accuracy of the simulation, we use a second-order ``bdf2`` scheme. Results are written every 500 time-steps. 

.. code-block:: text

  subsection simulation control
    set method           = bdf2
    set output name      = cylinder-Re3900
    set output path      = ./output/
    set time end         = 200                               
    set adapt            = true
    set max cfl          = 1
    set time step        = 0.002
    set output frequency = 500
  end


----------------------
Running the Simulation
----------------------

Launching the simulation is as simple as specifying the executable name and the parameter file. Assuming that the ``lethe-fluid-matrix-free`` executable are within your path, the matrix-free simulation can be launched by typing:

.. code-block:: text
  :class: copy-button

  mpirun -np n_proc lethe-fluid-matrix-free turbulent-cylinder.prm 

and choosing the number of processes ``n_proc`` according to the resources you have available.

.. note::

  THe simulation takes approximatively 10 hours on 16 cores of a AMD Ryzen 9 7950X 16-Core Processor.

----------------------
Results and Discussion
----------------------

In the following, results obtained with a box refinement of 2 as well a 3 and using 

The turbulent flow around a cylinder is quite complex. The following animation displays the evolution of the velocity magnitude on a slice of the domain as a function.


.. list-table::
   :widths: 20 20 20 20
   :header-rows: 1

   * - Study
     - :math:`C_D`
     - :math:`C_L`
     - :math:`S_t`
   * - Lethe example
     - 1.396 :math:`\pm` 0.048
     - -0.003 :math:`\pm` 0.72
     - 0.2
   * - Lethe Sharp [#barbeau2022]_
     - 1.395 :math:`\pm` 0.047
     - :math:`\pm` 0.71
     - 0.2
   * - Braza et al. [#braza1986]_
     - 1.400 :math:`\pm` 0.050
     - :math:`\pm` 0.75
     - 0.2


The flow patterns generated by the Taylor-Couette flow are quite complex. The following animation displays the evolution of velocity magnitude on the radial-vertical plane (left) and the Q-criterion iso-contours (right), illustrating the vortical structure as the vortex breaks down and generates smaller structures.

..
  +----------------------------------------------------------------------------------------------------------------------------------------------------+
  | .. raw:: html                                                                                                                                      |
  |                                                                                                                                                    |
  |    <iframe width="800" height="400" src="https://www.youtube.com/embed/bRa04yMDsXo?si=Q1ppAuakIsrNwFlw"  frameborder="0" allowfullscreen></iframe> |
  |                                                                                                                                                    |
  +----------------------------------------------------------------------------------------------------------------------------------------------------+

..
  +-------------------------------------------------------------------------------------------------------------------+
  |  .. figure:: images/enstrophy_comparison_Q3Q3_942k.png                                                            |
  |     :width: 620                                                                                                   |
  |                                                                                                                   |
  +-------------------------------------------------------------------------------------------------------------------+


----------------------------
Possibilities for Extension
----------------------------

- This case offers numerous options for postprocessing. Consider exploring alternative quantities such as vorticity and pressure and use the results to generate interesting animations. Feel free to share them with us!
- It could also be interesting to explore this case with an even higher Reynolds number

------------
References
------------

.. [#wang2021] \Z. J. Wang and E. Jourdan, “Benchmark for scale-resolving simulation with curved walls: the Taylor Couette flow,” Advances in Aerodynamics, vol. 3, no. 1, Jun. 2021, doi: `10.1186/s42774-021-00071-0 <https://doi.org/10.1186/s42774-021-00071-0>`_\.

.. [#wikipedia2024] \“Taylor–Couette flow,” *Wikipedia*. Feb. 15, 2024. Available: https://en.wikipedia.org/wiki/Taylor%E2%80%93Couette_flow\.
