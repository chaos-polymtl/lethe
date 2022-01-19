======================================
3D Flow over periodic hills
======================================

This example is a a well-established benchmark for Computational Fluid Dynamics softwares known as the `periodic hills flow <https://kbwiki.ercoftac.org/w/index.php?title=Abstr:2D_Periodic_Hill_Flow>`_. It includes complex flow features such as the generation of an unsteady shear layer, recirculation, strong pressure gradients, attached and detached boundary layers and turbulence recycling due to the periodicity assumption. 

Features
---------

- Solver: ``gls_navier_stokes_3d`` (with Q1-Q1) 
- Transient problem

Location of the example
------------------------

- Parameter file: ``/examples/incompressible_flow/3d_periodic_hills/periodic_hills.prm``

Description of the case
-----------------------
In this case a well-defined flow passes over a series of hills which repeat along a channel in a periodic fashion as it can be seen in the following figure (taken from ERCOFTAC [1]):

.. image:: images/geometry_description.jpg
    :alt: The geometry
    :align: center
    :name: geometry_description

Parameter file
--------------

All the sections of the parameter file used in this case have been already explained in previous examples. However, for the sake of completeness, the important sections are briefly explained.

Simulation control
~~~~~~~~~~~~~~~~~~~

This section controls the flow of the simulation. 

.. code-block:: text

 # --------------------------------------------------
 # Simulation and IO Control
 #---------------------------------------------------
 subsection simulation control
   set method                       = bdf2
   set output name                  = periodic_hills_output_
   set time step                    = 0.1
   set output frequency             = 1000
   set output path                  = ./output/
   set time end                     = 1000
 end

The ``method`` parameter specifies the time-stepping scheme chosen for this simulation. In this case is set to ``bdf2`` that corresponds to a second-order backward difference implicit scheme. The ``output name`` and ``output path``  are parameters that specity the name and the path for the ``.pvd`` files of the simulation, while the ``output frequency`` is relative to the number of time steps in the total time. The ratio between the total time of the simulation ``1000`` and the time step chosen ``0.1`` is 10000. Therefore, if the output frequency is set to ``1000`` a total of 11 ``.pvd`` file will be obtained including an initial output file at time 0. 

.. warning:: It is important to remember that the output path folder, in this case ``output`` must exist before running the simulation.


Physical properties
~~~~~~~~~~~~~~~~~~~

The physical properties section is used to target a specific Reynolds number:

.. code-block:: text

    #---------------------------------------------------
    # Physical Properties
    #---------------------------------------------------
    subsection physical properties
        set kinematic viscosity        = 1.78571E-04 # Re = 5600
    end

Recall the definition of the Reynolds number:

.. math::
 Re = \frac{u_B h}{\nu}

To be able to control it we set the bulk velocity to 1 :math:`m s^{-1}` and the height of the hill :math:`u_B` to 1 :math:`m`, and adjust the kinematic viscosity accordingly. In that sense, the Reynolds number set for this example is equal to 5600. 

Mesh 
~~~~~

The mesh subsection specifies the computational grid:

.. code-block:: text

    #---------------------------------------------------
    # Mesh
    #---------------------------------------------------
    subsection mesh
        set type                       = periodic_hills
        set initial refinement         = 5
        set grid arguments 		       = 1;1;4;2;1
    end

The standard geometry is included in the Lethe code with 6 polynomials designing the curve of the hill. Therefore the ``type`` parameter is set to ``periodic_hills``. The ``initial refinement`` parameter does a uniform refinement in the triangulation. The ``grid argument`` parameter requires five values that correspond to ``spacing y lines; alpha; repetitions x; repetitions y; repetitions z``. A brief explanation of each of them is given here:

* ``spacing y lines``: the periodic hill grid is generated with equally spaced horizontal lines or with gradually spaced line which means more lines near upper and lower walls. This parameter may be set in the range from 0 to 1, with 0 meaning equally spaced lines and 1 the maximum of shifting lines.

* ``alpha``: While the geometry of the benchmark is fixed, it can be elongated to compare the behaviour of the flow with the slope. Actually, this parameter elongates the slopes of the geometry and the flat region has the same length. This parameter should be set in the range from 0.5 to 3, meaning 1 has no effect on the geometry.

* ``repetitions x, repetitions y, repetitions z``: To get cells with an aspect ratio different on the domain, different number of subdivisions are set, given by repetition in different coordinate directions. The minimum number of subdivisions in each direction is 1. 

The following image displays a coarse mesh for this example. It can be seen that the horizontal lines are shifted with the associated parameter to get more lines near walls. Here, repetitions for x, y, z allow to get more cells in x and y directions.

.. image:: images/mesh.png
    :alt: Mesh
    :align: center
    :name: mesh


Flow control
~~~~~~~~~~~~

Since the flow is periodic and a specific Reynolds number is targeted for the simulation, the flow has to be controlled at each time step. To allow flow control, the subsection flow control has to be enabled. 

.. code-block:: text

 #---------------------------------------------------
 # Flow control
 #---------------------------------------------------
 subsection flow control
     set enable                    = true
     set boundary id    		   = 0
     set volumetric flow rate      = -9.1575 # bulk velocity = -1
     set flow direction 		   = 0
     set initial beta		       = 7.66
     set verbosity                 = verbose
 end

First we set the ``enable`` parameter to ``true`` in order to control the flow. The boundary id given, `0`, corresponds to the flow inlet where we want to control the flow. The ``volumetric flow rate`` has to be negative if the flow goes in x positive direction or ``flow direction = 0``. Therefore we adjust this parameter so that we obtain a bulk velocity :math:`u_B` equals to 1. The ``initial beta`` parameter is a coefficient calculated at each time step that speeds up the convergence of the flow rate targeted.

.. tip:: A good method to find a reasonable initial beta us to test two or three different initial beta parameters, write down the given flow rate at the first time step in the simulation and do a regression. The correlation is linear and giving a proper value will greatly speed up the convergence. 

Post-processing
~~~~~~~~~~~~~~~~~~~

The post-processing subsection allows the calculation of different quantities:

.. code-block:: text

 #---------------------------------------------------
 # Post-Processing
 #---------------------------------------------------
 subsection post-processing
     set calculate average velocities    = true
     set initial time 		             = 207
 end

In this example, we enable the calculation of average velocities through the parameter ``calculate average velocities`` after certain time of the simulation. In this case, this time is set to ``207`` as we allow for the flow to achieve some stability. The results of this calculations will be available in the ``.pvd`` files when open with a visualization software. 

Boundary conditions
~~~~~~~~~~~~~~~~~~~~
In this section, we specify the boundary conditions taking into account the IDs presented in the following scheme:

.. image:: images/boundary_conditions.png
    :alt: bcs
    :align: center
    :name: boundary_conditions

.. code-block:: text

 # --------------------------------------------------
 # Boundary Conditions
 #---------------------------------------------------
 subsection boundary conditions
   set number                      = 4
     subsection bc 0
         set type                  = periodic
         set id                    = 0
         set periodic_id           = 1
         set periodic_direction    = 0
     end
     subsection bc 1
         set id                    = 2
         set type                  = noslip
     end
     subsection bc 2
         set id                    = 3
         set type                  = noslip
     end
     subsection bc 3
         set type                  = periodic
         set id                    = 4
         set periodic_id           = 5
         set periodic_direction    = 2
     end 
 end

As it can be seen first a ``periodic`` boundary condition is set for both the inlet id ``0`` and outlet id ``1`` of the flow. For the bottom and top walls we set a ``noslip`` boundary conditions, while for the side walls id ``4`` and ``5`` we consider periodic boundary conditions, which allow to represent the bulk flow of the channel. All the boundary conditions are set to represent the actual benchmark case. 

FEM
~~~
The FEM subsection specifies the order of the elements used for both velocity and pressure.

.. code-block:: text

 #---------------------------------------------------
 # FEM
 #---------------------------------------------------
 subsection FEM
     set velocity order            = 1
     set pressure order            = 1
 end

This example uses Q1-Q1 for the periodic hills simulation,

Non-Linear Solver Control
~~~~~~~~~~~~~~~~~~~~~~~~~

The non-linear solver control section allows us to choose a method suitable for the problem that we are solving:

.. code-block:: text

 # --------------------------------------------------
 # Non-Linear Solver Control
 #---------------------------------------------------
 subsection non-linear solver
   set solver                      = inexact_newton
   set tolerance                   = 1e-5
   set max iterations              = 10
   set verbosity                   = verbose
 end

In this case, we use the ``inexact_newton`` method that reuses the Jacobian matrix between interations. This is a known strategy to reduce the cost of reassembling the Jacobian in every iteration. 

Running the simulation
----------------------
Launching the simulation is as simple as specifying the executable name and the parameter file. Assuming that the ``gls_navier_stokes_3d`` executable is within your path, the simulation can be launched by typing:

.. code-block:: text

  gls_navier_stokes_3d periodic_hills.prm

Lethe will generate a number of files. The most important one bears the extension ``.pvd``. It can be read by popular visualization programs such as `Paraview <https://www.paraview.org/>`_. 

Due to the complexity of this example we recommend that you run this example using a cluster or supercomputer if available. For this it is necessary to add ``mpirun -np 8`` command at the beginning of the line. The number of processes must be adjusted according to the machine. If you want to run this in a normal desktop we recommend that you set ``time end`` to ``5.0``, in order to observe how the simulation behaves initially.

Results
-------



Possibilities for extension
----------------------------

- **High-order elements**: It would be interesting to observe the effect of high-order elements in the simulation of the periodic hills flow. For example, Q2-Q2 elements. The only part of the parameter file that would need to be change would be the ``FEM`` section.

- **High Reynolds numbers**: The example can be run at higher Reynolds numbers. In fact, one can find experimental and numerical results in the literature for Reynolds numbers equal to 10600 or 37000.

- **Dynamic mesh adaptation:** This can be a good example to observe how dynamic mesh adaptation works in a complex transient turbulent simulation. 

References
----------
[1] ERCOFTAC. File: hill3d.jpg. 2010. URL https://www.kbwiki.ercoftac.org/w/index.580php?title=File:Hill3d.jpg#filelinks.

