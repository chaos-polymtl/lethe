==================================
Gas-Solid Spouted Bed
==================================


----------------------------------
Features
----------------------------------

- Solvers: ``lethe-particles`` and ``lethe-fluid-particles``
- Three-dimensional problem
- Displays the selection of models and physical properties
- Simulates a solid-gas spouted bed


---------------------------
Files Used in This Example
---------------------------

Both files mentioned below are located in the example's folder (``examples/unresolved-cfd-dem/gas-solid-spouted-bed``).

- Parameter file for particle generation and packing: ``packing-particles.prm``
- Parameter file for CFD-DEM simulation of the spouted bed: ``gas-solid-spouted-bed.prm``


-----------------------
Description of the Case
-----------------------

This example simulates the spouting of spherical particles in air. First, we use ``lethe-particles`` to fill the bed with particles. We enable check-pointing in order to write the DEM checkpoint files which will be used as the starting point of the CFD-DEM simulation. Then, we use the ``lethe-fluid-particles`` solver within Lethe to simulate the spouting of the particles by initially reading the checkpoint files from the DEM simulation.


-------------------
DEM Parameter File
-------------------

All parameter subsections are described in the `parameter section <../../../parameters/parameters.html>`_ of the documentation.

To set-up the rectangular spouted bed case, we first fill the bed with particles.

We introduce the different sections of the parameter file ``packing-particles.prm`` needed to run this simulation.

Mesh
~~~~~

In this example, we are simulating a rectangular spouted bed. In order to ensure that the flow is developed at the inlet of the bed and that the void fraction does not interfere with the Dirichlet condition of the inlet velocity imposed at the inlet of the spout, we introduce the flow through a channel that is connected to the inlet of the bed. The geometry of the bed was created using `GMSH <https://gmsh.info/>`_.  The following portion of the DEM parameter file shows the function called:

.. code-block:: text

    subsection mesh
      set type                                = gmsh
      set file name                           = ./mesh/spouted_structured.msh
      set expand particle-wall contact search = false
    end

where the file name includes the path to the mesh file.

Simulation Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another subsection, which is generally the one we put at the top of the parameter files, is the ``simulation control`` . ``time step``, ``end time``, ``log frequency``, and ``output frequency`` are defined here. Additionally, users can specify the output folder for the simulation results in this subsection. The ``log frequency`` parameter controls the frequency at which the iteration number is printed on the terminal. If ``log frequency = 1000`` the iteration number will be printed out every 1000 iterations. This is an easy way to monitor the progress of the simulation. A simulation time of 1 s was chosen with a time step of 0.000005. It is important to choose a long enough time as to allow all particles to come to rest. We store the output files generated in the folder ``output_dem``:


.. code-block:: text

    subsection simulation control
      set time step        = 0.00001
      set time end         = 0.8
      set log frequency    = 1000
      set output frequency = 10000
      set output path      = ./output_dem/
    end

Restart
~~~~~~~~~~~~~~~~~~~

The ``lethe-fluid-particles`` solver requires reading several DEM files to start the simulation. For this, we have to write the DEM simulation information. This is done by enabling the check-pointing option in the restart subsection. We give the written files a prefix "dem" set in the "set filename" option. The DEM parameter file is initialized exactly as the cylindrical packed bed example. The difference is in the number of particles, their physical properties, and the insertion box defined based on the new geometry. For more explanation about the individual subsections, refer to the `DEM parameters <../../../parameters/dem/dem.html>`_ and the `CFD-DEM parameters <../../../parameters/unresolved-cfd-dem/unresolved-cfd-dem.html>`_ .

.. code-block:: text

    subsection restart
      set checkpoint = true
      set frequency  = 10000
      set restart    = false
      set filename   = dem
    end



Model Parameters
~~~~~~~~~~~~~~~~~

The section on model parameters is explained in the DEM examples. We show the chosen parameters for this section:

.. code-block:: text

    subsection model parameters
      subsection contact detection
        set contact detection method                = dynamic
        set dynamic contact search size coefficient = 0.9
        set neighborhood threshold                  = 1.3
      end
      subsection load balancing
        set load balance method     = dynamic
        set threshold               = 0.5
        set dynamic check frequency = 10000
      end
      set particle particle contact force method = hertz_mindlin_limit_overlap
      set particle wall contact force method     = nonlinear
      set integration method                     = velocity_verlet
    end

We enable dynamic load balancing in order to fully take advantage of the parallelization of the code.


Lagrangian Physical Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The physical properties section of the particles allows us to specify the different parameters related to the particle such as its density, diameter, and the different coefficients that dictates the collision behavior of the particles. Also, in this section we define the total number of particles for the simulation. The gravitational acceleration as well as the physical properties of particles and walls are specified in the ``Lagrangian physical properties`` subsection. These properties include diameter and density of particles, Young's modulus, Poisson's ratio, restitution coefficient, friction and rolling friction coefficients. We insert 31,050 particles with a 2.5 mm diameter in the simulation.

.. code-block:: text

    subsection lagrangian physical properties
      set g                        = 0.0, -9.81, 0.0
      set number of particle types = 1
      subsection particle type 0
        set size distribution type            = uniform
        set diameter                          = 0.0025
        set number                            = 31050
        set density particles                 = 2526
        set young modulus particles           = 1e6
        set poisson ratio particles           = 0.25
        set restitution coefficient particles = 0.97
        set friction coefficient particles    = 0.4
        set rolling friction particles        = 0.3
      end
      set young modulus wall           = 1e6
      set poisson ratio wall           = 0.25
      set restitution coefficient wall = 0.33
      set friction coefficient wall    = 0.2
      set rolling friction wall        = 0.3
    end

Insertion Info
~~~~~~~~~~~~~~~~~~~

The ``insertion info`` subsection manages the insertion of particles. It allows us to control the insertion of particles at each time step. This section is already explained in the DEM examples. However, further information regarding the information box will be given. The volume of the insertion box should be large enough to fit all particles. Also, its bounds should be located within the mesh generated in the Mesh subsection.

.. code-block:: text

    subsection insertion info
      set insertion method                               = volume
      set inserted number of particles at each time step = 31050
      set insertion frequency                            = 2000
      set insertion box points coordinates               = -0.075, 0.0, 0 : 0.075, 0.3, 0.015
      set insertion distance threshold                   = 1.05
      set insertion maximum offset                       = 0.3
      set insertion prn seed                             = 19
    end


Floating Walls
~~~~~~~~~~~~~~~~~~~

We need to pack the particles in the bottom of the rectangular bed while preventing them from going down inside the inlet channel. Therefore, we create a stopper (floating wall) at the top of the channel. We chose the point with a y-coordinate of 0 to create the wall. We then define a normal to the wall at this point. Make sure that the end time of the floating wall is bigger than the simulation time to ensure that the particles remain outside the channel during the entire simulation time. This is shown in:

.. code-block:: text

    subsection floating walls
      set number of floating walls = 1
      subsection wall 0
        subsection point on wall
          set x = 0
          set y = 0
          set z = 0
        end
        subsection normal vector
          set nx = 0
          set ny = 1
          set nz = 0
        end
        set start time = 0
        set end time   = 50
      end
    end


---------------------------
Running the DEM Simulation
---------------------------
Launching the simulation is as simple as specifying the executable name and the parameter file. Assuming that the ``lethe-particles`` executable is within your path, the simulation can be launched in parallel as follows:

.. code-block:: text
  :class: copy-button

  mpirun -np 8 lethe-particles packing-particles.prm

.. note::
    Running the packing should take approximately 10-15 minutes on 8 cores.

After the particles have been packed inside the square bed, it is now possible to simulate the fluidization of particles.


-----------------------
CFD-DEM Parameter File
-----------------------

The CFD simulation is to be carried out using the packed bed simulated in the previous step. We will discuss the different parameter file sections. The mesh section is identical to that of the DEM so it will not be shown here.

Simulation Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simulation is run for 5 s with a time step of 0.0001 s. The time scheme chosen for the simulation is first order backward difference method (BDF1). The simulation control section is shown:

.. code-block:: text

    subsection simulation control
      set method               = bdf1
      set output frequency     = 50
      set time end             = 5
      set time step            = 0.0001
      set subdivision          = 1
      set log precision        = 10
      set output path          = ./output/
    end

Physical Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The physical properties subsection allows us to determine the density and viscosity of the fluid. We choose a density of 1 and a viscosity of 0.0000181 as to simulate the flow of air.

.. code-block:: text

    subsection physical properties
      subsection fluid 0
        set kinematic viscosity = 0.0000181
        set density             = 1
      end
    end


Initial Conditions
~~~~~~~~~~~~~~~~~~

For the initial conditions, we choose zero initial conditions for the velocity.

.. code-block:: text

    subsection initial conditions
      subsection uvwp
          set Function expression = 0; 0; 0; 0
      end
    end


Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the boundary conditions, we choose a slip boundary condition on all the walls of the bed and the channel except the inlet at the bottom of the channel and the bottom of the bed and the outlet on the top of the bed where an outlet boundary conditions was imposed.  At the base of the channel and bottom walls of the bed, we impose a Dirichlet boundary condition with an inlet velocity of 0.2 m/s and a background velocity of 1.25 respectively. For more information about the boundary conditions, please refer to the `Boundary Conditions Section <../../../parameters/cfd/boundary_conditions_cfd.html>`_

.. code-block:: text

    subsection boundary conditions
      set time dependent = false
      set number         = 4
      subsection bc 0
        set id   = 0
        set type = slip
      end
      subsection bc 1
        set id   = 2
        set type = outlet
      end
      subsection bc 2
        set id   = 1
        set type = function
        subsection u
          set Function expression = 0
        end
        subsection v
          set Function expression = 20
        end
        subsection w
          set Function expression = 0
        end
      end
      subsection bc 3
        set id   = 3
        set type = function
        subsection u
          set Function expression = 0
        end
        subsection v
          set Function expression = 1.25
        end
        subsection w
          set Function expression = 0
        end
      end
    end

The additional sections for the CFD-DEM simulations are the void fraction subsection and the CFD-DEM subsection. These subsections are described in detail in the `CFD-DEM parameters <../../../parameters/unresolved-cfd-dem/unresolved-cfd-dem.html>`_ .

Void Fraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since we are calculating the void fraction using the packed bed of the DEM simulation, we set the ``mode`` to ``dem``. For this, we need to read the dem files which we already wrote using check-pointing. We, therefore, set the ``read dem`` to ``true`` and specify the prefix of the dem files to be dem. We choose to use the quadrature centered method (QCM) to calculate the void fraction. This method does not require smoothing the void fraction as it is space and time continuous. For this simulation, we use a reference sphere having the same volume as the mesh elements as the averaging volume to calculate the void fraction.
For this, we specify the ``mode`` to be ``qcm``. We want the volume of the volume averaging sphere to be equal to the volume of the element. For this, we set the ``qcm sphere equal cell volume`` equals to ``true``.

.. code-block:: text

    subsection void fraction
      set mode                         = qcm
      set qcm sphere equal cell volume = true
      set read dem                     = true
      set dem file name                = dem
    end

CFD-DEM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We also enable grad-div stabilization in order to improve local mass conservation. The void fraction time derivative is enabled to account for the time variation of the void fraction.

.. note::
    For certain simulations, this parameter should be disabled to improve stability of the solver.

.. code-block:: text

    subsection cfd-dem
      set grad div                      = true
      set void fraction time derivative = true
      set drag force                    = true
      set buoyancy force                = true
      set shear force                   = true
      set pressure force                = true
      set saffman lift force            = false
      set drag model                    = rong
      set coupling frequency            = 100
      set implicit stabilization        = false
      set grad-div length scale         = 0.005
      set vans model                    = modelA
    end

We determine the drag model to be used for the calculation of particle-fluid forces. We enable buoyancy, drag, shear and pressure forces. For drag, we use the Rong model to determine the momentum transfer exchange coefficient. The VANS model we are solving is model A. Other possible option is model B.

Finally, the linear and non-linear solver controls are defined.

Non-linear Solver
~~~~~~~~~~~~~~~~~

.. code-block:: text

    subsection non-linear solver
      subsection fluid dynamics
      	set solver           = inexact_newton
      	set tolerance        = 1e-8
      	set max iterations   = 20
      	set verbosity        = verbose
      	set matrix tolerance = 0.75
      end
    end

We use the inexact_newton solver as to avoid the reconstruction of the system matrix at each Newton iteration. For more information about the non-linear solver, please refere to the `Non Linear Solver Section <../../../parameters/cfd/non-linear_solver_control.html>`_

Linear Solver
~~~~~~~~~~~~~

.. code-block:: text

    subsection linear solver
      subsection fluid dynamics
        set method                                = gmres
        set max iters                             = 1000
        set relative residual                     = 1e-3
        set minimum residual                      = 1e-10
        set preconditioner                        = ilu
        set ilu preconditioner fill               = 1
        set ilu preconditioner absolute tolerance = 1e-12
        set ilu preconditioner relative tolerance = 1
        set verbosity                             = verbose
        set max krylov vectors                    = 200
      end
    end

For more information about the linear solver, please refer to the `Linear Solver Section <../../../parameters/cfd/linear_solver_control.html>`_


------------------------------
Running the CFD-DEM Simulation
------------------------------

The simulation is run using the ``lethe-fluid-particles`` application. Assuming that the ``lethe-fluid-particles`` executable is within your path, the simulation can be launched as per the following command:

.. code-block:: text
  :class: copy-button

  lethe-fluid-particles spouted-bed.prm

--------
Results
--------

The results are shown in an animation below. We show the spouting of the particles as the gas is introduced from the channel at the base of the bed. Additionally, the void fraction profile is shown.
The bubble formation as well as the spouting strength are highly dependent on the drag model used. It would be interesting to try this case for different drag models.

.. raw:: html

    <iframe width="560" height="315" src="https://www.youtube.com/embed/KMVL2hPUbx8" frameborder="0" allowfullscreen></iframe>


