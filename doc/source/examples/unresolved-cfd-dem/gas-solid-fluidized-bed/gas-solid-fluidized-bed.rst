==================================
Gas-Solid Fluidized Bed
==================================

It is strongly recommended to visit `DEM parameters <../../../parameters/dem/dem.html>`_  and `CFD-DEM parameters <../../../parameters/unresolved-cfd-dem/unresolved-cfd-dem.html>`_ for more detailed information on the concepts and physical meaning of the parameters ind DEM and CFD-DEM.


----------------------------------
Features
----------------------------------

- Solvers: ``lethe-particles`` and ``lethe-fluid-particles``
- Three-dimensional problem
- Displays the selection of models and physical properties.
- Simulates a solid-gas fluidized bed.


---------------------------
Files Used in This Example
---------------------------

Both files mentioned below are located in the example's folder (``examples/unresolved-cfd-dem/gas-solid-fluidized-bed``).

- Parameter file for particle generation and packing: ``packing-particles.prm``
- Parameter file for CFD-DEM simulation of the gas-solid fluidized bed: ``gas-solid-fluidized-bed.prm``


-----------------------
Description of the Case
-----------------------

This example simulates the fluidization of spherical particles in air. First, we use ``lethe-particles`` to fill the bed with particles. We enable check-pointing in order to write the DEM checkpoint files which will be used as the starting point of the CFD-DEM simulation. Then, we use the ``lethe-fluid-particles`` solver within Lethe to simulate the fluidization of the particles by initially reading the checkpoint files from the DEM simulation.


-------------------
DEM Parameter File
-------------------

All parameter subsections are described in the `parameter section <../../../parameters/parameters.html>`_ of the documentation.

To set-up the square fluidized bed case, we first fill the bed with particles. 

We first introduce the different sections of the parameter file ``packing-particles.prm`` needed to run this simulation.

Mesh
~~~~~

In this example, we are simulating a squared fluidized bed that has a half length of 0.1 m, and a side of 0.04 m. We use the `subdivided_hyper_rectangle GridGenerator <https://www.dealii.org/current/doxygen/deal.II/namespaceGridGenerator.html#ac76417d7404b75cf53c732f456e6e971>`_  in order to generate the mesh. The square bed is divided 40 times in the y direction. The following portion of the DEM parameter file shows the function called:

.. code-block:: text

    subsection mesh
      set type               = dealii
      set grid type          = subdivided_hyper_rectangle
      set grid arguments     = 1,5,1:-0.02,-0.1,-0.02:0.02,0.1,0.02:true
      set initial refinement = 3
    end
    
Simulation Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another subsection, which is generally the one we put at the top of the parameter files, is the ``simulation control`` . ``time step``, end time, log, and ``output frequency`` are defined here. Additionally, users can specify the output folder for the simulation results in this subsection. The ``log frequency`` parameter controls the frequency at which the iteration number is printed on the terminal. If ``log frequency = 1000`` the iteration number will be printed out every 1000 iterations. This is an easy way to monitor the progress of the simulation. A simulation time of 1 s was chosen with a time step of 0.000005. It is important to choose a long enough time as to allow all particles to come to rest. We store the output files generated in the folder output_dem:


.. code-block:: text

    subsection simulation control
      set time step        = 0.000005
      set time end         = 0.5
      set log frequency    = 1000
      set output frequency = 1000
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
        set contact detection method = dynamic
        set neighborhood threshold   = 1.3
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

The physical properties section of the particles allows us to specify the different parameters related to the particle such as its density, diameter, and the different coefficients that dictates the collision behavior of the particles. Also, in this section we define the total number of particles for the simulation. The gravitational acceleration as well as the physical properties of particles and walls are specified in the ``Lagrangian physical properties`` subsection. These properties include diameter and density of particles, Young's modulus, Poisson's ratio, restitution coefficient, friction and rolling friction coefficients. We insert 30,000 particles in the simulation.

.. code-block:: text

    subsection lagrangian physical properties
      set g                        = 0.0, -9.81, 0.0
      set number of particle types = 1
      subsection particle type 0
        set size distribution type            = uniform
        set diameter                          = 0.001
        set number                            = 30000
        set density particles                 = 1500
        set young modulus particles           = 1e6
        set poisson ratio particles           = 0.3
        set restitution coefficient particles = 0.2
        set friction coefficient particles    = 0.1
        set rolling friction particles        = 0.2
      end
      set young modulus wall           = 1e6
      set poisson ratio wall           = 0.3
      set restitution coefficient wall = 0.2
      set friction coefficient wall    = 0.1
      set rolling friction wall        = 0.3
    end
    
Insertion Info
~~~~~~~~~~~~~~~~~~~

The ``insertion info`` subsection manages the insertion of particles. It allows us to control the insertion of particles at each time step. This section is already explained in the DEM examples. However, further information regarding the information box will be given. The volume of the insertion box should be large enough to fit all particles. Also, its bounds should be located within the mesh generated in the Mesh subsection.  

.. code-block:: text

  subsection insertion info
      set insertion method                               = volume
      set inserted number of particles at each time step = 30000
      set insertion frequency                            = 100000
      set insertion box points coordinates               = -0.018, -0.05, -0.018 : 0.018, 0.05, 0.018
      set insertion distance threshold                   = 1.5
      set insertion maximum offset                       = 0.2
      set insertion prn seed                             = 19
    end


Floating Walls
~~~~~~~~~~~~~~~~~~~

We need to pack the particles in the middle of the square bed. Therefore, we create a stopper (floating wall) somewhere below the center of the bed. We chose the point with a y-coordinate of -0.06 to create the wall. We then define a normal to the wall at this point. Make sure that the end time of the floating wall is bigger than the simulation time to ensure that the particles remain suspended. This is shown in:

.. code-block:: text

    subsection floating walls
      set number of floating walls = 1
      subsection wall 0
        subsection point on wall
          set x = 0
          set y = -0.06
          set z = 0
        end
        subsection normal vector
          set nx = 0
          set ny = 1
          set nz = 0
        end
        set start time = 0
        set end time   = 5
      end
    end


---------------------------
Running the DEM Simulation
---------------------------
Launching the simulation is as simple as specifying the executable name and the parameter file. Assuming that the ``lethe-particles`` executable is within your path, the simulation can be launched on a single processor by typing:

.. code-block:: text
  :class: copy-button

  lethe-particles packing-particles.prm

or in parallel (where 8 represents the number of processors)

.. code-block:: text
  :class: copy-button

  mpirun -np 8 lethe-particles packing-particles.prm

Lethe will generate a number of files. The most important one bears the extension ``.pvd``. It can be read by popular visualization programs such as `Paraview <https://www.paraview.org/>`_. 


.. note:: 
    Running the packing should take approximately 20 minutes on 8 cores.

After the particles have been packed inside the square bed, it is now possible to simulate the fluidization of particles.


-----------------------
CFD-DEM Parameter File
-----------------------

The CFD simulation is to be carried out using the packed bed simulated in the previous step. We will discuss the different parameter file sections. The mesh section is identical to that of the DEM so it will not be shown here.

Simulation Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simulation is run for 1 s with a time step of 0.002 s. The time scheme chosen for the simulation is first order backward difference method (BDF1). The simulation control section is shown:

.. code-block:: text

    subsection simulation control
      set method               = bdf1
      set output frequency     = 10
      set time end             = 1
      set time step            = 0.001
      set output path          = ./output/
    end

Physical Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The physical properties subsection allows us to determine the density and viscosity of the fluid. We choose a density of 1 and viscosity of 0.00001 as to simulate the flow of air. 

.. code-block:: text

    subsection physical properties
      subsection fluid 0
        set kinematic viscosity = 0.00001
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

For the boundary conditions, we choose a slip boundary condition on the walls of the square bed (IDs = 0, 1, 4, 5) and an inlet velocity of 0.2 m/s at the lower face of the bed (ID = 2).

.. code-block:: text

  subsection boundary conditions
    set number = 6
    subsection bc 0
      set id   = 0
      set type = slip
    end
    subsection bc 1
      set id   = 1
      set type = slip
    end
    subsection bc 2
      set id   = 2
      set type = function
      subsection u
        set Function expression = 0
      end
      subsection v
        set Function expression = 2
      end
      subsection w
        set Function expression = 0
      end
    end
    subsection bc 3
      set id   = 3
      set type = outlet
    end
    subsection bc 4
      set id   = 4
      set type = slip
    end
    subsection bc 5
      set id   = 5
      set type = slip
    end
  end

The additional sections for the CFD-DEM simulations are the void fraction subsection and the CFD-DEM subsection. These subsections are described in detail in the `CFD-DEM parameters <../../../parameters/unresolved-cfd-dem/unresolved-cfd-dem.html>`_ .

Void Fraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since we are calculating the void fraction using the packed bed of the DEM simulation, we set the mode to "dem". For this, we need to read the dem files which we already wrote using check-pointing. We, therefore, set the read dem to "true" and specify the prefix of the dem files to be dem. We now choose a smoothing factor for the void fraction to reduce discontinuity which can lead to oscillations in the velocity. The factor we choose is around the square of twice the particle's diameter. 
 
.. code-block:: text

    subsection void fraction
        set mode                = pcm
        set read dem            = true
        set dem file name       = dem
        set l2 smoothing factor = 0.000005
    end

CFD-DEM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We also enable grad_div stabilization in order to improve local mass conservation. The void fraction time derivative is enabled to account for the time variation of the void fraction. 

.. note:: 
    For certain simulations, this parameter should be disabled to improve stability of the solver.

.. code-block:: text

    subsection cfd-dem
        set grad div                      = true
        set void fraction time derivative = true
        set drag force                    = true
        set buoyancy force                = true
        set shear force                   = false
        set pressure force                = false
        set drag model                    = difelice
        set coupling frequency            = 100
        set vans model                    = modelB
    end
    
We determine the drag model to be used for the calculation of particle-fluid forces as the Di Felice model. Other optional forces that can be enabled are the buoyancy force, the shear force and the pressure force. We only decide to enable drag and buoyancy as for air, the other forces are considered to be negligible. The VANS model we are solving is model B. Other possible option is model A.

Finally, the linear and non-linear solver controls are defined.

Non-linear Solver
~~~~~~~~~~~~~~~~~

We use the inexact Newton non-linear solver to minimize the number of time the matrix of the system is assembled. This is used to increase the speed of the simulation, since the matrix assembly requires significant computations.

.. code-block:: text

  subsection non-linear solver
    subsection fluid dynamics
      set solver           = inexact_newton
      set tolerance        = 1e-7
      set max iterations   = 20
      set matrix tolerance = 0.2
      set verbosity        = verbose
    end
  end
    
Linear Solver
~~~~~~~~~~~~~

.. code-block:: text

    subsection linear solver
      subsection fluid dynamics
        set method                                = gmres
        set max iters                             = 5000
        set relative residual                     = 1e-3
        set minimum residual                      = 1e-11
        set preconditioner                        = ilu
        set ilu preconditioner fill               = 1
        set ilu preconditioner absolute tolerance = 1e-14
        set ilu preconditioner relative tolerance = 1.00
        set verbosity                             = verbose
        set max krylov vectors                    = 200
      end
    end


------------------------------
Running the CFD-DEM Simulation
------------------------------

The simulation is run using the ``lethe-fluid-particles`` application. Assuming that the ``lethe-fluid-particles`` executable is within your path, the simulation can be launched as per the following command:

.. code-block:: text
  :class: copy-button

  lethe-fluid-particles fluidized-bed.prm

--------
Results
--------

The results are shown in an animation below. We show the fluidization of the particles as the gas is introduced from the bottom of the bed.

.. raw:: html

    <iframe width="560" height="315" src="https://www.youtube.com/embed/ygJI42x4K5c" frameborder="0" allowfullscreen></iframe>
    

    
