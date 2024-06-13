==================================
Small Scale Rotating Drum
==================================

This example of Lethe-DEM simulates dry granular flow behaviour in a small scale rotating drum. The discrete element method (DEM) is responsible for describing the behaviour of particles.  More information regarding the DEM parameters are given in the Lethe-DEM documentation, i.e. `DEM parameters <../../../parameters/dem/dem.html>`_.


----------------------------------
Features
----------------------------------
- Solvers: ``lethe-particles``
- Three-dimensional problem
- Rotational boundary
- Load-balancing


----------------------------
Files Used in This Example
----------------------------

Both files mentioned below are located in the example's folder (``examples/dem/3d-small-scale-rotating-drum``).

- Parameters file for particle insertion: ``packing-rotating-drum.prm``
- Parameters file for drum rotation: ``small-rotating-drum-dem.prm``


-----------------------
Description of the Case
-----------------------

This example simulates a rolling regime in a small scale rotating drum. First, we use Lethe-DEM to fill the bed with 20000 particles. We enable check-pointing in order to write the DEM checkpoint files for the packing which then will be used as the starting point of the DEM simulation of the rotating drum. The solver ``lethe-particles`` is used to simulate the behaviour of dry granular flow within the rotating drum.


--------------
Parameter File
--------------

Mesh
~~~~~

In this example, we choose a ``cylinder`` grid type to create a cylinder. Grid arguments are the radius of the cylinder (0.056 m) and half-length (0.051 m), respectively.  The grid is refined 3 times using the ``set initial refinement`` parameters. The ``expand particle-wall contact search`` is used in concave geometries to enable extended particle-wall contact search with boundary faces of neighbor cells for particles located in each boundary cell (for more details see `Rotating Drum example <../rotating-drum/rotating-drum.html>`_).

.. code-block:: text

    subsection mesh
      set type               = dealii
      set grid type          = cylinder
      set grid arguments     = 0.056:0.051
      set initial refinement = 3
    end


Packing information
~~~~~~~~~~~~~~~~~~~~

An insertion box is defined inside the cylindrical domain, inserting 8000 particles every 0.5 seconds while the cylinder is at rest. It is important to note the size of the insertion box to make sure it is completely inside our geometry. Otherwise, particles will be lost during the insertion stage.

.. code-block:: text


    subsection insertion info
      set insertion method                               = volume
      set inserted number of particles at each time step = 8000
      set insertion frequency                            = 100000
      set insertion box points coordinates               = -0.05, 0., -0.04 : 0.05, 0.04, 0.04
      set insertion distance threshold                   = 1.1
      set insertion maximum offset                       = 0.05
      set insertion prn seed                             = 19
    end

Restart files are written once the packing ends. The restart files are used to start the DEM simulation with the imposed rotating boundary condition.

Lagrangian Physical Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The particles are mono-dispersed with a radius of 0.0015 m and a density of 2500 kg/m3, respectively. All other particles' physical parameters are taken arbitrary and should be changed based on the physical properties and the experimental values.

.. code-block:: text

    subsection lagrangian physical properties
        set g                        = 0.0, -9.81, 0.0
        set number of particle types = 1
            subsection particle type 0
                set size distribution type            = uniform
                set diameter                          = 0.003
                set number of particles               = 20000
                set density particles                 = 2500
                set young modulus particles           = 100000000
                set poisson ratio particles           = 0.24
                set restitution coefficient particles = 0.97
                set friction coefficient particles    = 0.3
                set rolling friction particles        = 0.1

        end
        set young modulus wall           = 100000000
        set poisson ratio wall           = 0.24
        set restitution coefficient wall = 0.85
        set friction coefficient wall    = 0.35
        set rolling friction wall        = 0.1
    end


Model Parameters
~~~~~~~~~~~~~~~~~

In this example, we use the ``dynamic`` load balancing method. This method checks frequently if load balancing should be applied based on a user inputted frequency. Load balancing is dynamically applied if a certain condition is applied. More details regarding load balancing are explained in the `Rotating Drum example <../rotating-drum/rotating-drum.html>`_. 

.. code-block:: text

    subsection model parameters
      subsection contact detection
        set contact detection method                = dynamic
        set dynamic contact search size coefficient = 0.8
        set neighborhood threshold                  = 1.3
      end
      subsection load balancing
        set load balance method                     = dynamic
        set threshold                               = 0.5
        set dynamic check frequency                 = 10000
      end
      set particle particle contact force method    = hertz_mindlin_limit_overlap
      set particle wall contact force method        = nonlinear
      set rolling resistance torque method          = constant_resistance
      set integration method                        = velocity_verlet
    end

DEM Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~

The rotation of the cylinder is applied using a rotational boundary condition with a value of 1 rad/s over the x axis. Based on `deal.II boundary colouring <https://www.dealii.org/current/doxygen/deal.II/namespaceGridGenerator.html>`_, the hull of the cylinder (rotating drum) has an id = 0.

.. code-block:: text

    subsection DEM boundary conditions
      set number of boundary conditions = 1
      subsection boundary condition 0
        set boundary id         = 0
        set type                = rotational
        set rotational speed    = 1
        set rotational vector   = 1, 0, 0
      end
    end


Simulation Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The packing ``lethe-particles`` simulation was run for 2 seconds in real time.

.. code-block:: text

    subsection simulation control
      set time step        = 5e-6
      set time end         = 2
      set log frequency    = 2000
      set output frequency = 2000
      set output path      = ./output_dem/
    end
    
The actual rotation of the drum is 3 seconds in real time. We set the time equal to 5 seconds as the simulation is restarted after the packing ``lethe-particles`` simulation.

.. code-block:: text

    subsection simulation control
      set time step        = 5e-6
      set time end         = 5
      set log frequency    = 2000
      set output frequency = 2000
      set output path      = ./output_dem/
    end


-----------------------
Running the Simulation
-----------------------

The simulation is launched in two steps: the first step packs the particle in the cylinder, while the second step rotates the drum and simulates the movement of the particles. 

.. code-block:: text
  :class: copy-button

   mpirun -np 8 lethe-particles packing-rotating-drum.prm;
   mpirun -np 8 lethe-particles small-rotating-drum-dem.prm


.. note::
 This example needs a simulation time of approximately 60 minutes on 8 processors using an 12th Gen Intel(R) Core(TM) i9-12900K


---------
Results
---------

The following movie displays the rolling regime inside the rotating drum obtained with a rotational velocity of 1 rad/s. 

.. raw:: html

    <iframe width="560" height="315" src="https://www.youtube.com/embed/F-uo2lzhObk" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>


----------------------------
Possibilities for Extension
----------------------------

- Use two types of particles with different radius to prove the Brazil-Nut effect.
- Perform an unresolved CFD-DEM simulation for wet granular flows to see the impact of the hydrodynamics of the fluid over the particles dynamics.


 
