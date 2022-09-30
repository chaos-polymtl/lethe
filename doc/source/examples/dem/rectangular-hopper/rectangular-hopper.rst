==================================
Rectangular hopper
==================================

This example simulates the filling and discharge of particles in a rectangular hopper.
We set up this simulation based on the simulation of Anand et al. `[1] <https://doi.org/10.1016/j.ces.2008.08.015>`_. It is recommended to visit `DEM parameters <../../../parameters/dem/dem.html>`_ for more detailed information on the concepts and physical meanings of the parameters in Lethe-DEM.

Features
----------------------------------
- Solvers: ``dem_3d``
- Floating walls
- Gmsh grids
- Python post-processing script using PyVista


Files used in this example
----------------------------
``/examples/dem/rectangular-hopper/hopper.prm``


Description of the case
-----------------------

This simulation consists of two stages: filling (0-4 s) and discharge (4-7.5 s) of particles. Anand et al. uses periodic boundaries in the z axis allowing to use a thin width for simulation. This case do not use periodic boundary in Lethe. Prior minimizing the impact of collision of particle with walls in z axis, the width and the number of particle were multiplied 6 times. This corresponds to a width of 15 times the particle diameter.

Parameter file
--------------

Mesh
~~~~~

The mesh is a hopper with 90° angle generated with GMSH having a top part for the filling and a bottom part which acts as a collector of the particle.

.. code-block:: text

    subsection mesh
        set type                                = gmsh
        set file name                           = hopper90.msh
        set initial refinement                  = 1
        set expand particle-wall contact search = false
        set check diamond cells                 = true
    end


Insertion info
~~~~~~~~~~~~~~~~~~~

An insertion box is defined inside and the top part of the hopper. The inserted number is choose to be a factor of the total number of particle. In this case, 14 insertion step is require to fill up the hopper with particles.

.. code-block:: text

    subsection insertion info
        set insertion method                               = non_uniform
        set inserted number of particles at each time step = 2910
        set insertion frequency                            = 25000
        set insertion box minimum x                        = -0.1030
        set insertion box minimum y                        = 0.10644
        set insertion box minimum z                        = .00224
        set insertion box maximum x                        = 0.1030
        set insertion box maximum y                        = 0.16020
        set insertion box maximum z                        = 0.03136
        set insertion distance threshold                   = 1.5
        set insertion random number range                  = 0.1
        set insertion random number seed                   = 20
    end


Lagrangian physical properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The total number of particles in this simulation is 40740 with a dismeter of 2.24 mm.

.. code-block:: text

    subsection lagrangian physical properties
        set gx                       = 0.0
        set gy                       = -9.81
        set gz                       = 0.0
        set number of particle types = 1
        subsection particle type 0
            set size distribution type            = uniform
            set diameter                          = 0.00224
            set number                            = 40740
            set density particles                 = 2500
            set young modulus particles           = 1000000
            set poisson ratio particles           = 0.3
            set restitution coefficient particles = 0.94
            set friction coefficient particles    = 0.2
            set rolling friction particles        = 0.09
        end
        set young modulus wall           = 1000000
        set poisson ratio wall           = 0.3
        set friction coefficient wall    = 0.2
        set restitution coefficient wall = 0.9
        set rolling friction wall        = 0.09
    end


Model parameters
~~~~~~~~~~~~~~~~~

.. code-block:: text

    subsection model parameters
        set contact detection method                = dynamic
        set dynamic contact search size coefficient = 0.9
        set load balance method                     = frequent
        set load balance frequency                  = 200000
        set neighborhood threshold                  = 1.3
        set particle particle contact force method  = hertz_mindlin_limit_overlap
        set particle wall contact force method      = nonlinear
        set integration method                      = velocity_verlet
    end


Simulation control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

    subsection simulation control
        set time step        = 1e-5
        set time end         = 7.5
        set log frequency    = 1000
        set output frequency = 1000
        set output path      = ./output/
        set output name      = hopper
    end


Floating walls
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Floating wall in this example is handled as explained in the `Silo example <../silo/silo.html>`_.

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
    set end time   = 4
    end
end


Running the simulation
----------------------
This simulation can be launched by (in parallel mode on 32 processes):

.. code-block:: text

  mpirun -np 8 dem_3d hopper.prm

Results post-processing
---------


Results
---------





Reference
---------
`[1] <https://doi.org/10.1016/j.ces.2008.08.015>`_ Anand, A., Curtis, J. S., Wassgren, C. R., Hancock, B. C., & Ketterhagen, W. R. (2008). Predicting discharge dynamics from a rectangular hopper using the discrete element method (DEM). Chemical Engineering Science, 63(24), 5821-5830.
