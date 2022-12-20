==================================
Rotating drum
==================================

This is the third example of Lethe-DEM. This example simulates a rotating drum. We setup this simulation according to the experiments of Alizadeh et al `[1] <https://doi.org/10.1002/aic.13982>`_. It is recommended to visit `DEM parameters <../../../parameters/dem/dem.html>`_ for more detailed information on the concepts and physical meanings of the parameters in Lethe-DEM.

Features
----------------------------------
- Solvers: ``dem_3d``
- Rotational boundary
- Load-balancing


Files used in this example
----------------------------
``/examples/dem/3d-rotating-drum/rotating-drum.prm``


Description of the case
-----------------------

226080 particles are first inserted into a cylindrical domain and then start rolling on the cylinder wall because of the rotation of the cylinder. The rotation of the cylinder is applied using a ``rotational`` boundary condition. Because of the large number of particles, this simulation should be launched in parallel mode with load-balancing. The concepts and different types of ``boundary condition`` and load-balancing are explained in this example.


Parameter file
--------------

Mesh
~~~~~

In this example, we choose a ``cylinder`` grid type to create a cylinder. Grid arguments are the radius and half-length, respectively. Therefore, the specified grid arguments create a cylinder with a diameter of 0.24 m and a length of 0.36 m. The grid is refined 4 times to reach the desired cell size to particle diameter ratio (see packing in ball example for more details). The ``expand particle-wall contact search`` is used in concave geometries to enable extended particle-wall contact search with boundary faces of neighbor cells for particles located in each boundary cell. Between two consecutive contact search steps, where the containing cells of the particles are updated (i.e. particles are mapped into cells), particles located close to the cell boundaries may change cells. If this situation occurs for a particle on a boundary, the particle-wall collision is calculated using the information (normal vector and location of the vertices) of the previous boundary cell. Hence, when the containing cell of the particle updates, the particle may already have a large overlap with the new cell that leads to a large contact force/velocity of the particles. We use ``expand particle-wall contact search`` to avoid this undesired situation.

.. code-block:: text

    subsection mesh
      set type                                = dealii
      set grid type                           = cylinder
      set grid arguments                      = 0.12:0.18
      set initial refinement                  = 4
      set expand particle-wall contact search = true
    end


Insertion info
~~~~~~~~~~~~~~~~~~~

An insertion box is defined inside the cylindrical domain. 75110 particles are inserted non-uniformly at each insertion step.

.. code-block:: text

    subsection insertion info
      set insertion method                               = non_uniform
      set inserted number of particles at each time step = 37555
      set insertion frequency                            = 150000
      set insertion box minimum x                        = -0.175
      set insertion box minimum y                        = -0.07
      set insertion box minimum z                        = 0
      set insertion box maximum x                        = 0.175
      set insertion box maximum y                        = 0.07
      set insertion box maximum z                        = 0.09
      set insertion distance threshold                   = 1.575
      set insertion random number range                  = 0.05
      set insertion random number seed                   = 19
    end


Lagrangian physical properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The particles (226080 particles) are monodispersed, their diameter and density are 0.003 m and 2500 kg/m3, respectively.

.. code-block:: text

    subsection lagrangian physical properties
      set gx                       = 0.0
      set gy                       = 0.0
      set gz                       = -9.81
      set number of particle types = 1
      subsection particle type 0
        set size distribution type            = uniform
        set diameter                          = 0.003
        set number                            = 226080
        set density particles                 = 2500
        set young modulus particles           = 100000000
        set poisson ratio particles           = 0.24
        set restitution coefficient particles = 0.97
        set friction coefficient particles    = 0.3
      end
      set young modulus wall           = 100000000
      set poisson ratio wall           = 0.24
      set restitution coefficient wall = 0.85
      set friction coefficient wall    = 0.35
    end


Model parameters
~~~~~~~~~~~~~~~~~

Load-balancing updates the distribution of the subdomains between the processes in parallel simulation to achieve better computational performance (less simulation time). Three load-balancing methods are available in Lethe-DEM: ``once``, ``frequent``, or ``dynamic``. Read `this article <https://www.mdpi.com/2227-9717/10/1/79>`_ for more information about different load-balancing methods and their performances in various types of DEM simulations.

Selecting ``repartition method = once``, requires defining the step at which the code calls load balancing (``load balance step``). ``Frequent`` ``repartition method`` requires defining ``load balance frequency``, and in ``dynamic`` ``repartition method``, we should define ``load balance threshold`` and ``dynamic load balance check frequency``. In ``dynamic`` load balancing, the code checks the distribution of particles among the processors, every ``dynamic load balance check frequency`` steps, and if

.. math::
    L_{max}-L_{min}>{\beta}\bar{L}

it calls load-balancing. :math:`{L}` and :math:`{\beta}` denote computational load on a process and ``load balance threshold``, respectively.

In the rotating drum simulation, we use a ``once`` load-balancing method, since particles occupy a constant region inside the rotating drum after reaching steady-state operation.

.. code-block:: text

    subsection model parameters
      set contact detection method                = dynamic
      set dynamic contact search size coefficient = 0.8
      set neighborhood threshold                  = 1.3
      set load balance method                     = once
      set load balance step                       = 150000
      set particle particle contact force method  = hertz_mindlin_limit_overlap
      set particle wall contact force method      = nonlinear
      set integration method                      = velocity_verlet
    end


Boundary Condition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this subsection, the boundary conditions of the DEM simulation are defined. First of all, the ``number of boundary conditions`` is specified. Then for each boundary condition, its information is defined. There are four boundary types: ``fixed_wall``, ``outlet``, ``rotational`` (around the center), and ``translational``. For ``rotational`` motion, ``rotational speed`` and ``rotational vector`` are required, while for ``translational`` motion, the ``speed`` should be defined in each direction.

``fixed_wall`` is a static wall, and particles collide with these static walls upon reaching them. The only way to move these walls is to move the entire triangulation. If the ``outlet`` condition is chosen for a boundary, particles can leave the simulation domain via this outlet. Using ``rotational`` or ``translational`` boundary conditions exerts imaginary rotational and translational velocities to that boundary. In other words, the boundary does not move, but the particles that have collisions with these walls feel a rotational or translational velocity from the wall. This feature is used in the rotating drum example. The boundary id of the ``cylinder`` side wall, defined with deal.ii grid generator is 4. We set the ``rotational speed`` equal to 11.6 rad/s, and the cylinder should rotate around its axis (`x` direction).

.. code-block:: text

    subsection DEM boundary conditions
      set number of boundary conditions = 1
      subsection boundary condition 0
        set boundary id         = 4
        set type                = rotational
        set rotational speed    = 11.6
        set rotational vector x = 1
        set rotational vector y = 0
        set rotational vector z = 0
      end
    end


Simulation control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

    subsection simulation control
      set time step        = 1e-6
      set time end         = 15
      set log frequency    = 1000
      set output frequency = 1000
    end

Running the simulation
----------------------
This simulation can be launched (in parallel mode on 64 processes) by:

.. code-block:: text

  mpirun -np 64 dem_3d rotating-drum.prm


.. warning::
	This example needs a simulation time of approximately 48 hours 64 cores. This high computational cost is because of the large number of particles.


Results
---------

Animation of the rotating drum simulation:

.. raw:: html

    <iframe width="560" height="315" src="https://www.youtube.com/embed/krM_rFIDHAA" frameborder="0" allowfullscreen></iframe>


Reference
---------

`[1] <https://doi.org/10.1002/aic.13982>`_ Alizadeh, E., Dubé, O., Bertrand, F. and Chaouki, J., 2013. Characterization of mixing and size segregation in a rotating drum by a particle tracking method. AIChE Journal, 59(6), pp.1894-1905.
