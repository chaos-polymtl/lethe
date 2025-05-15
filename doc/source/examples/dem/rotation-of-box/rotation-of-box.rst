==================================
Rotation of Box
==================================

This example simulates a triangulation (box) rotation. It is recommended to visit `DEM parameters <../../../parameters/dem/dem.html>`_ for more detailed information on the concepts and physical meanings of the parameters in Lethe-DEM.

----------------------------------
Features
----------------------------------
- Solvers: ``lethe-particles``
- Rotating a triangulation


----------------------------
Files Used in This Example
----------------------------

- Parameter file: ``examples/dem/3d-grid-rotation-in-box/grid-rotation-box.prm``


-----------------------
Description of the Case
-----------------------

4000 particles are inserted in a rotating box and rotate with the box. In this example, the whole triangulation is rotated.


--------------
Parameter File
--------------

Mesh
~~~~~

The ``grid type`` in this example is a ``hyper_cube``. Its dimensions are 0.04 m in every direction (from -0.02 m to 0.02 m), and it is refined 3 times.

.. code-block:: text

    subsection mesh
      set type               = dealii
      set grid type          = hyper_cube
      set grid arguments     = -0.02 : 0.02 : false
      set initial refinement = 3
    end


Insertion Info
~~~~~~~~~~~~~~~~~~~

An insertion box is defined inside the cubic domain. 4000 particles are inserted non-uniformly in the first iteration.

.. code-block:: text

    subsection insertion info
      set insertion method                               = volume
      set inserted number of particles at each time step = 4000
      set insertion frequency                            = 2000000
      set insertion box points coordinates               = -0.019, -0.019, -0.01 : 0.019, 0.019, 0.019
      set insertion distance threshold                   = 1.5
      set insertion maximum offset                       = 0.2
      set insertion prn seed                             = 19
    end


Lagrangian Physical Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``number`` of particles (4000) is equal to the specified ``inserted number of particles at each time step``. This means that all the particles are inserted during the first insertion iteration (if the inserted number of particles fits inside the specified insertion box).

.. code-block:: text

    subsection lagrangian physical properties
      set g                        = 0.0, 0.0, -9.81
      set number of particle types = 1
      subsection particle type 0
        set size distribution type            = uniform
        set diameter                          = 0.001
        set number of particles               = 4000
        set density particles                 = 1000
        set young modulus particles           = 1000000
        set poisson ratio particles           = 0.3
        set restitution coefficient particles = 0.3
        set friction coefficient particles    = 0.1
        set rolling friction particles        = 0.05
      end
      set young modulus wall           = 1000000
      set poisson ratio wall           = 0.3
      set restitution coefficient wall = 0.3
      set friction coefficient wall    = 0.1
      set rolling friction wall        = 0.05
    end


Model Parameters
~~~~~~~~~~~~~~~~~

.. code-block:: text

    subsection model parameters
      subsection contact detection
        set contact detection method                = dynamic
        set dynamic contact search size coefficient = 0.9
        set neighborhood threshold                  = 1.3
      end
      set particle particle contact force method    = hertz_mindlin_limit_overlap
      set particle wall contact force method        = nonlinear
      set rolling resistance torque method          = constant_resistance
      set integration method                        = velocity_verlet
    end


Simulation Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

    subsection simulation control
      set time step        = 1e-5
      set time end         = 5
      set log frequency    = 1000
      set output frequency = 1000
    end


----------------------
Running the Simulation
----------------------
This simulation can be launched by:

.. code-block:: text
  :class: copy-button

  lethe-particles grid-rotation-box.prm


---------
Results
---------

Animation of the rotating box simulation:

.. raw:: html

    <p align="center"><iframe width="560" height="315" src="https://www.youtube.com/embed/zGjEVskObIc" frameborder="0" allowfullscreen></iframe>
