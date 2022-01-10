Post-Processing
---------------------
This subsection controls the post-processing other than the forces and torque on the boundary conditions. Default values are

.. code-block:: text

  subsection post-processing
  set calculate enstrophy           = false
  set calculate kinetic energy      = false
  set calculate apparent viscosity  = false
  set calculate pressure drop       = false
  set calculate average velocities  = false
  set calculate tracer statistics   = false
  set initial time                  = 0.0
  set inlet boundary id           = 0
  set outlet boundary id          = 1
  set enstrophy name              = enstrophy
  set kinetic energy name         = kinetic_energy
  set apparent viscosity name     = apparent_viscosity
  set pressure drop name          = pressure_drop
  set tracer statistics name      = tracer_statistics
  set verbosity                   = quiet
  set calculation frequency       = 1
  set output frequency            = 1
  end
 

* The ``calculate enstrophy`` parameter enables calculation of total enstrophy.

* The ``calculate kinetic energy`` parameter enables calculation of total kinetic energy.

* The ``calculate apparent viscosity`` parameter enables calculation of an apparent viscosity when using a non Newtonian flow. This is mainly used to define  Reynolds number a posteriori. 

* The ``calculate pressure drop`` parameter enables calculation of the pressure drop from the inlet boundary to the outlet boundary.

* The ``calculate average velocities`` parameter enables calculation of time-averaged velocities.

* The ``calculate tracer statistics`` parameter enables calculation of tracer statistics.

* The ``initial time`` parameter sets the initial time to start calculations of time-averaged velocities.

* The ``inlet boundary id`` parameter sets the inlet boundary id for pressure drop calculations. 

* The ``outlet boundary id`` parameter sets the outlet boundary id for pressure drop calculations. 

* The ``enstrophy name`` parameter sets output file name for enstrophy calculations.

* The ``kinetic energy name`` parameter sets output file name for kinetic energy calculations.

* The ``apparent viscosity name`` parameter sets output file name for apparent viscosity calculations.

* The ``pressure drop name`` parameter sets output file name for pressure drop calculations.

* The ``tracer statistics name`` parameter sets output file name for tracer statistics calculations.

* The ``verbosity`` parameter specifies if the post-processing values are printed to the terminal at the end of every <output frequency> iteration. This does not affect the printing of output files.

* The ``calculation frequency`` parameter sets the frequency at which the enstrophy, kinetic energy, apparent viscosity, pressure drop, average velocities and tracer statistics are calculated (for the ones that are set to ``true``). 

* The ``output frequency`` parameter sets the frequency at which the enstrophy, kinetic energy, apparent viscosity, pressure drop and tracer statistics are output to their respective files (for the ones that are set to ``true``). 
