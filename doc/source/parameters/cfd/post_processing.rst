Post-Processing
---------------------
This subsection controls the post-processing other than the forces and torque on the boundary conditions. Default values are

.. code-block:: text

  subsection post-processing
	
	set calculation frequency       = 1
	set output frequency            = 1
	set verbosity                   = quiet

	#---------------------------------------------------
  	# Fluid dynamic post-processing
  	#---------------------------------------------------
	# Kinetic energy calculation
	set calculate kinetic energy    = false
	set kinetic energy name         = kinetic_energy

	# Average velocities calculation	
	set calculate average velocities  = false
	set initial time                  = 0.0

	# Pressure drop calculation
	set calculate pressure drop     = false
	set pressure drop name          = pressure_drop
	set inlet boundary id           = 0
	set outlet boundary id          = 1

	# Enstrophy calculation
	set calculate enstrophy         = false
	set enstrophy name              = enstrophy

	#---------------------------------------------------
  	# Physical properties post-processing
  	#---------------------------------------------------
	set calculate apparent viscosity  = false
	set apparent viscosity name       = apparent_viscosity

	#---------------------------------------------------
  	# Multiphysics post-processing
  	#---------------------------------------------------
	set calculate tracer statistics   = false
	set tracer statistics name        = tracer_statistics
	set calculate temperature range   = false

  end
 

* ``calculation frequency``: frequency at which the enabled post-processing is calculated.

* ``output frequency``: frequency at which the enabled post-processing is outputted in the respective file and in the terminal (if ``set verbosity = true``). For ``output frequency = 1`` (default value), results are outputted at each iteration.

* ``verbosity``: enables the display of the post-processing values in the terminal. This does not affect the printing of output files. Choices are: ``quiet`` (default, no output) or ``verbose``. Results are outputted at every ``output frequency``.

* ``calculate kinetic energy``: controls if calculation of total kinetic energy is enabled. 
	* ``kinetic energy name``: output filename for kinetic energy calculations.

* ``calculate average velocities``: controls if calculation of time-averaged velocities is enabled.
	* ``initial time``: initial time used for the average velocities calculations.

* ``calculate pressure drop``: controls if calculation of the pressure drop from the inlet boundary to the outlet boundary is enabled. 
	* ``inlet boundary id`` and ``outlet boundary id``: defines the id for inlet and outlet boundaries, respectively. 
	* ``pressure drop name``: output filename for pressure drop calculations.

* ``calculate enstrophy``: controls if calculation of total enstrophy, which corresponds to dissipation effects in the fluid, is enabled. 
	* ``enstrophy name``: output filename for enstrophy calculations.

* ``calculate apparent viscosity``: controls if parameter calculation of an apparent viscosity is enabled, when using a non Newtonian flow (see section `Rheological models <https://lethe-cfd.github.io/lethe/parameters/cfd/physical_properties.html#rheological-models>`_). This is mainly used to define the Reynolds number `a posteriori`. 
	* ``apparent viscosity name``: output filename for apparent viscosity calculations.

* ``calculate tracer statistics``: controls if calculation of tracer statistics is enabled. Statistics include: minimum, maximum, average and standard-deviation.
	* ``tracer statistics name``: output filename for tracer statistics calculations.

.. warning::

	Do not forget to ``set tracer = true`` in the :doc:`multiphysics` subsection of the ``.prm``.

* ``calculate heat transfer statistics``: controls if calculation of tracer statistics is enabled. Statistics include: minimum, maximum, average and standard-deviation on the temperature.
	* ``heat transfer statistics name``: output filename for tracer statistics calculations.

.. warning::

	Do not forget to ``set heat transfer = true`` in the :doc:`multiphysics` subsection of the ``.prm``.
