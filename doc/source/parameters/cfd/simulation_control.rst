==================
Simulation Control
==================

This subsection is the most important in a simulation and therefore, the most commonly modified in a parameter file. It controls the general parameters of the simulation, such as, the time integration method, end time of a simulation and output settings for paraview files. 

.. tip::
	A standard convention in Lethe is to keep this section at the top of the parameter file, since it is generally the most accessed one.

.. code-block:: text

  subsection simulation control
    # Type of solver or time-stepping scheme
    set method = steady 
  
    #---------------------------------------------------
    # Steady-state simulation parameters
    #---------------------------------------------------
    # Number of mesh adaptation
    set number mesh adapt = 0
  
    # Tolerance at which the simulation is stopped
    set stop tolerance = 1e-10
  
    #---------------------------------------------------
    # BDF scheme parameters
    #---------------------------------------------------
    # Method used to startup high order BDF methods
    set bdf startup method = multiple step bdf
  
    # Scaling factor used in the iterations necessary to startup the BDF schemes
    set startup time scaling = 0.4
  
    #---------------------------------------------------
    # Transient simulations parameters
    #---------------------------------------------------
    # End time value of the simulation
    set time end = 1
  
    # Time step value
    set time step = 1
  
    # Adaptative time-stepping
    set adapt = false

    # True if the time-step should be overriden upon restart
    set override time step on restart = false
  
    # Maximum CFL value
    set max cfl = 1
  
    # Maximum time step value
    set max time step = 1e6
  
    # Adaptative time step scaling
    set adaptative time step scaling = 1.1
  
    # Time step independent of end time
    set time step independent of end time = true
  
    #---------------------------------------------------
    # Log file parameters
    #---------------------------------------------------
    # Log frequency
    set log frequency = 1
  
    # Display precision when writing to log
    set log precision = 6
  
    #---------------------------------------------------
    # Output file parameters
    #---------------------------------------------------
    # File output path
    set output path = ./
  
    # File output prefix
    set output name = out
  
    # The control for the output of the simulation results
    set output control = iteration
  
    # Output iteration frequency
    set output frequency = 1
  
    # Output time frequency
    set output time frequency = -1
  
    # List of specific output times
    set output times = -1
  
    # Output time interval
    set output time interval = 0, 1.7e308
  
    # Maximum number of vtu output files
    set group files = 1
  
    # Output the boundaries of the domain along with their ID
    set output boundaries = false
  
    # Subdivision of mesh cell in postprocessing
    set subdivision = 1
  end

* ``method``: time-stepping method used. The available options are: 
	* ``steady``: steady-state simulation
	* ``steady_bdf``: steady-state simulation using adjoint time stepping with a bdf1 scheme
	* ``bdf1``: 1st order backward differentiation
	* ``bdf2``: 2nd order backward differentiation
	* ``bdf3``: 3rd order backward differentiation
	* ``sdirk22``: 2nd order 2 stages singly diagonally implicit Runge-Kutta
	* ``sdirk33``: 3rd order 3 stages singly diagonally implicit Runge-Kutta

.. warning::
	For now, the SDIRK schemes are not supported by any physics other than Fluid Dynamics.

-----------------------------------
Steady-state simulation parameters
-----------------------------------

* ``number mesh adapt``: number of mesh adaptations during the steady-state simulation.

* ``stop tolerance``: tolerance at which the adjoint time stepping steady-state simulation (``method = steady_bdf``) stops. The adjoint time stepping will stop when the :math:`\mathcal{L}_2` norm of the initial residual is lower than ``stop tolerance`` at the beginning of a time step.

----------------------
BDF scheme parameters
----------------------

* ``bdf startup method``: scheme used to start a high order BDF scheme (2nd order and above). The available options are: 
	* ``multiple step bdf``:  A lower order BDF scheme is used to start the simulation. For example, in the case of ``bdf3``, the first step is done using ``bdf1``, the second with ``bdf2`` and the third and onward are done with ``bdf3``.
	* ``initial solution``: In this case, a time-dependent initial solution is provided and that initial solution is used to start the time stepping. This is mostly useful when using the method of manufactured solutions to establish the formal accuracy of the BDF time stepping schemes. 

* ``startup time scaling``: scaling factor used in the iterations necessary to startup the BDF schemes.

.. note::
	SDIRK schemes are self-starting and do not require any additional parameter.

---------------------------------
Transient simulations parameters
---------------------------------

* ``time end``: value of the time at which the simulation ends.

* ``time step``: value of the time step.

* ``adapt``: controls if adaptive time-stepping is enabled. If set to ``true``, the time-step will evolve to ensure that the ``max cfl`` value is reached.

* ``override time step on restart``: controls if the time step should be overridden by the set value upon restart. If set to ``true``, the time-step will be set to the value of ``time step`` and the time-step value recorded at the last checkpoint will be overridden at the start of the simulation.

* ``max cfl``: maximum value of the :math:`\text{CFL}` condition number that can be reached during the simulation. This parameter is only used when ``set adapt = true``. The :math:`N_{\mathrm{CFL}}` is calculated as:

  .. math::
    N_{\mathrm{CFL}} = \max_q \frac{|\mathbf{u}_q| \Delta t} {h}

  where :math:`q` the Gauss points and :math:`|\mathbf{u}_q|` is the velocity at the Gauss points. Essentially, the maximum CFL is the max of the CFL evaluated at every Gauss point in the mesh.

* ``max time step``: maximum time step value that can be reached during the simulation. This parameter is only used when ``set adapt = true``. It is useful when the problem of interest has an additional time step constraint such as the capillary time step limit described in :doc:`../../examples/multiphysics/capillary-wave/capillary-wave`.

* ``adaptative time step scaling``: rate of increase of the time step value. The new time step value is fixed by ``adaptative time step scaling`` * ``previous value of the time step``.

* ``time step independent of end time``: this variable ensures that the time step of the simulation is always consistent at the end of the simulation. If one uses a time step that eventually leads exactly to the end time of the simulation this variable does not do anything. However, if adaptive time stepping is used or the end time is not exactly reached when using certain fixed time step, this flag ensures that the simulation does not change the last time step to reach the end time. For example, if your end time is 20, and you have a time step that leads to a last iteration until 20.1, all your results will be outputted until 20.1. If you wish to have exactly 20, you need to set this flag to ``false``. 

--------------------
Log file parameters
--------------------

* ``log frequency``: frequency at which information is written in the terminal.

* ``log precision``: number of significant digits used when writting in the terminal.

--------------------------------
Paraview output file parameters
--------------------------------

* ``output path``: directory for the output files.

* ``output name``: prefix for the Paraview output files (``.pvd`` / ``.vtu`` / ``.pvtu``)

  .. important::
	  Lethe saves the simulation results in the Paraview format. For every iteration, one or more ``.vtu`` are produced, which are indexed by a single ``.pvtu`` file. A single ``.pvd`` file linking all iterations together is also generated. Use the open-source software `Paraview <https://www.paraview.org/>`_ to visualize them.

* ``output control``: control for the output of the simulation results. The available options are: 
	* ``iteration``: results will be outputted at constant iteration frequency. 
	* ``time`` : results will be outputted based on time parameters (specific times or time frequency). The results can also be outputted for certain time interval.

* ``output frequency``: controls after which number of iterations the results are written. This parameter is only used when ``set output control = iteration``.

  .. tip::
	  If ``set output frequency = 0``, no output file will be generated. This is the only way to prevent the generation of output files.

* ``output time frequency``: controls the time frequency when the results are written, e.g., if set to 1, paraview files will be outputted every unit of time. This parameter is only used when ``set output control = time``.

* ``output times``: allows to specify specific times for the output of ``.pvd`` / ``.vtu`` files. This parameter is only used when ``set output control = time``. As an example, one can output files only at 5 seconds, by setting ``set output times = 5`` or at multiple specific times separating the values with commas: ``set output times = 5, 14``.

  .. warning::
	  Since it is possible that the times specified in the interval or in specific output times do not correspond to the time of specific iterations, Lethe will always write the Paraview files before and after the time specified. Furthermore, it is **not possible** to use both ``output times`` and ``output time frequency`` at the same time. For the ``output times`` to work, the value of ``output time frequency`` must be set to -1, which is the default value for the parameter. 

* ``output time interval``: Only writes the ``.pvd`` / ``.vtu`` files when the simulation time is within the closed interval defined by the ``output time interval``. Default values are 0s and 1.7e308s. Used for both ``iteration`` and ``time`` output control.

* ``group files``: number of ``.vtu`` files generated in a parallel simulation. 

  .. tip::
	  This parameter reduces the number of files generated when the simulation is run with a large number of processors. ``set group files = 1`` ensures that a single ``.vtu`` file will be generated. In this case, the file is written using MPI IO functionalities.

	  The value for this parameter should always be a compromise between keeping a low number of files but preventing excessive MPI communications. We have found that the default value of 1 does not have a significant impact except in very large simulations.

	  .. warning::
	  	As soon as the size of the output ``.vtu`` files reaches 3 Gb, it is preferable to start splitting them into multiple smaller files as this may lead to corrupted files on some file systems.

* ``output boundaries``: controls if the boundaries of the domain are written to a file. This will write additional ``.vtu`` files made of the contour of the domain. 

  .. tip::
	  This is particularly useful for the visualisation of 3D flows with obstacles or objects.

* ``subdivision``: sub-division of the mesh cells to enable visualisation of high-order elements with Paraview. 

  .. tip::
	  Generally, we advise to use a subdivision level of :math:`(n)` for interpolation order of :math:`n`. For example, a Q2-Q1 interpolation could be visualized with ``set subdivision = 2``.
