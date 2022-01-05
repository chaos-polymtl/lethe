Simulation Control 
------------------

This subsection contains the general information of the simulation, including the time integration method and the general output names. It is the most commonly modified section for a simulation. A standard convention in Lethe is to keep this section at the top of the parameter file, since it is generally the most accessed one.

.. code-block:: text

    # Adaptative time-stepping <true|false>
    set adapt                        = false
  
    # Adaptative time step scaling
    set adaptative time step scaling = 1.1
  
    # The kind of method used to startup high order bdf methods. Choices are
    # <multiple step bdf|sdirk step|initial solution>.
    set bdf startup method           = multiple step bdf
  
    # Maximal number of vtu output files
    set group files                  = 1
  
    # Log frequency
    set log frequency                = 1
  
    # Display precision when writing to log
    set log precision                = 6
  
    # Maximum CFL value
    set max cfl                      = 1
  
    # The kind of solver for the linear system. Choices are
    # <steady|steady_bdf|bdf1|bdf2|bdf3|sdirk2|sdirk3>.
    set method                       = steady
  
    # Number of mesh adaptation (for steady simulations)
    set number mesh adapt            = 0
  
    # Output the boundaries of the domain along with their ID
    set output boundaries            = false
  
    # The control for the output of the simulation results. Results can be either
    # outputted at constant iteration frequency or at constant time
    set output control               = iteration
  
    # Output frequency
    set output frequency             = 1
  
    # File output prefix
    set output name                  = out
  
    # File output path
    set output path                  = ./
  
    # Output time
    set output time                  = 1
  
    # Scaling factor used in the iterations necessary to start-up the BDF
    # schemes.
    set startup time scaling         = 0.4
  
    # Tolerance at which the simulation is stopped
    set stop tolerance               = 1e-10
  
    # Subdivision of mesh cell in postprocessing
    set subdivision                  = 1
  
    # End time value of the simulation
    set time end                     = 1
  
    # Time step value
    set time step                    = 1.
  end

* The ``group files`` parameter controls the number of vtu files generated in a parallel simulation. This parameter allows to reduce the number of files generated when the simulation is run with a large number of processors. Setting group files to 1 ensures that there is a single vtu file generated. In this case, the file is written using MPI IO functionalities. The value for this parameter should always be a compromise between keeping a low number of files but preventing excessive MPI communications. We have found that the default value of 1 does not have a significant impact on performance on Compute Canada clusters.

* The ``method`` parameter controls the time-stepping method used. The available options are: 
    * ``steady`` (steady-state simulation)
    * ``steady-bdf`` (steady-state simulation using adjoint time stepping with a bdf1 scheme)
    * ``bdf1`` (1st order backward differentiation)
    * ``bdf2`` (2nd order backward differentiation)
    * ``bdf3`` (3rd order backward differentiation)
    * ``sdirk2`` (2nd order singly diagonally implicit Runge Kutta)
    * ``sdirk3`` (3rd order singly diagonally implicit Runge Kutta)

* The ``bdf startup method`` parameter controls the scheme used to start a high order bdf scheme.

* The ``output frequency`` controls after which number of iterations the ``.pvd`` / ``.vtu`` results are written. If ``output frequency = 0``, no ``.pvd`` / ``.vtu`` file will be written.

* The ``output name`` and ``output path`` control the name of the output files and their directory.

* The ``output boundaries`` controls if the boundaries of the domain are written to a file. This will write additional VTU files made of the contour of the domain. This is particularly useful for the visualisation of 3D flows with obstacles or objects.

* The ``subdivision`` parameter enables sub-division of the mesh cells to enable visualisation of high-order elements with Paraview. Generally, we advise to use a subdivision level of (n-1) for interpolation order of n. For example, a Q2-Q2 interpolation could be visualized with ``subdivision=1``.

* The ``time end`` parameter indicates at which value of the time to end the simulations if they are transient.

* The ``time step`` parameter controls the value of the time step.

* The ``adapt`` parameter controls if adaptive time-stepping is enabled. If adaptive time-stepping is enabled, the time-step will evolve to ensure that the 'max cfl' value is reached.

* The ``max cfl`` parameter controls the maximal value of the CFL condition that can be reached during the simulation. This parameter is only used when `adapt` is set to true.

* The ``adaptative time step scaling`` parameter controls the rate of increase of the time step value. The new time step value is fixed by ``adaptative time step scaling`` * ``previous value of the time step``

* The ``stop tolerance`` parameter controls the tolerance at which the adjoint time stepping steady state simulations (steady_bdf) stops. The adjoint time stepping will stop when the L2 norm of the initial residual is lower than ``stop tolerance`` at the start of a non-linear solution step.

* The ``log frequency`` parameter controls the frequency at which information is written to the log (the terminal).

* The ``log precision`` parameter controls the number of significant digits used when writing to the log (the terminal).
