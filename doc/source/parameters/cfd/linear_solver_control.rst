=============
Linear Solver
=============

In this subsection, the control options of the linear solvers are specified. The default values for the linear solver parameters are given in the text box below. Lethe supports different physics (``fluid dynamics``, ``VOF``, ``heat transfer``, ``cahn hilliard`` and ``tracer``) and it is possible to specify linear solver parameters for each of them. In the example below, only ``fluid dynamics`` is used, however, the same block can be used for other physics.

.. seealso::
	For further understanding about the linear solvers used, the preconditioners supported and all parameters, see the :doc:`../../theory/fluid_dynamics/linear_solvers` theory section.

.. code-block:: text

  subsection linear solver
    subsection fluid dynamics
      # Iterative solver for the linear system of equations
      set method                           = gmres

      # State whether information from the linear solver should be printed
      set verbosity                        = verbose

      # Linear solver minimum residual
      set minimum residual                 = 1e-12

      # Linear solver residual
      set relative residual                = 1e-3

      # Maximum solver iterations
      set max iters                        = 1000

      # Force the linear solver to continue even if it fails
      set force linear solver continuation = false

      # Maximum number of krylov vectors for GMRES solver
      set max krylov vectors               = 100

      # Set type of preconditioner for the iterative solver
      set preconditioner                   = ilu
    end
  end


* The ``method`` parameter enables to choose an iterative solver for the linear system of equations. Lethe currently supports the following core solution strategies:
	* ``gmres`` (default parameter value), a GMRES iterative solver with ILU preconditioning or AMG preconditioning with an ILU coarsener and smoother.
	* ``bicgstab``, a BICGSTAB iterative solver with ILU preconditioning.
	* ``direct``, a direct solver using `TrilinosWrappers::SolverDirect <https://www.dealii.org/current/doxygen/deal.II/classTrilinosWrappers_1_1SolverDirect.html>`_.

	.. hint::
		Which solver is the most efficient for a given problem?
		
		* For steady-state problems (2D and 3D):
			* coarse meshes: ``gmres``/``bicgstab`` solver with ``ilu`` preconditioner.
			* fine meshes (from :math:`\sim 1` M cells): ``gmres`` solver with ``amg`` preconditioner.
		* For transient problems, this is much more problem dependent, but generally:
			* relatively low :math:`\text{CFL}` condition: ``gmres`` solver with ``ilu`` preconditioner with a low fill level, for instance ``set ilu preconditioner fill = 0``. This applies even if the mesh is very large (:math:`>10` M cells), because the transient version of the system of equations is much easier to solve.
			* large :math:`\text{CFL}` condition (:math:`\text{CFL}>10`) and/or very large mesh: ``gmres`` solver with ``amg`` preconditioner may become preferable.

	.. caution:: 
		The use of ``direct`` solver should be in general avoided as it is not efficient for large problems. It can be, however, used for debugging purposes or the development of new features.


* The ``verbosity`` option enables the display of the residual at each non-linear iteration, to monitor the progress of the linear iterations.

.. admonition:: Example of a ``set verbosity = verbose`` output:

	.. code-block:: text

		-Tolerance of iterative solver is : 0.0429541
		-Iterative solver took : 11 steps 
		-Tolerance of iterative solver is : 3.62082e-05
		-Iterative solver took : 16 steps 
		-Tolerance of iterative solver is : 1.05775e-08
		-Iterative solver took : 17 steps 
		-Tolerance of iterative solver is : 1.00205e-12
		-Iterative solver took : 16 steps 
		-Tolerance of iterative solver is : 1e-13
		-Iterative solver took : 5 steps 
		-Tolerance of iterative solver is : 1e-13
		-Iterative solver took : 0 steps 


* The ``minimum residual`` for the linear solver.

* The ``relative residual`` for the linear solver.

.. tip::
	A good rule of thumb is to set the linear solver ``minimum residual`` at least :math:`10` times (preferably :math:`100` times) smaller than the :doc:`non-linear_solver_control`, ``tolerance`` parameter, and keep the relative residual reasonable, for instance ``set relative residual = 1e-3``. To lower the computational cost for more complex simulations, it can be lowered to ``set relative residual = 1e-4``.

* The ``max iters`` puts a hard stop on the number of solver iterations (number of steps printed when ``set verbosity = verbose``).

.. tip::
	If ``max iters`` is reached, the code will throw this type of message: 
	
	.. code-block:: text
	
		GMRES solver failed! Trying with a higher preconditioner fill level.

	meaning that the code increases the preconditioner fill (see tip on default values below) in order to converge within the number of solver iterations. If you encounter this, consider increasing the ``max iters`` or adjusting other parameters, for example increasing ``max krylov vectors``.

* ``force linear solver continuation`` when set to ``true``, forces the linear solver to continue, even if the ``minimum residual`` is not reached. Only available for ``gmres`` and ``bicgstab`` solvers within the ``lethe-fluid`` application.

.. warning::
	With this mode on, errors on the linear solver convergence are not thrown. Forcing the solver to continue can be useful for debugging purposes if a given iteration is hard to pass, but use it with caution!

* ``max krylov vectors`` sets the maximum number of krylov vectors for ``gmres`` solver with ``ilu`` and ``amg`` preconditioners.

.. tip::
	Consider using ``set max krylov vectors = 200`` for complex simulations with convergence issues. 

* ``preconditioner`` sets the type of preconditioning used for the linear solver. It can be either ``ilu`` for an Incomplete LU decomposition or ``amg`` for an Algebraic Multigrid. 

.. warning::
    Currently, the ``lethe-fluid-sharp`` solver makes it almost impossible to reach convergence with the ``amg`` preconditioner. Therefore, it is recommended to use ``ilu`` instead, even for fine meshes. In addition, the ``VOF``, ``heat transfer``, ``cahn hilliard`` and ``tracer`` physics only support ``ilu``.

.. caution:: 
		Be aware that the setup of the ``amg`` preconditioner is very expensive and does not scale linearly with the size of the matrix. As such, it is generally preferable to minimize the number of assembly of such preconditioner. This can be achieved by using the ``inexact newton`` for the nonlinear solver (see :doc:`non-linear_solver_control`).

In addition to the method parameters, one can also set specific parameters for each of the preconditioners by adding specific lines inside of the specific physics subsection:

-------------------
ILU preconditioner
-------------------

.. code-block:: text

    # ILU preconditioner fill
    set ilu preconditioner fill               = 0

    # ILU preconditioner tolerance
    set ilu preconditioner absolute tolerance = 1e-12

    # ILU relative tolerance
    set ilu preconditioner relative tolerance = 1.00

.. tip::
	The default values for these parameters are good starting values. 

	For each iteration of the linear solver (at the beginning of which the tolerance of the iterative solver is computed, as printed if ``set verbosity = verbose``), the chosen solver starts by using the ``preconditioner fill`` given in the parameter file. If for any reason the linear solver would have crashed, it will restart with a fill level increased by 1. This restart process will happen up to a maximum of 3 times, after which it will let the solver crash. 

	Hence, for complex simulations, if you get at almost every linear iteration the message:

	.. code-block:: text
	
		GMRES solver failed! Trying with a higher preconditioner fill level. New fill = ...

	and it does not disappear when increasing ``max iters``, increasing the ``ilu preconditioner fill`` in the ``.prm`` file will make the computation slightly faster.

-------------------
AMG preconditioner
-------------------

.. code-block:: text

    # AMG preconditioner ILU smoother/coarsener fill
    set amg preconditioner ilu fill               = 0

    # AMG preconditioner ILU smoother/coarsener absolute tolerance
    set amg preconditioner ilu absolute tolerance = 1e-12

    # AMG preconditioner ILU smoother/coarsener relative tolerance
    set amg preconditioner ilu relative tolerance = 1.00

    # AMG aggregation threshold
    set amg aggregation threshold                 = 1e-14

    # AMG number of cycles
    set amg n cycles                              = 1

    # AMG w cycling. If this is set to true, W cycling is used. Otherwise, V cycling is used.
    set amg w cycles                              = false

    # AMG smoother sweeps
    set amg smoother sweeps                       = 2

    # AMG smoother overlap
    set amg smoother overlap                      = 1

.. seealso::
	For more information about the ``amg`` preconditioner parameters, the reader is referred to the deal.II documentation for the `AMG preconditioner <https://www.dealii.org/current/doxygen/deal.II/classTrilinosWrappers_1_1PreconditionAMG.html>`_ and its `Additional Data <https://www.dealii.org/current/doxygen/deal.II/structTrilinosWrappers_1_1PreconditionAMG_1_1AdditionalData.html>`_.