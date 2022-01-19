Linear Solver Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In these subsections, the control options of linear solvers are specified. These control options may include solution method, maximum number of iterations, tolerance, residual precision and preconditioner information. The default values for these parameters are given in the text box below.

.. code-block:: text

	subsection linear solver

	  # Iterative solver for the linear system of equations
	  set method		= gmres

	  # State whether information from the linear solver should be printed
	  set verbosity		= verbose

	  # Linear solver minimum residual
	  set minimum residual  = 1e-12

	  # Linear solver residual
	  set relative residual = 1e-3

	  # Maximum solver iterations
	  set max iters         = 1000

	  # Force the linear solver to continue even if it fails
	  set force linear solver continuation = false

	  #---------------------------------------------------
	  # GMRES solver parameters
	  #---------------------------------------------------
	  # Maximum number of krylov vectors for GMRES solver
	  set max krylov vectors = 100

	  #---------------------------------------------------
	  # ILU parameters for GMRES and BICGSTAB
	  #---------------------------------------------------
	  # Ilu preconditioner tolerance
	  set ilu preconditioner absolute tolerance = 1e-12

	  # Ilu preconditioner fill
	  set ilu preconditioner fill = 0

	  # Ilu relative tolerance
	  set ilu preconditioner relative tolerance = 1.00

	  #---------------------------------------------------
	  # AMG solver parameters
	  #---------------------------------------------------
	  # amg preconditioner ILU smoother/coarsener fill
	  set amg preconditioner ilu fill = 0

	  # amg preconditioner ilu smoother/coarsener absolute tolerance
	  set amg preconditioner ilu absolute tolerance = 1e-12

	  # amg preconditioner ilu smoother/coarsener relative tolerance
	  set amg preconditioner ilu relative tolerance = 1.00

	  # amg aggregation threshold
	  set amg aggregation threshold = 1e-14

	  # amg number of cycles
	  set amg n cycles              = 1

	  # amg w cycling. If this is set to true, W cycling is used. Otherwise, V cycling is used.
	  set amg w cycles              = false

	  # amg smoother sweeps
	  set amg smoother sweeps       = 2

	  # amg smoother overlap
	  set amg smoother overlap      = 1
	end


* The ``method`` parameter enables to choose an iterative solver for the linear system of equations. Lethe currently supports these core solution strategies :
	* ``gmres`` (default parameter value), a GMRES iterative solver with ILU preconditioning.
	* ``amg``, a GMRES iterative solver with AMG preconditioning and an ILU coarsener and smoother.
	* ``bicgstab``, a BICGSTAB iterative solver with ILU preconditioning.
	* ``direct``, a direct solver using `TrilinosWrappers::SolverDirect <https://www.dealii.org/current/doxygen/deal.II/classTrilinosWrappers_1_1SolverDirect.html>`_. 

	.. tip:: 
		On coarse meshes, the ``gmres``/``bicgstab`` solver with ILU preconditioning is more efficient. 

		As the number of mesh elements increase, the ``amg`` solver is the most efficient. Generally, at 1M elements, the ``amg`` solver always outperforms the ``gmres`` or ``bicgstab`` for 2D steady-state problems.
		
		The use of ``direct`` solver should be avoided for 3D problems.

.. note::
	For further understanding about iterative solver methods, the interested reader is referred to the Theory Documentation.


* The ``verbosity`` option enables the display of the residual at each non-linear iteration, to monitor the progress of the linear iterations.

.. note::
	Example of a ``set verbosity = verbose`` output:

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


For all methods, a number of parameters to control the preconditioner can be set.
