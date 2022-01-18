Linear Solver Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In these subsections, the control options of linear solvers are specified. These control options may include solution method, maximum number of iterations, tolerance, residual precision and preconditioner information. The default values for these parameters are given in the text box below.

.. code-block:: text

	subsection linear solver
	  # State whether output from solver runs should be printed. Choices are
	  # <quiet|verbose>.
	  set verbosity                                 = verbose

	  # The iterative solver for the linear system of equations. Choices are
	  # <gmres|bicgstab|amg|direct>. gmres is a GMRES iterative solver with ILU
	  # preconditioning. bicgstab is a BICGSTAB iterative solver with ILU
	  # preconditioning. amg is GMRES + AMG preconditioning with an ILU coarsener
	  # and smoother. On coarse meshes, the gmres/bicgstab solver with ILU
	  # preconditioning is more efficient. As the number of mesh elements
	  # increase, the amg solver is the most efficient. Generally, at 1M elements,
	  # the amg solver always outperforms the gmres or bicgstab
	  set method                                    = gmres



	  # amg aggregation threshold
	  set amg aggregation threshold                 = 1e-14

	  # amg number of cycles
	  set amg n cycles                              = 1

	  # amg preconditioner ilu smoother/coarsener absolute tolerance
	  set amg preconditioner ilu absolute tolerance = 1e-12

	  # amg preconditioner ilu smoother/coarsener fill
	  set amg preconditioner ilu fill               = 1

	  # amg preconditioner ilu smoother/coarsener relative tolerance
	  set amg preconditioner ilu relative tolerance = 1.00

	  # amg smoother overlap
	  set amg smoother overlap                      = 1

	  # amg smoother sweeps
	  set amg smoother sweeps                       = 2

	  # amg w cycling. If this is set to true, W cycling is used. Otherwise, V
	  # cycling is used.
	  set amg w cycles                              = false

	  # Ilu preconditioner tolerance
	  set ilu preconditioner absolute tolerance     = 1e-6

	  # Ilu preconditioner fill
	  set ilu preconditioner fill                   = 1

	  # Ilu relative tolerance
	  set ilu preconditioner relative tolerance     = 1.00

	  # Maximum solver iterations
	  set max iters                                 = 1000

	  # Maximum number of krylov vectors for GMRES solver
	  set max krylov vectors                        = 30


	  # Linear solver minimum residual
	  set minimum residual                          = 1e-8

	  # Linear solver residual
	  set relative residual                         = 1e-3
	end

Lethe currently supports three core solution strategy or ``method`` : 
* GMRES solver + ILU preconditioning. This is set by ``method = gmres``
* GMRES solver + AMG preconditioning. This is set by ``method = amg``
* BiCGStab solver + ILU preconditioning. This is set by ``method = bicgstab``

For all methods, a number of parameters to control the preconditioner can be set.

On coarse meshes, the gmres/bicgstab solver with ILU preconditioning is more efficient. As the number of mesh elements increase, the amg solver is the most efficient. Generally, at 1M elements,
the amg solver always outperforms the gmres or bicgstab for 2D steady-state problems.


