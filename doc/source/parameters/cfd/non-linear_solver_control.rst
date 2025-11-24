=================
Non-linear Solver
=================

The Navier-Stokes equations (and others) are non-linear equations. The parameters in ``subsection non-linear solver``, whose default values are given in the text block below, control the non-linear solver used within Lethe. Lethe supports different physics (``fluid dynamics``, ``VOF``, ``heat transfer``, ``cahn hilliard`` and ``tracer``) and it is possible to specify non-linear solver parameters for each of them. The ``VOF algebraic interface reinitialization`` subequation has also its own subsection.

In the example below, only ``fluid dynamics`` is shown but the same block can be used for other physics.

.. code-block:: text

  subsection non-linear solver
    subsection fluid dynamics
      # Non-linear solver that will be used
      set solver                       = newton

      # State whether information from the non-linear solver should be printed
      set verbosity                    = verbose

      # For the inexact_newton solver, sets the tolerance to re-assemble the Jacobian matrix
      set matrix tolerance             = 0.1

      # For the inexact_newton solver, carry jacobian matrix over to the new non-linear problem
      set reuse matrix                 = false

      # For the newton solver, it allows to reuse the preconditioner for the non linear iterations
      set reuse preconditioner         = false

      # For the kinsol_newton solver
      set kinsol_strategy              = line_search

      # Newton solver tolerance
      set tolerance                    = 1e-6

      # Tolerance for the acceptation of the step
      set step tolerance               = 0.9

      # Maximum number of Newton Iterations
      set max iterations               = 10

      # Number of digits displayed when showing residuals
      set residual precision           = 4

      # State if the RHS must be calculated at the beginning of every newton iteration
      set force rhs calculation        = false

      # Force the simulation to stop and throw an error if a non-linear solution has failed to converge
      set abort at convergence failure = false
    end
  end

* The ``solver`` parameter enables to choose the nonlinear solver used. Currently, Lethe supports three non-linear solvers:
	* ``newton`` solver (default parameter value), a Newton-Raphson solver which recalculates the Jacobian matrix at every iteration (see the Theory Documentation).
	* ``inexact_newton`` solver, a Newton-Raphson solver where the Jacobian matrix is reused between iterations.
		*  ``matrix tolerance`` parameter sets the tolerance to re-assemble the Jacobian matrix. If the residual after a newton step :math:`<` ``matrix tolerance`` :math:`\times` the previous residual, that iteration is considered sufficient and the Newton iteration will keep using the same jacobian matrix.
		* Setting ``reuse matrix = true`` enables the usage of the same Jacobian matrix for the following non-linear problem.

	.. tip::
		The ``inexact_newton`` solver, along with ``reuse matrix = true`` can be worthwhile in transient simulations with a small time-step. The goal is to seek a compromise between the cost of assembling the matrix and the preconditioner versus the cost of solving the linear system of equations.

	* ``kinsol_newton`` solver, that uses the Newton-Raphson solver through deal.II, as implemented in the `Sundials library <https://computing.llnl.gov/projects/sundials/kinsol>`_. This solver has an internal algorithm that decides whether to reassemble the Jacobian matrix or not. This non-linear solver is still being tested.
		* ``kinsol strategy`` parameter enables to choose the strategy that will be used by the kinsol newton solver, and can be ``line_search`` (default value), ``normal_newton``, ``fixed_point`` or ``picard``.
* The ``verbosity`` option enables the display of the residual at each non-linear iteration, to monitor the progress of the non-linear iterations.

.. note::
	The residual should decrease rapidly between Newton iterations.
	Example of a ``set verbosity = verbose`` output:
	
	.. code-block:: text

		Newton iteration: 0  - Residual:  429.541
				alpha = 1      res = 0.3705
		Newton iteration: 1  - Residual:  0.3705
				alpha = 1      res = 0.0001044
		Newton iteration: 2  - Residual:  0.0001044
				alpha = 1      res = 1.12e-08
		Newton iteration: 3  - Residual:  1.12e-08
				alpha = 1      res = 1.547e-12

* The ``step tolerance`` parameter controls how much the L2 norm of the residual must decrease to proceed to the next non-linear step. If the ``new_residual``:math:`<` ``old_residual``:math:`\times` ``step tolerance``, then a Newton iteration is accepted. If this condition is not reached, then a relaxation of the step is applied (increasing the ``alpha`` parameter, as printed on the terminal if ``set verbosity = verbose``) until this condition is reached.
* The ``tolerance`` parameter controls the value under which the residual must be to proceed to the next iteration.

.. hint::
	The ``tolerance`` parameter is directly linked to the numerical convergence of the simulation, but also to the computational cost (number of Newton iteration).

	For simple simulations, the tolerance can be set quite low, for instance ``set tolerance = 1e-12``. However, such a tolerance can be impossible to attain for more complex simulations : the step tolerance of the non-linear solver can be increased, for instance ``set tolerance = 1e-4``

.. warning::

    The residual of the non-linear solver is affected by the ``rescale residual`` parameter in the ``subsection linear solver``. When the parameter is activated, the L2 norm of the residual is divided by the square root of the domain's volume. Hence, the tolerance of the non-linear solver must be adapted accordingly. For further details, please refer to the `Linear Solver parameters <./linear_solver_control.html>`_ section of this guide.

* The ``max iterations`` parameter sets a hard limit to the number of Newton iterations, even if the ``tolerance`` is not reached.

.. warning::
	Be careful to always set an absolute tolerance for the linear solver that is below the tolerance of the non-linear solver. Otherwise, you might find that it is impossible to converge because the linear system of equation is solved with insufficient accuracy.

* The ``residual precision`` parameter enables to change the number of digits displayed when showing residuals (with ``set verbosity = verbose``).
* The ``force_rhs_calculation``: Force RHS recalculation at the beginning of every non-linear steps, This is required if there is a fixed point component to the non-linear solver that is changed at the beginning of every newton iteration. This is notably the case of the sharp edge method. The default value of this parameter is false.
* The ``abort at convergence failure`` allows the user to stop the simulation and throw an error if the non-linear solver has failed to converge. Setting ``abort at convergence failure = true`` will enable this feature. This is generally useful when launching a large batch of simulation to quickly identify which one have failed.
* The ``reuse preconditioner = true`` allows the simulation to use the same preconditioner between Newton iterations when using the Newton solver. This can reduce the overall time depending on the problem, and it is especially useful for the ``lethe-fluid-matrix-free`` application.
