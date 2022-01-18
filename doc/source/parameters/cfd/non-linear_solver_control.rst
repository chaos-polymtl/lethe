Non-linear Solver Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The Navier-Stokes equations (and other) are non-linear equations. The parameters in ``subsection non-linear solver``, whose default values are given in the text block below, control the non-linear solver used within Lethe.

.. code-block:: text

	subsection non-linear solver
	  # Non-linear solver that will be used. Choices is newton.
	  set solver             = newton

	  # State whether from the non-linear solver should be printed Choices are
	  # <quiet|verbose>.
	  set verbosity          = verbose

	  # Newton solver tolerance
	  set tolerance          = 1e-6

	  # Tolerance for the acceptation of the step
	  set step tolerance     = 0.9

	  # Maximum number of Newton Iterations
	  set max iterations     = 10
	end

* The ``solver`` parameter enables to choose the nonlinear solver used. Currently, the only solver implemented is ``newton``, a Newton-Raphson solver which recalculates the jacobian matrix at every iteration (see the Theory Documentation).
* The ``verbosity`` options enables to display the residual at each non-linear iteration, to monitor how the non-linear iterations are progressing.
.. note::
	The residual should decrease rapidly between Newton iterations.
	Exemple of a ``set verbosity = verbose`` output:
	
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

	For simple simulations, the step tolerance can be quite low, for instance ``set tolerance = 1e-12``. However, such a tolerance can be impossible to attain for more complex simulations : the step tolerance of the non-linear solver can be increased, for instance ``set tolerance = 1e-4``
* The ``max iterations`` parameter sets a hard limit to the number of Newton iterations, even if the ``tolerance`` is not reached.
