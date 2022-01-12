
Analytical solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the problem being simulated has a known analytical solution, or an exact solution is imposed (manufactured solutions) it can be added in this section. The default parameters are:

.. code-block:: text

   subsection analytical solution
    set enable                = false
    set verbosity             = quiet
    set filename              = L2Error
    subsection uvwp
      set Function expression = 0; 0; 0
        or
      set Function expression = 0; 0; 0; 0
    end
    subsection temperature
      set Function expression = 0
    end
    subsection tracer
      set Function expression = 0
    end
    subsection phase
      set Function expression = 0
    end
   end

* The ``enable`` parameters sets if the problem has an analytical solution and enables the calculation of the analytical solution and the norm of the L2 error.

* The ``verbosity`` parameter enables printing the L2 error after each mesh refinement if set to ``verbose``, if ``enable`` is ``true``.

* The ``filename`` parameter sets the file name to output the L2 error norm if ``enable`` is ``true``.

If there is an analytical solution for velocity and pressure, enter the ``uvwp`` subsection.

* The ``Function expression`` parameters sets the expression for the analytical solution in regards to u, v and p for a 2D simulation and to u, v, w and p for a 3D simulation.

If there is an analytical solution for the fluid's temperature, enter the ``temperature`` subsection.

* The ``Function expression`` parameters sets the expression of the temperature.

If there is an analytical solution for a tracer, enter the ``tracer`` subsection.

* The ``Function expression`` parameters sets the expression of the tracer.

If there is an analytical solution for the fluid's phasee, enter the ``phase`` subsection.

* The ``Function expression`` parameters sets the expression of the phase.

.. note:: 
    The variables *x*, *y*, *z* (3D) and *t* (time-dependant) can be used in the function expressions.

In all four last subsections, you can add a ``Function constant`` parameter that will act as a constant in the ``Function expression``.

You can add a ``Function constants`` parameter that will act as a constant in the ``Function expression``. 

Ex.

.. code-block:: text

   subsection analytical solution
    set enable                = true
    set verbosity             = true
    set filename              = L2Error
    subsection uvwp
      set Function constants = A=2.0
      set Function expression = A*y; -A*x; 0
    end
   end
   
.. note:: 
    A :math:`\pi` variable is already defined as both ``pi`` and ``Pi``

Your function expression can also contain common functions such as :math:`\sin`, :math:`\cos` as well as ``if`` statements.

Ex.

.. code-block:: text

   subsection analytical solution
    set enable                = true
    set verbosity             = true
    set filename              = L2Error
    subsection phase
      set Function expression = if(sin(x) > pi, 1, 0)
    end
   end

