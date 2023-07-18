===================
Analytical Solution
===================

If the problem being simulated has a known analytical solution, or an exact solution is imposed (manufactured solutions) it can be added in this section. The default parameters are:

.. code-block:: text

   subsection analytical solution
    set enable                = false
    set verbosity             = quiet
    set filename              = L2Error
    subsection uvwp
      set Function expression = 0; 0; 0  # In 2D : u;v;p
        or
      set Function expression = 0; 0; 0; 0  #In 3D u;v;w;p
    end
    subsection temperature
      set Function expression = 0
    end
    subsection tracer
      set Function expression = 0
    end
    subsection VOF
      set Function expression = 0
    end
    subsection cahn hilliard
      set Function expression = 0; 0 # phase order parameter; chemical potential
    end
   end

* The ``enable`` parameter is set to true if the problem has an analytical solution. This enables the calculation of the analytical solution and the :math:`L^2` norm of the error.

The :math:`L^2` norm of the error is calculated as

.. math::
 L^2 = \int_\Omega (u-u_a)^2 \mathrm{d} \Omega

where :math:`u` is the numerical solution and  :math:`u_a` is the analytical solution.



* The ``verbosity`` parameter enables printing the L2 error after each mesh refinement if it is set to ``verbose`` and if ``enable`` is ``true``.

* The ``filename`` parameter sets the file name to output the L2 error norm if ``enable`` is ``true``.

If there is an analytical solution for velocity and pressure, enter the ``uvwp`` subsection.

* The ``Function expression`` parameter sets the expression for the analytical solution in regards to :math:`u`, :math:`v` and :math:`p` for a 2D simulation and to :math:`u`, :math:`v`, :math:`w` and :math:`p` for a 3D simulation.

If there is an analytical solution for the fluid's temperature, enter the ``temperature`` subsection.

* The ``Function expression`` parameter sets the expression of the temperature.

If there is an analytical solution for a tracer, enter the ``tracer`` subsection.

* The ``Function expression`` parameter sets the expression of the tracer.

If there is an analytical solution for the VOF physics, enter the ``VOF`` subsection.

* The ``Function expression`` parameter sets the expression of the VOF indicator.

If there is an analytical solution for the Cahn-Hilliard physics, enter the ``cahn hilliard`` subsection.

* The ``Function expression`` parameter sets the expression of the phase order parameter and the chemical potential.

.. note:: 
    The variables *x*, *y*, *z* (3D) and *t* (time-dependant) can be used in the function expressions.

In all four last subsections, you can add a ``Function constants`` parameter that will act as a constant in the ``Function expression``.


Ex:

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

.. note:: 
   The first parameter in the ``if()`` function is the statement. If this statement is *true*, then the function expression takes the second parameter as value. If this statement is *false*, the function expression takes the third parameter as value. In this example, the analytical phase will vary within the calculation domain.

