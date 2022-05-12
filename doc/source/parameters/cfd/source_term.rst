Source term
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the problem being simulated has a source, it can be added in this section. The default parameters are:

.. code-block:: text

   subsection source term
    set enable                = true
    subsection xyz
      set Function expression = 0; 0; 0
        or
      set Function expression = 0; 0; 0; 0
    end
    subsection heat transfer
      set Function expression = 0
    end
    subsection tracer
      set Function expression = 0
    end
   end

* The ``enable`` parameter is set to true if the problem has a source term and enables its calculation.

If there is a Navier-Stokes source term, enter the ``xyz`` subsection.

* The ``Function expression`` parameters sets the expression for the source term. In 2D, the first two terms are the source terms for  the *x*, *y* component of the momentum equation. The third term is the mass source term. In 3D, the first three terms are for the *x*, *y* and *z* component of the momentum equation and the fourth term is for the mass source term.

If there is a heat source term, enter the ``heat transfer`` subsection.

* The ``Function expression`` parameter sets the expression for the heat transfer.

If there is an source term for a tracer, enter the ``tracer`` subsection.

* The ``Function expression`` parameter sets the expression for the tracer term.

.. note:: 
    The variables *x*, *y*, *z* (3D) and *t* (time-dependant) can be used in the function expressions.

You can add a ``Function constants`` parameter that will act as a constant in the ``Function expression``. 

Ex.

.. code-block:: text

   subsection source term
    set enable                = true
    subsection xyz
      set Function constants = A=2.0
      set Function expression = A*y; -A*x; 0
    end
   end
   
.. note:: 
    A :math:`\pi` variable is already defined as both ``pi`` and ``Pi``

Your function expression can also contain common functions such as :math:`\sin`, :math:`\cos` as well as ``if`` statements.

Ex.

.. code-block:: text

   subsection source term
    set enable                = true
    subsection heat transfer
      set Function expression = if(sin(x) > pi, 1, 0)
    end
   end

.. note:: 
   The first parameter in the ``if()`` function is the statement. If this statement is *true*, then the function expression takes the second parameter as value. If this statement is *false*, the function expression takes the third parameter as value. In this example, the heat source term will vary within the calculation domain.
