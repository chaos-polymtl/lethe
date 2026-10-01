===========
Source Term
===========

If the problem being simulated has a source, it can be added in this section. The default parameters are:

.. code-block:: text

  subsection source term
    subsection fluid dynamics
      set Function expression = 0; 0; 0 #In 2D
      set Function expression = 0; 0; 0; 0 #In 3D
      set enable              = true
    end

    subsection heat transfer
      set Function expression = 0
    end

    subsection tracer
      set Function expression = 0
    end

    subsection cahn hilliard
      set Function expression = 0; 0
    end
  end

.. tip:: 
  ``Function expression``, used in this subsection (but also in :doc:`./initial_conditions`, :doc:`./analytical_solution`, namely), give access to several tools:
  
  * define ``Function constants``
  * use :math:`\pi` variable, as ``pi`` or ``Pi``
  * use common functions such as :math:`\sin`, :math:`\cos` 
  * use ``if`` statements

  Check the :ref:`ex function` for further help.

* ``subsection fluid dynamics``: defines the parameters for a Navier-Stokes source term. This source term is defined by a ``Function expression`` and can depend on both space and time.

  * In 2D, the first two terms are the source terms for  the :math:`x`, :math:`y` component of the momentum equation. The third term is the mass source term. 
  * In 3D, the first three terms are for the :math:`x`, :math:`y` and :math:`z` component of the momentum equation and the fourth term is for the mass source term.

  .. tip::

	For ``subsection fluid dynamics``, each term can depend on both space (``x``, ``y`` and, if 3D, ``z``) and time (``t``). See :ref:`ex function`.

  .. tip::

	If you are using the ``lethe-fluid-matrix-free`` application the usage of a source term significantly affects performance. If you are not using it, we advice you to disable it explicitly by setting ``enable = false``.

* ``subsection heat transfer``: defines the parameters for a heat source term. This source term is defined by a ``Function expression`` and can depend on both space (``x``, ``y`` and, if 3D, ``z``) and time (``t``). See :ref:`ex function`.

* ``subsection tracer``: defines the parameters for the a source term for a tracer. This source term is defined by a ``Function expression`` and can depend on both space (``x``, ``y`` and, if 3D, ``z``) and time (``t``). See :ref:`ex function`.

* ``subsection cahn hilliard``: defines the parameters for a source term in the Cahn-Hilliard equations. This source term is defined by a ``Function expression`` and can depend on both space (``x``, ``y`` and, if 3D, ``z``) and time (``t``). Both the phase order parameter (first component) and chemical potential (second component) can have source terms, hence the two components. See :ref:`ex function`.


.. _ex function:

Examples of Function Expression
--------------------------------

CFD source term with ``Function constants``:

.. code-block:: text

    subsection fluid dynamics
      set Function constants = A=2.0
      set Function expression = A*y; -A*x; 0
    end

CFD source term varying in time:

.. code-block:: text

    subsection fluid dynamics
        set Function expression = 0; -10*cos(2*pi*t); 0
    end

Heat transfer source term with ``if()`` condition:

.. code-block:: text

    subsection heat transfer
      set Function expression = if(sin(x) > pi, 1, 0)
	# if ( condition , value if true , value if false )
    end

.. note:: 
  The first parameter in the ``if()`` function is the statement. If this statement is :
    * ``true``, then the function expression takes the second parameter as value
    * ``false``, the function expression takes the third parameter as value. 

  In this example, the heat source term will vary within the calculation domain.

CFD source term with ``Function constants``:

.. code-block:: text

    subsection fluid dynamics
      set Function constants = A=2.0, B=1.0
      set Function expression = A*y; -B*x; 0
    end

