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

    subsection electromagnetics
      subsection current density
        set Function expression = 0; 0; 0; 0; 0; 0 #In 3D
      end
      subsection surface current density
        set Function expression = 0; 0; 0; 0; 0; 0 #In 3D
      end
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

* ``subsection electromagnetics``: defines the imposed electric currents of the time-harmonic Maxwell equations. Both currents are defined by a ``Function expression`` and can depend on both space (``x``, ``y`` and ``z``) and time (``t``). See :ref:`ex function`. Their six components are the real parts of the :math:`x`, :math:`y` and :math:`z` components, followed by their imaginary parts. The currents are phasors: the physical current is :math:`\mathbf{j}(\mathbf{x},t) = \Re(\mathbf{J}(\mathbf{x}) e^{-i\omega t})`, so the imaginary part sets the phase of a current relative to the other sources (e.g., the waveguide ports). It can be left to zero for a single source. Both currents are zero by default and are imposed as given in the non-dimensional system solved by the time-harmonic Maxwell solver.

  * ``subsection current density``: volume current density :math:`\mathbf{J}`, which enters the Ampère law as :math:`\nabla \times \mathbf{H} + i\omega\varepsilon\mathbf{E} = \mathbf{J}`.
  * ``subsection surface current density``: surface current density :math:`\mathbf{K}` of a current sheet, which imposes the jump :math:`\hat{\mathbf{n}} \times (\mathbf{H}_2 - \mathbf{H}_1) = \mathbf{K}` of the tangential magnetic field across the faces of the sheet, where :math:`\hat{\mathbf{n}}` points from side 1 to side 2. Only the tangential part of :math:`\mathbf{K}` is used. The surface current density is only applied on the interior faces of the mesh, so the current sheet must be aligned with the faces of the cells and the function must only be non-zero on these faces, e.g., ``if(abs(z-1)<1e-8, 1, 0); 0; 0; 0; 0; 0`` for a sheet in the plane :math:`z=1`. Sources on the boundaries of the domain are defined with the ``magnetic field``, ``impedance boundary`` or ``waveguide port`` boundary conditions instead.

  .. warning::

	Non-zero imposed currents have not been validated yet.


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

