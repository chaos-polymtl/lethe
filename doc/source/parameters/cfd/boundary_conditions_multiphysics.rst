==================================
Boundary Conditions - Multiphysics
==================================

This subsection's purpose is defining the boundary conditions associated to multiphysic problems. 

.. _heat transfer bc:

Heat Transfer
^^^^^^^^^^^^^

For heat transfer boundary conditions, the possible ``types`` are ``noflux`` (default), ``temperature`` and ``convection-radiation``.
The default parameters for ``temperature`` and ``convection-radiation`` are shown: 

.. code-block:: text

  subsection boundary conditions heat transfer
    set number         = 2
    set time dependent = false
    subsection bc 0
      set id    = 0
      set type  = temperature
      subsection value
        set Function expression = 0
      end
    end
    subsection bc 1
      set id         = 1
      set type       = convection-radiation
      subsection h
        set Function expression = 0
      end
      subsection Tinf
        set Function expression = 0
      end
      subsection emissivity
        set Function expression = 0
      end
    end
    set Stefan-Boltzmann constant = 0.000000056703
  end

* ``number``: This is the number of boundary conditions of the problem.

* ``time dependent`` specifies if a boundary condition is time-dependent (``true``) or steady (``false``). By default, this parameter is set to ``false``. This is here to improve the computational efficiency for transient cases in which the boundary conditions do not change.

.. warning::
    The ``number`` of boundary conditions must be specified explicitly. This is often a source of error.

.. note::
    The index in ``subsection bc ..`` must be coherent with the ``number`` of boundary conditions set: if ``number = 2``, ``bc 0`` and ``bc 1`` are created but ``bc 2`` does not exist. 

    Likewise, if ``number = 2`` and there is no ``subsection bc 0`` explicitly stated, the boundary is still created, with ``noflux`` by default.

* ``id`` is the number associated with the boundary condition. By default, Lethe assumes that the id is equivalent to the number of the bc.

* ``type``: type of boundary condition being imposed. At the moment, choices are:
    * ``noflux`` (default) so that there is no heat transfer boundary condition,
    * ``temperature`` (Dirichlet BC), to impose a given temperature ``value`` at the boundary,
    * ``convection-radiation`` (Robin BC) for cooling/heating, depending on the environment temperature at the boundary ``Tinf``, with a given heat transfer coefficient ``h`` and emissivity of the boundary :math:`\mathbf{\epsilon}` following Newton's law of cooling (and heating) and Stefan-Boltzmann law of radiation. Note that the expressions for ``h``, ``Tinf`` and ``emissivity`` can be time-dependent, but the current implementation doesn't allow for space dependence (the expressions are evaluated at the origin).

.. math::
    \frac{ \partial T}{\partial \mathbf{n}} = h (T - T_{inf}) + \epsilon \sigma (T^4 - T_{inf}^4)


where :math:`\mathbf{\sigma}` is the Stefan-Boltzmann constant.

.. seealso::

  The :doc:`../../examples/multiphysics/warming-up-a-viscous-fluid/warming-up-a-viscous-fluid` example uses heat transfer boundary conditions.


Tracer
^^^^^^

For tracer boundary conditions, the defaults parameters are:

.. code-block:: text

  subsection boundary conditions tracer
    set number         = 1
    set time dependent = false
    subsection bc 0
      set id   = 0
      set type = dirichlet
      subsection dirichlet
        set Function expression = 0
      end
    end
  end

* ``number``: This is the number of boundary conditions of the problem. 

* ``time dependent`` specifies if a boundary condition is time-dependent (``true``) or steady (``false``). By default, this parameter is set to ``false``. This improves the computational efficiency for transient cases in which the boundary conditions do not change.

* ``id`` is the number associated with the boundary condition. By default, Lethe assumes that the id is equivalent to the number of the bc.

* ``type``: This is the type of boundary condition being imposed. At the moment, only dirichlet boundary conditions can be imposed for tracer.


VOF
^^^

For VOF boundary conditions (multiphase flow), the possible ``types`` are ``none`` (default) and ``dirichlet``, as shown below.

.. code-block:: text

  subsection boundary conditions VOF
    set number         = 2
    set time dependent = false
    subsection bc 0
      set id   = 0
      set type = none
    end
    subsection bc 1
      set id   = 1
      set type = dirichlet
      subsection dirichlet
        set Function expression = 0
      end
    end
  end

.. warning::
    The ``number`` of boundary conditions must be specified explicitly. This is often a source of error.

.. note::
    The index in ``subsection bc ..`` must be coherent with the ``number`` of boundary conditions set: if ``number = 2``, ``bc 0`` and ``bc 1`` are created but ``bc 2`` does not exist. 

    Likewise, if ``number = 2`` and there is no ``subsection bc 0`` explicitly stated, the boundary is still created, with ``none`` by default.

* ``number``: This is the number of boundary conditions of the problem.

* ``time dependent`` specifies if a boundary condition is time-dependent (``true``) or steady (``false``). By default, this parameter is set to ``false``. This improves the computational efficiency for transient cases in which the boundary conditions do not change.

* ``id`` is the number associated with the boundary condition. By default, Lethe assumes that the id is equivalent to the number of the bc.

* ``type``: This is the type of boundary condition being imposed. At the moment, choices are:
    * ``none`` for which nothing happens.
    * ``dirichlet`` for inlet and outlet boundary conditions, to specify which fluid should be at the selected boundary.

.. note::
    For periodic boundary conditions, there is no need to specify anything in the ``boundary conditions VOF`` subsection. The periodic boundary condition must be specified in the ``boundary conditions`` subsection (see :doc:`boundary_conditions_cfd`).
    
    
Cahn-Hilliard
^^^^^^^^^^^^^^

For Cahn-Hilliard boundary conditions, the available ``types`` are ``noflux`` (default), ``dirichlet``, ``angle_of_contact``, and ``free_angle``. The parameters for each type of Cahn-Hilliard boundary conditions are:

.. code-block:: text

  subsection boundary conditions cahn hilliard
    set number         = 3
    set time dependent = false
    subsection bc 0
        set id            = 0
        set type          = dirichlet
        subsection phi
            set Function expression = 0
        end
     end
    subsection bc 1
      set id              = 1
      set type            = angle_of_contact
      set angle value     = 90 # The angle is given in degrees (°)
    end
    subsection bc 2
      set id              = 2
      set type            = free_angle
    end
  end
    

* ``number``: This is the number of boundary conditions of the problem. 

* ``time dependent`` specifies if a boundary condition is time-dependent (``true``) or not (``false``). By default, this parameter is set to ``false``. It is used to improve the computational efficiency of transient cases in which the boundary conditions do not change.

* ``id`` is the number associated with the boundary condition. By default, Lethe assumes that the id is equivalent to the number of the bc.

* ``type``: Type of boundary condition being imposed. At the moment, the choices are:
    * ``noflux`` (default): no phase leaves the simulation domain.
    * ``dirichlet``: Imposes a given phase order parameter function on the boundary. This function can depend on position (:math:`x,y,z`) and on time (:math:`t`).
    * ``angle_of_contact``: Imposes a given angle of contact ``angle value`` between the two phases at the boundary. It refers to the inner angle of contact, in degrees (°).
    * ``free_angle``: Leaves the angle as a free variable to be solved.

