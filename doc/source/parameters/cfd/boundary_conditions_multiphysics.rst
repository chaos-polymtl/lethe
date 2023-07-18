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
    set number = 2
    subsection bc 0
      set id    = 0
      set type  = temperature
      set value = 0
    end
    subsection bc 1
      set id         = 1
      set type       = convection-radiation
      set h          = 0
      set Tinf       = 0
      set emissivity = 0
    end
    set Stefan-Boltzmann constant = 0.000000056703
  end

* ``number``: This is the number of boundary conditions of the problem. 

.. warning::
    The ``number`` of boundary conditions must be specified explicitly. This is often a source of error.

.. note::
    The index in ``subsection bc ..`` must be coherent with the ``number`` of boundary conditions set: if ``number = 2``, ``bc 0`` and ``bc 1`` are created but ``bc 2`` does not exist. 

    Likewise, if ``number = 2`` and there is no ``subsection bc 0`` explicitly stated, the boundary is still created, with ``noflux`` by default.

* ``id`` is the number associated with the boundary condition. By default, Lethe assumes that the id is equivalent to the number of the bc.

* ``type``: type of boundary condition being imposed. At the moment, choices are:
    * ``noflux`` (default) so that there is no heat transfer boundary condition,
    * ``temperature`` (Dirichlet BC), to impose a given temperature ``value`` at the boundary,
    * ``convection-radiation`` (Robin BC) for cooling/heating, depending on the environment temperature at the boundary ``Tinf``, with a given heat transfer coefficient ``h`` and emissivity of the boundary :math:`\mathbf{\epsilon}` following Newton's law of cooling (and heating) and Stefan-Boltzmann law of radiation.

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
    set number = 1
    subsection bc 0
      set id   = 0
      set type = dirichlet
      subsection dirichlet
        set Function expression = 0
      end
    end
  end

* ``number``: This is the number of boundary conditions of the problem. 

* ``id`` is the number associated with the boundary condition. By default, Lethe assumes that the id is equivalent to the number of the bc.

* ``type``: This is the type of boundary condition being imposed. At the moment, only dirichlet boundary conditions can be imposed for tracer.


VOF
^^^

For VOF boundary conditions (multiphase flow), the possible ``types`` are ``none`` (default), ``peeling/wetting``, and ``dirichlet``, as shown below.

.. code-block:: text

  subsection boundary conditions VOF
    set number = 3
    subsection bc 0
      set id   = 0
      set type = none
    end
    subsection bc 1
      set id   = 1
      set type = peeling/wetting
    end
    subsection bc 2
      set id   = 2
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

* ``id`` is the number associated with the boundary condition. By default, Lethe assumes that the id is equivalent to the number of the bc.

* ``type``: This is the type of boundary condition being imposed. At the moment, choices are:
    * ``none`` for which nothing happens.
    * ``peeling/wetting`` for the fluid can attach to (`wet`) or detach from (`peel`) the boundary. The parameters for peeling/wetting are defined in the :doc:`volume_of_fluid` subsection of the parameter file.
    * ``dirichlet`` for inlet and outlet boundary conditions, to specify which fluid should be at the selected boundary.

.. note::
    For periodic boundary conditions, there is no need to specify anything in the ``boundary conditions VOF`` subsection. The periodic boundary condition must be specified in the ``boundary conditions`` subsection (see :doc:`boundary_conditions_cfd`).
    
    
Cahn-Hilliard
^^^^^^^^^^^^^^

For Cahn-Hilliard boundary conditions, the possible ``types`` are ``noflux`` (default), ``dirichlet`` and ``angle_of_contact``. The default parameters for ``dirichlet`` and ``angle_of_contact`` are shown:

.. code-block:: text

    subsection boundary conditions cahn hilliard
    set number                  = 2
    set time dependent          = false
        subsection bc 0
            set id 		= 0
            set type            = dirichlet
            subsection phi
                set Function expression = 0
            end
         subsection bc 1 
            set id              = 1
            set type            = angle_of_contact
            set angle value     = 90 # The angle is given in degrees (°) 
         end
    end
    

* ``number``: This is the number of boundary conditions of the problem. 

* ``time dependent`` specifies if a  boundary condition is time dependent (``true``) or steady (``false```). By default, this parameter is set to ``false``. This is there to improve the computational efficiency for transient cases in which the boundary conditions do not change. 

* ``id`` is the number associated with the boundary condition. By default, Lethe assumes that the id is equivalent to the number of the bc.

* ``type``: This is the type of boundary condition being imposed. At the moment, choices are:
    * ``noflux`` (default) so that no phase leave the simulation domain.
    * ``dirichlet`` to impose a given phase order parameter function on the boundary. This function can depend on position (:math:`x,y,z`) and on time (:math:`t`).
    * ``angle_of_contact`` to impose a given angle of contact ``angle value`` between the two phases at the boundary. It refers to the inner angle of contact, in degrees (°).

