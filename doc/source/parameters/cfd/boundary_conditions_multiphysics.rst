==================================
Boundary Contitions - Multiphysics
==================================

This subsection's purpose is defining the boundary conditions associated to multiphysic problems. 

Heat Transfer
^^^^^^^^^^^^^

For heat transfer boundary conditions, the possible ``types`` are ``temperature`` (default) and ``convection-radiation``.
The default parameters of each are shown: 

.. code-block:: text

    subsection boundary conditions heat transfer
    set number                  = 2
        subsection bc 0
	    set id 		= 0
            set type	        = temperature
            set value	        = 0
        end
        subsection bc 1
	    set id 		= 1
            set type		= convection-radiation
            set h 		= 0
            set Tinf 		= 0
            set emissivity  	= 0
        end
        set Stefan-Boltzmann constant = 0.000000056703
    end

* ``number``: This is the number of boundary conditions of the problem. 

.. warning::
    The number of boundary conditions must be specified explicitly as the ParameterHandler is not able to deduce the number of boundary conditions from the number of ``bc`` subsections. This is often a source of error.

* ``id`` is the number associated with the boundary condition. By default, Lethe assumes that the id is equivalent to the number of the bc.

* ``type``: This is the type of boundary condition been imposed. At the moment, choices are:
    * ``temperature`` (Dirichlet BC), to impose a given temperature ``value`` at the boundary.
    * ``convection-radiation`` (Robin BC) for cooling/heating, depending on the environment temperature at the boundary ``Tinf``, with a given heat transfer coefficient ``h`` and emissivity of the boundary :math:`\mathbf{\epsilon}` following Newton's law of cooling (and heating) and Stefan-Boltzmann law of radiation.

.. math::
    \frac{ \partial T}{\partial \mathbf{n}} = h (T - T_{inf}) + \epsilon \sigma (T^4 - T_{inf}^4)


where :math:`\mathbf{\sigma}` is the Stefan-Boltzmann constant.

.. note::
    By default, the boundary condition type is ``temperature`` with ``value = 0``

Tracer
^^^^^^

For tracer boundary conditions, the defaults parameters are:

.. code-block:: text

    subsection boundary conditions tracer
    set number                  = 1
        subsection bc 0
	        set id 		= 0
                set type        = dirichlet
                subsection dirichlet
                    set Function expression = 0
                end
        end
    end

* ``number``: This is the number of boundary conditions of the problem. 

* ``id`` is the number associated with the boundary condition. By default, Lethe assumes that the id is equivalent to the number of the bc.

* ``type``: This is the type of boundary condition been imposed. At the moment, only dirichlet boundary conditions can be imposed for tracer.


VOF
^^^

For VOF boundary conditions (multiphase flow), the possible ``types`` are ``none`` (default) and ``peeling/wetting``, as shown below.

.. code-block:: text

    subsection boundary conditions VOF
    set number                  = 2
        subsection bc 0
	    set id 		= 0
            set type	        = none
        end
        subsection bc 1
	    set id 		= 1
            set type		= peeling/wetting
        end
    end


* ``number``: This is the number of boundary conditions of the problem. 

* ``id`` is the number associated with the boundary condition. By default, Lethe assumes that the id is equivalent to the number of the bc.

* ``type``: This is the type of boundary condition been imposed. At the moment, choices are:
    * ``none`` for which nothing happens.
    * ``peeling/wetting`` for the fluid can attach to (`wet`) or detach from (`peel`) the boundary. The parameters for peeling/wetting are defined in the :doc:`volume_of_fluid` subsection of the parameter file.


