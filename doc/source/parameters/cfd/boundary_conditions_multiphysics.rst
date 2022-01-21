==================================
Multiphysics - Boundary Contitions
==================================

This subsection's purpose is defining the boundary conditions associated to multiphysic problems. 

Heat Transfer
^^^^^^^^^^^^^

For heat transfer boundary conditions, the possible ``types`` are ``temperature`` (default) and ``convection``.
The default parameters of each are shown: 

.. code-block:: text

    subsection boundary conditions heat transfer
    set number                  = 2
        subsection bc 0
            set type	        = temperature
            set value	        = 0
        end
        subsection bc 1
            set type		    = convection
            set h 		    = 0
            set Tinf	   	    = 0
        end
    end

* ``number``: This is the number of boundary conditions of the problem. 

.. warning::
    The number of boundary conditions must be specified explicitly as the ParameterHandler is not able to deduce the number of boundary conditions from the number of ``bc`` subsections. This is often a source of error.

* ``type``: This is the type of boundary condition been imposed. At the moment, choices are
    * ``temperature`` (Dirichlet BC), to impose a given temperature ``value`` at the boundary 
    * ``convection`` (Robin BC) for cooling/heating, depending on the environment temperature at the boundary ``Tinf``, with a given heat transfer coefficient ``h`` following Newton's law of cooling (and heating)

.. math::
    \frac{ \partial T}{\partial \mathbf{n}} = h (T - T_\textit{inf})


.. note::
    By default, the boundary condition type is ``temperature`` with ``value = 0``

Tracer
^^^^^^

For tracer boundary conditions, the defaults parameters are:

.. code-block:: text

    subsection boundary conditions tracer
    set number                  = 1
        subsection bc 0
                set type              = dirichlet
                subsection dirichlet
                    set Function expression = 0
                end
        end
    end

* ``number``: This is the number of boundary conditions of the problem. 

* ``type``: This is the type of boundary condition been imposed. At the moment, only dirichlet boundary conditions can be imposed for tracer.

