Boundary conditions - CFD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This subsection's defines the boundary conditions associated with fluid dynamics physics. Lethe supports slip, no-slip, periodic, and dirichlet boundary conditions which are set from arbitrary functions. These functions can be used to define all sort of steady-state and transient velocity boundary conditions such as rotating walls.

.. code-block:: text

   subsection boundary conditions
     set number                  = 2
     set time dependent          = false
     subsection bc 0
           set type              = function
          
           subsection u
               set Function expression = -y
           end
           subsection v
               set Function expression = x
           end
           subsection w
               set Function expression = 0
           end
   
           # Center of rotation used for torque calculation
           subsection center of rotation
              set x = 0
              set y = 0
              set z = 0
           end
     end
     subsection bc 1
           set type              = noslip
     end
   end

* ``number`` specifies the number of boundary conditions of the problem. Periodicity between 2 boundaries counts as 1 condition even if it requires two distinct boundary ids.

.. warning::
    The number of boundary conditions must be specified explicitely as the ParameterHandler is not able to deduce the number of boundary conditions from the number of ``bc`` subsections. This is often a source of error.

* ``time dependent`` specifies if a  boundary condition is time dependent (``true``) or steady (``false```). This parameter is only there to improve the computational efficient for transient cases in which the boundary conditions do not change.

* Each fluid dynamics boundary condition is stored in a ``bc no`` subsection :
    * ``type`` is the type of the boundary condition. Onlyslip, no-slip, periodic and function are enabled.
    
    * The subsections ``u``, ``v`` and ``w`` are used to specify the individual components of the velocity at the boundary using function expressions. These function can depend on position (:math:`x,y,z`) and on time (:math:`t`).

    * The ``center of rotation`` subsection is only necessary when calculating the torque applied on a boundary. See  See :doc:`force_and_torque` for more information.
