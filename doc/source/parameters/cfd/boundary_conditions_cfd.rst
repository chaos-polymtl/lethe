Boundary conditions - CFD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This subsection defines the boundary conditions associated with fluid dynamics physics. Lethe supports the following boundary conditions:

* ``noslip`` boundary conditions strongly impose the velocity on a boundary to be :math:`\mathbf{u}=[0,0]^T` and :math:`\mathbf{u}=[0,0,0]^T` in 2D and 3D respectively.
* ``slip`` boundary conditions impose :math:`\mathbf{u} \cdot \mathbf{n}=0`, with :math:`\mathbf{n}` the normal vector of the boundary. Imposing slip boundary conditions strongly is not trivial in FEM. We refer the reader to the deal.II `documentation <https://www.dealii.org/current/doxygen/deal.II/group__constraints.html>`_ for explanations on how this is achieved.
* ``periodic`` boundary conditions, in which fluid exiting the domain will reenter on the opposite side. 
* ``function`` where a Dirichlet boundary condition is set from an arbitrary function. These functions can be used to define all sorts of steady-state and transient velocity boundary conditions such as rotating walls.
* Finally, Lethe also supports not imposing a boundary condition on an ID. Not imposing a boundary condition is equivalent to the *do nothing* boundary condition, which results in a zero net traction on a boundary. This, in fact, imposes :math:`\int_{\Gamma}(-p\mathcal{I} + \mathbf{\tau}) \cdot \mathbf{n}=0` where :math:`p` is the pressure, :math:`\mathcal{I}` is the identity tensor, :math:`\mathbf{\tau}` is the deviatoric stress tensor  and :math:`\Gamma` is the boundary. 


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
        set periodic_id         = 1
        set periodic_direction  = 0
     end
     subsection bc 1
           set type              = noslip
     end
   end

* ``number`` specifies the number of boundary conditions of the problem. Periodicity between 2 boundaries counts as 1 condition even if it requires two distinct boundary ids.

.. warning::
    The number of boundary conditions must be specified explicitly as the ParameterHandler is not able to deduce the number of boundary conditions from the number of ``bc`` subsections. This is often a source of error.

* ``time dependent`` specifies if a  boundary condition is time dependent (``true``) or steady (``false```). By default, this parameter is set to ``false``. This is there to improve the computational efficiency for transient cases in which the boundary conditions do not change.

* Each fluid dynamics boundary condition is stored in a ``bc no`` subsection :
    * ``type`` is the type of the boundary condition. Only slip, no-slip, periodic and function are enabled.
    
    * The subsections ``u``, ``v`` and ``w`` are used to specify the individual components of the velocity at the boundary using function expressions. These functions can depend on position (:math:`x,y,z`) and on time (:math:`t`).

    * The ``center of rotation`` subsection is only necessary when calculating the torque applied on a boundary. See  See :doc:`force_and_torque` for more information.
    * ``periodic id`` and ``periodic_direction`` specify the id and direction of the matching periodic boundary condition. For example, if boundary id 0 (located at xmin) is matched with boundary id 1 (located at xmax), we would set ``ìd=0``, ``periodic_id=1`` and ``periodic_direction=0``.
