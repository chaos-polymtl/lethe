=========================
Boundary Conditions - CFD
=========================

This subsection defines the boundary conditions associated with fluid dynamics physics. Lethe supports the following boundary conditions:

* ``none`` boundary condition (default).
* ``noslip`` boundary conditions strongly impose the velocity on a boundary to be :math:`\mathbf{u}=[0,0]^T` and :math:`\mathbf{u}=[0,0,0]^T` in 2D and 3D respectively.
* ``slip`` boundary conditions impose :math:`\mathbf{u} \cdot \mathbf{n}=0`, with :math:`\mathbf{n}` the normal vector of the boundary. Imposing slip boundary conditions strongly is not trivial in FEM. We refer the reader to the deal.II `documentation <https://www.dealii.org/current/doxygen/deal.II/group__constraints.html>`_ for explanations on how this is achieved.
* ``partial slip`` boundary condition simulates an intermediary between ``slip`` and ``noslip`` boundary conditions, in which the fluid feels an attenuated stress due to the walls. The attenuation is controlled by the  boundary layer thickness (m). The ``partial slip`` boundary condition introduces the "penalization factor" ``beta`` to the :math:`\mathbf{n}` normal vector of the boundary, and the ``boundary layer thickness`` (m) as a parameter to calculate the shear stress at the boundaries.
* ``periodic`` boundary conditions, in which fluid exiting the domain will reenter on the opposite side. 
* ``function`` where a Dirichlet boundary condition is set from an arbitrary function. These functions can be used to define all sorts of steady-state and transient velocity boundary conditions such as rotating walls. It is also possible to weakly impose a Dirichlet boundary condition. In this case, the type should be set to ``function weak``. This will result in Nitsche method being used to weakly impose the boundary condition instead of it being strongly imposed by overwriting the values of the degrees of freedom. The ``function weak`` boundary type should only be used in very specific cases where the problem is very stiff, for example when it is fully enclosed and a non-trivial velocity profile is imposed. It can also be used to impose an outbound dirichlet boundary condition.
* ``outlet`` where a *do nothing* boundary condition, which is a zero traction, is imposed when the fluid is leaving the domain (:math:`\mathbf{u} \cdot \mathbf{n}>0`) and a penalization is imposed when the fluid is inbound. This is useful when turbulent structures or vortices are leaving the domain since it prevents the re-entry of the fluid. The boundary condition imposed is thus:

.. math::

  \nu \nabla \mathbf{u} \cdot \mathbf{n} - p \mathcal{I} \cdot \mathbf{n} - \beta (\mathbf{u}\cdot \mathbf{n})_{-} \mathbf{u} = 0

or in Einstein notation:

.. math::
    \nu \partial_i u_j n_j  - p n_i - \beta ( u_k n_k)_{-} u_i = 0

where :math:`\beta` is a constant  and :math:`(\mathbf{u}\cdot \mathbf{n})_{-}` is :math:`\min (0,\mathbf{u}\cdot \mathbf{n})`. We refer the reader to the work of `Arndt et al 2015 <https://www.mathsim.eu/~darndt/files/ENUMATH_2015.pdf>`_  for more detail.

* Finally, Lethe also supports not imposing a boundary condition on an ID. Not imposing a boundary condition is equivalent to the *do nothing* boundary condition (``none``), which results in a zero net traction on a boundary. This, in fact, imposes :math:`\int_{\Gamma}(-p\mathcal{I} + \mathbf{\tau}) \cdot \mathbf{n}=0` where :math:`p` is the pressure, :math:`\mathcal{I}` is the identity tensor, :math:`\mathbf{\tau}` is the deviatoric stress tensor  and :math:`\Gamma` is the boundary. 


.. code-block:: text

  subsection boundary conditions
    set number                = 2
    set time dependent        = false
    set fix pressure constant = false
    subsection bc 0
      set id                 = 0
      set type               = function

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

      set periodic_id        = 1
      set periodic_direction = 0
      set beta               = 1
    end
    subsection bc 1
      set type = noslip
    end
  end

* ``number`` specifies the number of boundary conditions of the problem. Periodicity between 2 boundaries counts as 1 condition even if it requires two distinct boundary ids.

.. warning::
    The ``number`` of boundary conditions must be specified explicitly. This is often a source of error.

.. note::
    The index in ``subsection bc ..`` must be coherent with the ``number`` of boundary conditions set: if ``number = 2``, ``bc 0`` and ``bc 1`` are created but ``bc 2`` does not exist. 

    Likewise, if ``number = 2`` and there is no ``subsection bc 0`` explicitly stated, the boundary is still created, with ``none`` by default.

* ``time dependent`` specifies if a boundary condition is time-dependent (``true``) or steady (``false``). By default, this parameter is set to ``false``. This is here to improve the computational efficiency for transient cases in which the boundary conditions do not change.

* ``fix pressure constant`` specifies if a zero pressure constraint should be applied on a single node of the coarse grid solver when using geometric multigrid preconditioning combined with the **lethe-fluid-matrix-free** solver. Essentially, this condition should be set to true whenever a user is using the **lethe-fluid-matrix-free** solver and simulating the flow within a closed domain (that is a domain on which all boundaries are either periodic or Dirichlet boundary conditions).

* Each fluid dynamics boundary condition is stored in a ``bc #`` subsection :
    * ``id``  is the number associated with the boundary condition. By default, Lethe assumes that the id is equivalent to the number ``#`` of the bc. 
    
    * ``type`` is the type of the boundary condition.
    
    * The subsections ``u``, ``v`` and ``w`` are used to specify the individual components of the velocity at the boundary using function expressions. These functions can depend on position (:math:`x,y,z`) and on time (:math:`t`).

    * The ``center of rotation`` subsection is only necessary when calculating the torque applied on a boundary. See  See :doc:`force_and_torque` for more information.

    * ``periodic id`` and ``periodic_direction`` specify the id and direction of the matching periodic boundary condition. For example, if boundary id 0 (located at xmin) is matched with boundary id 1 (located at xmax), we would set ``id = 0``, ``periodic_id = 1`` and ``periodic_direction = 0``.

    * ``beta`` is a penalization parameter used for both the ``outlet``, ``partial slip``, and ``function weak`` boundary conditions. For the outlet boundary conditions ``beta`` should be close to unity, whereas ``beta`` of 10 or a 100 can be appropriate for the ``function weak`` boundary condition. For the ``partial slip`` condition, use high values of ``beta`` (`i.e.` > 50).

    * ``boundary layer thickness`` (:math:`d_w`) is the parameter applied to the ``partial slip`` boundary condition. It is used to estimate the tangential shear stress :math:`\tau_t = -\mu \frac{u}{d_w}`. For very high ``boundary layer thicknes``, the boundary layer should behave exactly like the ``slip`` condition.

.. caution::
	While using the ``lethe-fluid-sharp`` solver, it is wise to assign a weak type of boundary (``outlet``, ``partial slip``, or ``function weak``) to at least one boundary. The presence of particle(s) has a non-null contribution to the divergence of the problem, making it much harder for the linear solver to converge unless it is given some flexibility through of boundaries.

.. caution::
  The ``lethe-fluid-matrix-free`` application does not support the ``pressure`` and ``partial slip`` boundary conditions. 
