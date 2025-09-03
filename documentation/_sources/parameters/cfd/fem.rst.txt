==================================
Finite Element Interpolation (FEM)
==================================

This subsection specifies the characteristics of the finite element method used in the simulations. The interpolation orders of velocity, pressure and other physics are specified. An additional option is also present to enable the use of high-order mapping throughout the entire mesh. At the moment, all variables are interpolated using Lagrange elements. On quad/hex meshes, `FE_Q <https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html>`_ elements are used whereas on simplex meshes (triangles, tetrahedron), `FE_P <https://www.dealii.org/current/doxygen/deal.II/classFE__SimplexP.html>`_ elements are used.


.. code-block:: text

  subsection FEM
    # interpolation order pressure
    set pressure order     = 1

    # interpolation order velocity
    set velocity order     = 1

    # interpolation order temperature
    set temperature order  = 1

    # interpolation order tracer
    set tracer order       = 1
    set tracer uses dg     = false

    # interpolation order vof
    set VOF order          = 1

    # interpolation order cahn hilliard
    set phase cahn hilliard order     = 1
    set potential cahn hilliard order = 1

    # bubble enrichment function
    set enable bubble function velocity = false
    set enable bubble function pressure = false
  end


* ``velocity order`` specifies the interpolation order for velocity.

* ``pressure order`` specifies the interpolation order for pressure. For the **lethe-fluid** family of solvers, the interpolation order for pressure can be equal or lower than that of velocity. For the **lethe-fluid-block** family of solvers, the interpolation order for pressure has to be one degree lower than that of velocity.

* ``temperature order`` specifies the interpolation order for temperature for the heat transfer physics.

* ``tracer order`` specifies the interpolation order for the tracer physics.

* ``tracer uses dg`` specifies if the Discontinuous Galerkin (DG) formulation is used instead of the Continuous Galerkin (CG) that is used by default, for the tracer physics. This formulation allows discontinuities between elements and at boundaries, and only penalizes them; this prevents oscillations and enforces local conservation of mass, at the cost of more numerous degrees of freedom.

* ``VOF order`` specifies the interpolation for the VOF phase indicator. It is not recommended to use higher order interpolation for the VOF method as this may conflict with the bounding and the sharpening mechanism used therein.

* ``phase cahn hilliard order`` and ``potential cahn hilliard order`` specify the interpolation order for the phase order parameter and the chemical potential in the Cahn-Hilliard equations. The orders chosen should be equal. They are left as two separate parameters for debugging purposes.

* ``enable bubble function velocity`` and ``enable bubble function pressure`` specifies if the bubble enrichment function is used in the velocity and pressure fields, respectively. This is a polynomial enrichment function centered at the mid-point of the cell and that vanishes at the element boundary. It can be used to improve accuracy and stability in Galerkin FEM, similarly to SUPG stabilization; we refer the reader to the work of  `Franca and Farhat 1995 <https://www.sciencedirect.com/science/article/abs/pii/004578259400721X>`_ and `Brezzi et al 1992 <https://www.sciencedirect.com/science/article/abs/pii/004578259290102P>`_ for more detail.

.. warning::
  
  The bubble enrichment function `FE_Q_Bubbles <https://www.dealii.org/current/doxygen/deal.II/classFE__Q__Bubbles.html>`_ is not compatible with simplex meshes.