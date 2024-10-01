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

    # interpolation order vof
    set VOF order          = 1

    #interpolation order cahn hilliard
    set phase cahn hilliard order     = 1
    set potential cahn hilliard order = 1
  end


* ``velocity order`` specifies the interpolation order for velocity.

* ``pressure order`` specifies the interpolation order for pressure. For the **lethe-fluid** family of solvers, the interpolation order for pressure can be equal or lower than that of velocity. For the **lethe-fluid-block** family of solvers, the interpolation order for pressure has to be one degree lower than that of velocity.

* ``temperature order`` specifies the interpolation order for temperature for the heat transfer physics.

* ``tracer order`` specifies the interpolation order for the tracer physics.

* ``VOF order`` specifies the interpolation for the VOF phase indicator. It is not recommended to use higher order interpolation for the VOF method as this may conflict with the bounding and the sharpening mechanism used therein.

* ``phase cahn hilliard order`` and ``potential cahn hilliard order`` specify the interpolation order for the phase order parameter and the chemical potential in the Cahn-Hilliard equations. The orders chosen should be equal. They are left as two separate parameters for debugging purposes.


