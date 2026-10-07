=============
Stabilization
=============

To solve the Navier-Stokes equations (and other), Lethe uses stabilization techniques to formulate a Petrov-Galerkin strategy in which the test function is not strictly equal to the interpolation. The stabilization provided by Lethe are relatively robust and do not require any manual tinkering. However, stabilization of FEM schemes remain an active area a research, notably when it comes to variational multi-scales (VMS) methods. This is a field in which we are doing active research. As such, Lethe possesses some advanced parameters to control the stabilization techniques used to solve the Navier-Stokes. These are *advanced* parameters and, in general, the defaults value should be used.


.. code-block:: text

  subsection stabilization
    set use default stabilization        = true

    set stabilization                    = pspg_supg     # <pspg_supg|gls|grad_div>.

    # DCDD stabilization
    set heat transfer dcdd stabilization = false
    set cls dcdd stabilization           = true
    set cls dcdd diffusion factor        = 0.5

    # Pressure scaling factor
    set pressure scaling factor          = 1
    
    # Scalar limiter
    set scalar limiter                   = none #<none|moe|kuzmin>
  end
  

The ``use default stabilization`` indicates that the solver should use the default stabilization strategy, which is generally the most adequate strategy. To use an alternative strategy, this parameter must be set to false and the strategy to be used must be manually specified using the ``stabilization`` parameter.

There are three choices of stabilization strategy:

* ``stabilization=pspg-supg`` assembles a PSPG/SUPG stabilization for the Navier-Stokes equations. This stabilization should only be used with the monolithic solver for the Navier-Stokes equations (``lethe-fluid`` or ``lethe-fluid-matrix-free``).

* ``stabilization=gls`` assembles a full GLS stabilization for the Navier-Stokes equations which adds two Least-Squares terms (for more details see :doc:`../../../../theory/multiphysics/fluid_dynamics/stabilization`). This stabilization should only be used with the monolithic solver for the Navier-Stokes equations (``lethe-fluid`` or ``lethe-fluid-matrix-free``).

* ``stabilization=grad_div`` assembles a grad-div penalization term in the momentum equation to ensure mass conservation. This is not a stabilization method per-say and should not be used with elements that are not LBB stable. This stabilization should only be used with the grad-div block Navier-Stokes solver (``lethe-fluid-block``).

* ``heat transfer dcdd stabilization`` applies the Discontinuity-Capturing Directional Dissipation (DCDD) stabilization term on the heat transfer equation. For more information, see `Tezduyar, T. E. (2003) <https://doi.org/10.1002/fld.505>`_\.

* ``cls dcdd stabilization`` applies the DCDD stabilization term on the :doc:`CLS equation<../../theory/multiphase/cfd/cls>`. For more information, see `Tezduyar, T. E. (2003) <https://doi.org/10.1002/fld.505>`_\.

* ``cls dcdd diffusion factor`` is the diffusion coefficient applied to the DCDD stabilization term in the :doc:`CLS equation<../../theory/multiphase/cfd/cls>`.

* ``pressure scaling factor`` used as a multiplier for the pressure in the momentum equation; the inverse of the factor is applied to the pressure after solving. It helps the convergence of the linear solver by decreasing the condition number for cases where pressure and velocity have very different scales.

* ``scalar limiter`` applies a scalar limiter to the solution of the tracer equation when a Discontinuous Galerkin method (DG) is used. This is useful to prevent oscillations in the solution of the tracer equation, especially in cases with sharp gradients and when the diffusion coefficient is very small or zero. The limiters are applied after each time step and they do not change the average of the solution within a cell, hence they conserve the tracer. The available options are:

  * ``none`` disables the limiter.

  * ``moe`` applies a monotone upstream-centered scheme for conservation laws (MUSCL) type of limiter (see `Moe et al. (2015) <https://doi.org/10.48550/arXiv.1507.03024>`_). The deviation of the solution from its cell average is rescaled so that the solution remains between the minimum and the maximum of the solution in the neighboring cells. This limiter also reduces the smooth extrema of the solution.

  * ``kuzmin`` applies the hierarchical vertex-based limiter of `Kuzmin (2010) <https://doi.org/10.1016/j.cam.2009.05.028>`_. The linear part of the solution and its higher-order part are rescaled separately, using bounds that are established at the vertices of the cell from the cells that share them. The linear part is never limited more than the higher-order part, which preserves the smooth extrema of the solution. This property requires ``tracer degree`` to be at least 2. With a degree of 1, this limiter behaves as a classical slope limiter and it reduces the smooth extrema. This limiter does not have any parameter. The vertices located on the boundaries of the domain are not used to establish the bounds. Consequently, the solution is not limited in the cells of which all the vertices are on a boundary, which is the case of every cell of a mesh that is a single cell thick.

.. warning::

  The limiters only modify the solution within a cell and they preserve its average. Consequently, they cannot remove the oscillations that are generated by the time integration scheme when the time step is large, since these oscillations are present in the cell averages. With the ``bdf2`` scheme, a CFL number of the order of 0.25 or below is recommended when a limiter is used.
