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
    set vof dcdd stabilization           = true
    set vof dcdd diffusion factor        = 0.5

    # Pressure scaling factor
    set pressure scaling factor          = 1
    
    # Scalar limiter
    set scalar limiter                   = none #<none|moe>
  end
  

The ``use default stabilization`` indicates that the solver should use the default stabilization strategy, which is generally the most adequate strategy. To use an alternative strategy, this parameter must be set to false and the strategy to be used must be manually specified using the ``stabilization`` parameter.

There are three choices of stabilization strategy:

* ``stabilization=pspg-supg`` assembles a PSPG/SUPG stabilization for the Navier-Stokes equations. This stabilization should only be used with the monolithic solver for the Navier-Stokes equations (``lethe-fluid`` or ``lethe-fluid-matrix-free``).

* ``stabilization=gls`` assembles a full GLS stabilization for the Navier-Stokes equations which adds two Least-Squares terms (for more details see :doc:`../../../../theory/multiphysics/fluid_dynamics/stabilization`). This stabilization should only be used with the monolithic solver for the Navier-Stokes equations (``lethe-fluid`` or ``lethe-fluid-matrix-free``).

* ``stabilization=grad_div`` assembles a grad-div penalization term in the momentum equation to ensure mass conservation. This is not a stabilization method per-say and should not be used with elements that are not LBB stable. This stabilization should only be used with the grad-div block Navier-Stokes solver (``lethe-fluid-block``).

* ``heat transfer dcdd stabilization`` applies the Discontinuity-Capturing Directional Dissipation (DCDD) stabilization term on the heat transfer equation. For more information, see `Tezduyar, T. E. (2003) <https://doi.org/10.1002/fld.505>`_\.

* ``vof dcdd stabilization`` applies the DCDD stabilization term on the :doc:`VOF equation<../../theory/multiphase/cfd/vof>`. For more information, see `Tezduyar, T. E. (2003) <https://doi.org/10.1002/fld.505>`_\.

* ``vof dcdd diffusion factor`` is the diffusion coefficient applied to the DCDD stabilization term in the :doc:`VOF equation<../../theory/multiphase/cfd/vof>`.

* ``pressure scaling factor`` used as a multiplier for the pressure in the momentum equation; the inverse of the factor is applied to the pressure after solving. It helps the convergence of the linear solver by decreasing the condition number for cases where pressure and velocity have very different scales.

* ``scalar limiter`` applies a scalar limiter to the solution of the tracer equation when a Discontinuous Galerkin method (DG) is used. The available option are ``none`` and  ``moe``.  The ``none`` option disables the limiter. The ``moe`` limiter applies a monotone upstream-centered scheme for conservation laws (MUSCL) type of limiter (see `Moe et al. (2015) <https://doi.org/10.48550/arXiv.1507.03024>`). This is useful to prevent oscillations in the solution of the tracer equation, especially in cases with sharp gradients and when the diffusion coefficient is very small or zero.