=============
Stabilization
=============

To solve the Navier-Stokes equations (and other), Lethe uses stabilization techniques to formulate a Petrov-Galerkin strategy in which the test function is not strictly equal to the interpolation. The stabilization provided by Lethe are relatively robust and do not require any manual tinkering. However, stabilization of FEM schemes remain an active area a research, notably when it comes to variational multi-scales (VMS) methods. This is a field in which we are doing active research. As such, Lethe possesses some advanced parameters to control the stabilization techniques used to solve the Navier-Stokes. These are *advanced* parameters and, in general, the defaults value should be used.


.. code-block:: text

  subsection stabilization
    set use default stabilization        = true

    set stabilization                    = pspg_supg     # <pspg_supg|gls|grad_div>.

    # Stabilization parameter (tau) definition (lethe-fluid-matrix-free only)
    set stabilization parameter definition = tezduyar     # <tezduyar|metric_tensor>
    set vms c_i                            = 36
    set vms c_t                            = 4

    # DCDD stabilization
    set heat transfer dcdd stabilization = false
    set cls dcdd stabilization           = true
    set cls dcdd diffusion factor        = 0.5

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

* ``stabilization parameter definition`` controls how the stabilization parameter :math:`\tau` is computed (only supported by ``lethe-fluid-matrix-free``). The ``tezduyar`` option (default) uses the classical scalar definition based on a single isotropic element size :math:`h`, which is robust and accurate on relatively homogeneous meshes. The ``metric_tensor`` option uses the residual-based variational multiscale (VMS) definition built from the covariant element metric tensor :math:`\boldsymbol{G} = 4\, \boldsymbol{J}^{-T} \boldsymbol{J}^{-1}` (with :math:`\boldsymbol{J}` the mapping Jacobian), giving :math:`\tau = \left( C_t/\Delta t^2 + \boldsymbol{u} \cdot \boldsymbol{G} \boldsymbol{u} + C_I \nu^2\, \boldsymbol{G}:\boldsymbol{G} \right)^{-1/2}`. This definition is direction-aware and provides better conditioning on stretched, sheared or unstructured meshes (see `Bazilevs et al. (2007) <https://doi.org/10.1016/j.cma.2007.07.016>`_ and `Tezduyar and Osawa (2000) <https://doi.org/10.1016/S0045-7825(00)00211-5>`_).

* ``vms c_i`` is the constant :math:`C_I` scaling the viscous term of the ``metric_tensor`` stabilization parameter. It is degree-dependent; the default value of ``36`` targets linear elements.

* ``vms c_t`` is the constant :math:`C_t` scaling the transient term of the ``metric_tensor`` stabilization parameter.

* ``heat transfer dcdd stabilization`` applies the Discontinuity-Capturing Directional Dissipation (DCDD) stabilization term on the heat transfer equation. For more information, see `Tezduyar, T. E. (2003) <https://doi.org/10.1002/fld.505>`_\.

* ``cls dcdd stabilization`` applies the DCDD stabilization term on the :doc:`CLS equation<../../theory/multiphase/cfd/cls>`. For more information, see `Tezduyar, T. E. (2003) <https://doi.org/10.1002/fld.505>`_\.

* ``cls dcdd diffusion factor`` is the diffusion coefficient applied to the DCDD stabilization term in the :doc:`CLS equation<../../theory/multiphase/cfd/cls>`.

* ``pressure scaling factor`` used as a multiplier for the pressure in the momentum equation; the inverse of the factor is applied to the pressure after solving. It helps the convergence of the linear solver by decreasing the condition number for cases where pressure and velocity have very different scales.

* ``scalar limiter`` applies a scalar limiter to the solution of the tracer equation when a Discontinuous Galerkin method (DG) is used. The available option are ``none`` and  ``moe``.  The ``none`` option disables the limiter. The ``moe`` limiter applies a monotone upstream-centered scheme for conservation laws (MUSCL) type of limiter (see `Moe et al. (2015) <https://doi.org/10.48550/arXiv.1507.03024>`). This is useful to prevent oscillations in the solution of the tracer equation, especially in cases with sharp gradients and when the diffusion coefficient is very small or zero.