=============================================================
Finite Element Method for Time-Harmonic Maxwell's Equations
=============================================================

This section presents two finite-element formulations of the time-harmonic
Maxwell equations. We first introduce the conventional weak form (the primal
formulation), then the ultraweak formulation used by Lethe's DPG solver.

.. seealso::
    For additional details on the DPG setting, see :doc:`DPG formulation of time-harmonic Maxwell problems <./dpg_time_harmonic_maxwell>`.

Primal Formulation
------------------

Starting from the electric-field wave equation shown in
:doc:`Time Harmonic Maxwell's Equations <time_harmonic>`:

.. math::
    \nabla \times \left( \mu_{\mathrm{em}}^{-1} \nabla \times \mathbf{E} \right) -\omega^2 \varepsilon_{\mathrm{em,eff}} \mathbf{E} = i \omega \mathbf{J}_{\mathrm{ext}},

we consider a domain :math:`\Omega` with boundary :math:`\Gamma`.

.. note::
    In general, it is not necessary to solve both the electric- and
    magnetic-field wave equations, since one can be derived from the other.
    The field to solve for can be chosen based on the problem at hand.

Boundary conditions depend on the problem (for example Sommerfeld,
absorbing/radiation, or impedance conditions), but they can be grouped into
three main categories:

- First-type (Dirichlet) boundary conditions: :math:`\mathbf{n} \times \mathbf{E} = \mathbf{E}_{\mathrm{D}}`;
- Second-type (Neumann) boundary conditions: :math:`\mu_{\mathrm{em}}^{-1} \mathbf{n} \times (\nabla \times \mathbf{E})  = i \omega \mathbf{J}_{\mathrm{ext,N}}`;
- Robin boundary conditions: :math:`\mathbf{n} \times (\mu_{\mathrm{em}}^{-1} \nabla \times \mathbf{E}) + i \omega Z_{\mathrm{s}}^{-1} \mathbf{n} \times ( \mathbf{E} \times \mathbf{n}) = \mathbf{E}_\mathrm{R}`,

In the expressions above, :math:`\mathbf{n}` is the outward unit normal and
:math:`Y_\mathrm{s}` is the boundary surface admittance. The Robin condition
generalizes Dirichlet and Neumann conditions, and can represent impedance,
absorbing, and related boundary models.

For simplicity, we consider a perfect electric conductor (PEC), so
:math:`\mathbf{E}_{\mathrm{D}} = 0`. Multiplying the strong form by a complex
test function :math:`\mathbf{v}` satisfying
:math:`\mathbf{v} \times \mathbf{n}=0`, and integrating over :math:`\Omega`,
gives:

.. math::
    \begin{align*}
     \int_{\Omega} i \omega \mathbf{J}_{\mathrm{ext}} \cdot \mathbf{v^*} \mathrm{d}\Omega = &
     \int_{\Omega}  \mu_{\mathrm{em}}^{-1} (\nabla \times \mathbf{E}) \cdot (\nabla \times \mathbf{v^*}) \mathrm{d}\Omega - \int_{\Omega} \omega^2 \varepsilon_{\mathrm{em},\mathrm{eff}} \mathbf{E} \cdot \mathbf{v^*} \mathrm{d}\Omega \\ & + \int_{\Gamma} \mu_{\mathrm{em}}^{-1} (\nabla \times \mathbf{E}) \cdot (\mathbf{v^*} \times \mathbf{n}) \mathrm{d}\Gamma .
    \end{align*}


Because the fields are complex-valued, the usual :math:`L^2` inner product is
replaced by :math:`\int_{\Omega} \mathbf{u}\mathbf{v^*} \, \mathrm{d}\Omega`,
where :math:`\mathbf{v^*}` denotes the complex conjugate of
:math:`\mathbf{v}`. For these integrals to be well-defined,
:math:`\mathbf{E}` and :math:`\mathbf{v}` must belong to the Sobolev space
:math:`H(\mathrm{curl}, \Omega)`:

.. math::
    H(\mathrm{curl}, \Omega) = \{ \mathbf{v} \in [L^2(\Omega)]^3 : \nabla \times \mathbf{v} \in [L^2(\Omega)]^3 \},

and the Dirichlet boundary condition imposes the homogeneous curl space:

.. math::
    H_0(\mathrm{curl}, \Omega) = \{ \mathbf{v} \in H(\mathrm{curl}, \Omega) : \mathbf{n} \times \mathbf{v}|_{\Gamma} = 0 \}.
    
Thus, the boundary term vanishes
(:math:`\int_{\Gamma} \mu_{\mathrm{em}}^{-1} (\nabla \times \mathbf{E}) \cdot (\mathbf{v^*} \times \mathbf{n}) \mathrm{d}\Gamma = 0`),
and the weak form of the time-harmonic Maxwell equations becomes:

.. math::
    B(\mathbf{E}, \mathbf{v}) = &\int_{\Omega}  \mu_{\mathrm{em}}^{-1} (\nabla \times \mathbf{E}) \cdot (\nabla \times \mathbf{v^*}) \mathrm{d}\Omega - \int_{\Omega} \omega^2 \varepsilon_{\mathrm{em},\mathrm{eff}} \mathbf{E} \cdot \mathbf{v^*} \mathrm{d}\Omega  \\
    L(\mathbf{v}) = &\int_{\Omega} i \omega \mathbf{J}_{\mathrm{ext}} \cdot \mathbf{v^*} \mathrm{d}\Omega .

Formally, :math:`\mathbf{E}` should also satisfy Gauss's law
(:math:`\nabla \cdot \mathbf{D} = \rho_f`). In this setting, it is
implicitly enforced by the electromagnetic wave equation and is therefore
contained in the weak form above.

Ultraweak Formulation
---------------------

The formulation above is the conventional weak form, used here as a reference
for comparison with other implementations of the time-harmonic Maxwell
equations. In Lethe's DPG setting, the ultraweak formulation is used instead.
The main practical difference is that the primal formulation above eliminates
the magnetic field to obtain a single second-order equation for the electric
field, whereas the ultraweak formulation solves for the electric and magnetic
fields simultaneously by weakening both equations without eliminating either
one.

Retaining both fields as independent unknowns is especially useful when the
electromagnetic power dissipated in the medium must be evaluated, since some
materials exhibit non-negligible magnetic losses in addition to dielectric
losses; recovering :math:`\mathbf{H}` directly, rather than reconstructing it from
:math:`\mathbf{E}` after the fact, avoids an additional loss of accuracy in
that calculation. Additionally, studies have shown that the ultraweak
formulation is more robust than the primal formulation in the context of the
DPG method.

The starting point is the time-harmonic Maxwell system written in the form

.. math::

   \begin{align}
   \nabla \times \mathbf{E} - i\omega \mu_{r} \mathbf{H}       &= 0, \\
   \nabla \times \mathbf{H} + i\omega \varepsilon_{r,\mathrm{eff}} \mathbf{E} &= \mathbf{J},
   \end{align}

where :math:`\varepsilon_{r,\mathrm{eff}} = \varepsilon_r + i \frac{\sigma_e}{\omega \varepsilon_0}`. The system is supplemented with boundary conditions of three kinds:

.. math::

   \begin{align}
     \mathbf{n} \times \mathbf{E} &= \mathbf{n} \times \mathbf{E}_D; \\
     \mathbf{n} \times \mathbf{H} &= \mathbf{n} \times \mathbf{H}_N - \mathbf{J}_{\mathrm{s}}; \\
     \mathbf{n} \times \mathbf{H} + Z_{\mathrm{s}}^{-1}\mathbf{n} \times (\mathbf{E}\times \mathbf{n}) &= \mathbf{E}_\mathrm{R}.
   \end{align}

To obtain the ultraweak formulation, Faraday's law and Ampère-Maxwell equations are multiplied by independent test functions :math:`\mathbf{I}` and :math:`\mathbf{F}`,
respectively, and integrated by parts over the mesh :math:`\Omega_h`:

.. math::

   \begin{align}
      (\nabla \times\mathbf{I},\mathbf{E})_{\Omega_h} - ( \mathbf{I}, i\omega \mu_r\mathbf{H})_{\Omega_h} + \langle \mathbf{I}, \mathbf{n}\times \hat{\mathbf{E}} \rangle_{\partial \Omega_h}
        &= 0 , \\
      (\nabla \times \mathbf{F},\mathbf{H})_{\Omega_h} + (\mathbf{F}, i \omega \varepsilon_{r,\text{eff}} \mathbf{E})_{\Omega_h}
      + \langle \mathbf{F} , \mathbf{n}\times \hat{\mathbf{H}}\rangle_{\partial \Omega_h \backslash \Gamma_R }
      - \langle \mathbf{F} ,Y_\mathrm{s} \mathbf{E}\rangle_{ \Gamma_R }
        &=(\mathbf{F},\mathbf{J})_{\Omega_h} -  \langle \mathbf{F} , \mathbf{E_\mathrm{R}} \rangle_{ \Gamma_R },
   \end{align}

where the Robin boundary condition is applied to the electric field without
loss of generality, and :math:`\hat{\mathbf{E}}`, :math:`\hat{\mathbf{H}}` are
the trace unknowns on the mesh skeleton. Each product above is understood to
be complex-conjugate in its first argument:

.. math::

   \begin{align}
   (\mathbf{v},\mathbf{u})_{\Omega_h} &:= \sum_{K\in\Omega_h} \int_{K} \mathbf{v}^* \cdot \mathbf{u}\,\mathrm{d}x, &
   \langle \mathbf{v},\mathbf{u} \rangle_{\partial\Omega_h} &:= \sum_{K\in\Omega_h} \int_{\partial K} \mathbf{v}^* \cdot \mathbf{u}\, \mathrm{d}s,
   \end{align}
   
with :math:`\mathbf{v}` and :math:`\mathbf{u}` being arbitrary complex-valued vector fields. From the above ultraweak system, the trial and test spaces are defined as :

.. math::

   \begin{align}
     \mathbf{E},\mathbf{H} &\in \mathcal{X}_{\mathrm{em}},
       & \mathcal{X}_{\mathrm{em}} &:= (L^2(\Omega_h))^3, \\
     \hat{\mathbf{E}},\hat{\mathbf{H}} &\in \hat{\mathcal{X}}_{\mathrm{em}},
       & \hat{\mathcal{X}}_{\mathrm{em}} &:= \bigl\{ \hat{\mathbf{x}} \in H^{-1/2}(\mathrm{curl},\Omega_h) : \mathbf{n}\times \hat{\mathbf{x}} = \mathbf{n}\times \mathbf{x}_\mathrm{D} \text{ on } \Gamma_D \bigr\}, \\
     \mathbf{F},\mathbf{I} &\in \mathcal{Y}_{\mathrm{em}},
       & \mathcal{Y}_{\mathrm{em}} &:= \bigl\{ \mathbf{y} \in H(\mathrm{curl},\Omega_h) \bigr\}.
   \end{align}

For an ultraweak formulation, the Neumann boundary condition is, in fact, a
Dirichlet condition on the flux variable (here, the magnetic field
:math:`\mathbf{H}`). This is why only a Dirichlet-type condition needs to be
imposed on the trace space :math:`\hat{\mathcal{X}}_{\mathrm{em}}`, which
covers both the electric and magnetic trace unknowns.

To complete the discrete DPG formulation, the test space is equipped with the
adjoint graph norm induced by the operator associated with the bilinear form
of the system:

.. math::

   \begin{align}
   \|(\mathbf{F}, \mathbf{I})\|^2_V &= \|  \nabla \times \mathbf{F} - i \omega \mu_r \mathbf{I} \|^2_{\Omega_h} + \| \nabla \times \mathbf{I} + i \omega \varepsilon_{r,\mathrm{eff}} \mathbf{F} \|^2_{\Omega_h} \\
     &\quad + \alpha\left(\|\mathbf{F}\|^2_{\Omega_h} + \|\mathbf{I}\|^2_{\Omega_h}\right)
      + \| \mathbf{n} \times \mathbf{I} + Y_\mathrm{s} \mathbf{n} \times (\mathbf{F} \times \mathbf{n}) \|^2_{\Gamma_R},
   \end{align}

where :math:`\alpha \sim 1` is a strictly positive constant (:math:`\alpha=1`
in practice) introduced to ensure that the test norm remains localizable once
the test space is broken. The boundary term on
:math:`\Gamma_R` enforces the Robin boundary condition and is required to
obtain convergence of the method when this type of boundary condition is
present.


.. tip::

    The general Dirichlet, Neumann, and Robin conditions introduced above cover a
    range of standard boundary condition types encountered in electromagnetic
    simulations, each corresponding to a particular choice of :math:`\mathbf{E}_\mathrm{D}`,
    :math:`\mathbf{H}_\mathrm{N}`, :math:`Y_\mathrm{s}`, and :math:`\mathbf{E}_\mathrm{R}`. Writing
    :math:`\mathbf{E}_\parallel = \mathbf{n} \times (\mathbf{E} \times \mathbf{n})`
    for the tangential component of a field (and similarly
    :math:`\mathbf{H}_\parallel`), the list below is not exhaustive but covers the
    conditions commonly used in practice. 

    * **Perfect electric conductor (PEC)**: :math:`\mathbf{n} \times \mathbf{E} = 0`.
    * **Perfect magnetic conductor (PMC)**: :math:`\mathbf{n} \times \mathbf{H} = 0`.
    * **Prescribed electric field (Dirichlet)**: :math:`\mathbf{n} \times \mathbf{E} = \mathbf{n} \times \mathbf{E}_D`.
    * **Prescribed magnetic field (Neumann)**: :math:`\mathbf{n} \times \mathbf{H} = \mathbf{n} \times \mathbf{H}_N - \mathbf{J}_{\mathrm{s}}`, with :math:`\mathbf{J}_{\mathrm{s}}` a prescribed surface current density.
    * **Impedance condition (Robin)**: :math:`\mathbf{n} \times \mathbf{H} + Y_\mathrm{s} \mathbf{E}_\parallel = \mathbf{E_\mathrm{R}}`, with :math:`Y_\mathrm{s}` the surface admittance of the boundary.
    * **Absorbing / radiation condition (Silver–Müller)**: :math:`\mathbf{n} \times \mathbf{H} + \sqrt{\varepsilon_{r,\mathrm{eff}}/\mu_r}\, \mathbf{E}_\parallel = 0`. This is the impedance condition above with :math:`Y_\mathrm{s} = \sqrt{\varepsilon_{r,\mathrm{eff}}/\mu_r}` (the admittance of the exterior medium) and no source term.
    * **Lossy / imperfect conductor**: :math:`\mathbf{n} \times \mathbf{H} + \sqrt{\varepsilon_{r,\mathrm{eff},2}/\mu_{r,2}}\, \mathbf{E}_\parallel = 0`, the same form as the radiation condition above, using the (possibly complex) effective properties :math:`\varepsilon_{r,\mathrm{eff},2}`, :math:`\mu_{r,2}` of the adjacent conducting medium.
    * **Waveguide inlet port, TE**\ :math:`_{mn}` **mode**: :math:`\mathbf{n} \times \mathbf{H} + \dfrac{\mathbf{k} \cdot \mathbf{n}}{\omega \mu_r} \mathbf{E}_\parallel = \mathbf{n} \times \mathbf{H}_{inc} + \dfrac{\mathbf{k} \cdot \mathbf{n}}{\omega \mu_r} \mathbf{E}_{inc,\parallel}`.
    * **Waveguide inlet port, TM**\ :math:`_{mn}` **mode**: :math:`\mathbf{n} \times \mathbf{H} + \dfrac{\omega \varepsilon_{r,\mathrm{eff}}}{\mathbf{k} \cdot \mathbf{n}} \mathbf{E}_\parallel = \mathbf{n} \times \mathbf{H}_{inc} + \dfrac{\omega \varepsilon_{r,\mathrm{eff}}}{\mathbf{k} \cdot \mathbf{n}} \mathbf{E}_{inc,\parallel}`.

    In each waveguide, :math:`\mathbf{H}_{inc}` and
    :math:`\mathbf{E}_{inc,\parallel}` denote the incident field of the excited
    mode, and :math:`k = \omega \sqrt{\varepsilon_{r,\mathrm{eff}} \mu_r}` its wavevector amplitude.

.. note::
    When solving time-harmonic problems numerically, the mesh resolution should be chosen to ensure that the dimensionless wavenumber :math:`k h / (2 \pi)` is sufficiently small, where :math:`h` is the characteristic mesh size. A common rule of thumb is to have at least 10 degree 1 elements per wavelength, i.e., :math:`k h / (2 \pi) \leq 0.1`.