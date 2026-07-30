=============================================================
Finite Element Method for Time-Harmonic Maxwell's Equations
=============================================================

This section describes the FEM formulation of the time-harmonic Maxwell's equations. Here, we will present the common weak form which is referred to as the primal formulation, and then we will present the ultra-weak formulation used by Lethe's DPG solver. For more details, go see the section on the :doc:`DPG formulation of time-harmonic Maxwell problems <./dpg_time_harmonic_maxwell>` section. 

Primal Formulation
------------------

Starting from the combined strong form of the combined electromagnetic field showed in the :doc:`Time Harmonic Maxwell's Equations <time_harmonic>` section:

.. math::
    \nabla \times \left( \mu_{\mathrm{em}}^{-1} \nabla \times \mathbf{E} \right) -\omega^2 \varepsilon_{\mathrm{em}_{eff}} \mathbf{E} = i \omega \mathbf{J}_{\mathrm{ext}}, \\

we consider a domain :math:`\Omega` with boundary :math:`\Gamma`.

.. note::
    In general, it is not necessary to solve both the magnetic and the electric field wave equation since one results from the other. One can choose the field to solve that suits the best the problem under consideration.

The choice of boundary conditions can vary greatly depending on the problem to be solved (e.g., Sommerfeld condition, absorbing boundary condition, impedance condition, etc.), but they can be classified into three main categories:

- First-type (Dirichlet) boundary conditions: :math:`\mathbf{n} \times \mathbf{E} = \mathbf{E}_{\mathrm{D}}`;
- Second-type (Neumann) boundary conditions: :math:`\mu_{\mathrm{em}}^{-1} \mathbf{n} \times (\nabla \times \mathbf{E})  = -i \omega \mathbf{J}_{\mathrm{ext,N}}`;
- Robin boundary conditions: :math:`\mathbf{n} \times (\mu_{\mathrm{em}}^{-1} \nabla \times \mathbf{E}) + i \omega Y_\mathrm{s} \mathbf{n} \times ( \mathbf{E} \times \mathbf{n}) = \mathbf{E}_\mathrm{R}`,

In the above, the :math:`\mathbf{n}` represent the unit normal vector, and :math:`Y_\mathrm{s}` is the surface admittance of the boundary. The Robin condition is a generalization of the Dirichlet and Neumann conditions, and can be used to model impedance boundaries, absorbing boundaries, and other types of boundary conditions.

For simplicity, in the current derivation of the weak form, a perfect electric conductor (PEC) is considered, which implies that :math:`\mathbf{E}_{\mathrm{D}} = 0`. Now, multiplying the strong form by a complex test function :math:`\mathbf{v}` that satisfies :math:`\mathbf{v} \times \mathbf{n}=0` and integrating over the domain :math:`\Omega`:

.. math::
    \begin{align*}
     \int_{\Omega} -i \omega \mathbf{J}_{\mathrm{ext}} \cdot \mathbf{v^*} \mathrm{d}\Omega = &
     \int_{\Omega}  \mu_{\mathrm{em}}^{-1} (\nabla \times \mathbf{E}) \cdot (\nabla \times \mathbf{v^*}) \mathrm{d}\Omega - \int_{\Omega} \omega^2 \varepsilon_{\mathrm{em}_{eff}} \mathbf{E} \cdot \mathbf{v^*} \mathrm{d}\Omega \\ & + \int_{\Gamma} \mu_{\mathrm{em}}^{-1} (\nabla \times \mathbf{E}) \cdot (\mathbf{v^*} \times \mathbf{n}) \mathrm{d}\Gamma .
    \end{align*} \\


Note that the terms above involve complex values, so the usual :math:`L^2` inner product is replaced by :math:`\int_{\Omega} \mathbf{u}\mathbf{v^*} \mathrm{d}\Omega`, where :math:`\mathbf{v^*}` is the complex conjugate of :math:`\mathbf{v}`. Finally, for the above integral to be well defined, :math:`\mathbf{E}` and :math:`\mathbf{v}` must belong to the Sobolev space :math:`H(curl, \Omega)` defined by:

.. math::
    H(curl, \Omega) = \{ \mathbf{v} \in [L^2(\Omega)]^3 : \nabla \times \mathbf{v} \in [L^2(\Omega)]^3 \},

and the Dirichlet boundary condition imposed the homogeneous curl space:

.. math::
    H_0(curl, \Omega) = \{ \mathbf{v} \in H(curl, \Omega) : \mathbf{n} \times \mathbf{v}|_{\Gamma} = 0 \}.
    
Thus the boundary term vanishes (:math:`\int_{\Gamma} \mu_{\mathrm{em}}^{-1} (\nabla \times \mathbf{E}) \cdot (\mathbf{v^*} \times \mathbf{n}) \mathrm{d}\Gamma = 0`) and the weak form of the time-harmonic Maxwell's equations is:

.. math::
    B(\mathbf{E}, \mathbf{v}) = &\int_{\Omega}  \mu_{\mathrm{em}}^{-1} (\nabla \times \mathbf{E}) \cdot (\nabla \times \mathbf{v^*}) \mathrm{d}\Omega + \int_{\Omega} \omega^2 \varepsilon_{\mathrm{em}_{eff}} \mathbf{E} \cdot \mathbf{v^*} \mathrm{d}\Omega  \\
    L(\mathbf{v}) = &\int_{\Omega} -i \omega \mathbf{J}_{\mathrm{ext}} \cdot \mathbf{v^*} \mathrm{d}\Omega .

Formally, :math:`\mathbf{E}` should also satisfy Gauss's law (:math:`\nabla \cdot \mathbf{D} = \rho_f`), but it is implicitly taken into account by the electromagnetic wave equation and holds in the weak form presented above.

Ultraweak Formulation
---------------------

The formulation above is the conventional weak form used as a reference for comparison with other other implementations of the time-harmonic Maxwell equations. In the DPG setting used in Lethe, the ultraweak formulation is used. The main practical difference is that the primal formulation above eliminates
the magnetic field to obtain a single second-order equation for the electric
field, whereas the ultraweak formulation solves for the electric and
magnetic fields simultaneously by weakening both equations without eliminating either one. Retaining both fields as independent unknowns is especially useful when the electromagnetic power dissipated in the medium must be evaluated, since some materials exhibit non-negligible magnetic losses in addition to dielectric losses; recovering :math:`\mathbf{H}` directly, rather than reconstructing it from
:math:`\mathbf{E}` after the fact, avoids an additional loss of accuracy in
that calculation. Additionally, studies have shown that the ultraweak formulation is more robust than the primal formulation in the context of the DPG method.


The starting point is the time-harmonic Maxwell system written in the form

.. math::

   \begin{align}
   \nabla \times \mathbf{E} - i\omega \mu \mathbf{H}       &= 0, \\
   \nabla \times \mathbf{H} + i\omega \varepsilon_{eff} \mathbf{E} &= \mathbf{J},
   \end{align}

where :math:`\mathbf{E}` is the electric field, :math:`\mathbf{H}` the
magnetic field, :math:`\mu` the permeability, :math:`\varepsilon_{eff}` the
effective permittivity, and :math:`\mathbf{J}` a prescribed current density.
After non-dimensionalization, one writes

.. math::

   \begin{align}
   \nabla \times \mathbf{E} - i\tilde{\omega} \mu_r \mathbf{H} &= 0, \\
   \nabla \times \mathbf{H} + i \tilde{\omega} (\varepsilon_r + i \sigma_r) \mathbf{E}
   &= \tilde{\mathbf{J}}.
   \end{align}

In what follows, the tildes are dropped to simplify the notation, but every
quantity remains dimensionless. The combination :math:`\varepsilon_{r,\text{eff}}
= \varepsilon_r + i \sigma_r` is the effective relative permittivity. The
system is supplemented with boundary conditions of three kinds:

.. math::

   \begin{align}
     \mathbf{n} \times \mathbf{E} &= \mathbf{n} \times \mathbf{E}_D
       &&\text{on } \Gamma_D, \\
     \mathbf{n} \times \mathbf{H} &= \mathbf{n} \times \mathbf{H}_N - \mathbf{J}_{\mathrm{s}}
       &&\text{on } \Gamma_N, \\
     \mathbf{n} \times \mathbf{H} + Z_{\mathrm{s}}^{-1}\mathbf{n} \times (\mathbf{E}\times \mathbf{n}) &= \mathbf{g}
       &&\text{on } \Gamma_R,
   \end{align}

where :math:`\Gamma_D`, :math:`\Gamma_N`, and :math:`\Gamma_R` are the
Dirichlet, Neumann, and Robin boundaries, :math:`\mathbf{n}` is the outward
unit normal, :math:`Z_{\mathrm{s}}` is the surface impedance on
:math:`\Gamma_R`, :math:`\mathbf{E}_D` and :math:`\mathbf{H}_N` are prescribed
fields on :math:`\Gamma_D`/:math:`\Gamma_N`, :math:`\mathbf{J}_{\mathrm{s}}` is
a prescribed surface current density, and :math:`\mathbf{g}` is a source term.
A catalog of the boundary condition types commonly used in electromagnetic
simulations, expressed in this dimensionless convention, is given in
:ref:`common-boundary-conditions` below.

Rather than eliminating the magnetic field to obtain a single second-order
equation for the electric field (as in the conventional primal formulation,
see :doc:`Finite Element Method for Time-Harmonic Maxwell's Equations
<fem_time_harmonic_maxwell>`), the DPG solver is built on the ultra-weak
formulation, which retains both :math:`\mathbf{E}` and :math:`\mathbf{H}` as
independent unknowns. This choice makes it straightforward to recover both
field solutions, for example to evaluate the electromagnetic power dissipated
in materials with non-negligible magnetic losses, and it gives the discrete
formulation the lowest inter-element regularity requirements among DPG
formulations, which is beneficial for high-frequency wave problems.

To obtain the ultra-weak formulation, both equations are multiplied by
independent test functions :math:`\mathbf{I}` and :math:`\mathbf{F}`,
respectively, and integrated by parts over the mesh :math:`\Omega_h`:

.. math::

   \begin{align}
      (\nabla \times\mathbf{I},\mathbf{E})_{\Omega_h} - ( \mathbf{I}, i\omega \mu_r\mathbf{H})_{\Omega_h} + \langle \mathbf{I}, \mathbf{n}\times \hat{\mathbf{E}} \rangle_{\partial \Omega_h}
        &= 0 , \\
      (\nabla \times \mathbf{F},\mathbf{H})_{\Omega_h} + (\mathbf{F}, i \omega \varepsilon_{r,\text{eff}} \mathbf{E})_{\Omega_h}
      + \langle \mathbf{F} , \mathbf{n}\times \hat{\mathbf{H}}\rangle_{\partial \Omega_h \backslash \Gamma_R }
      - \langle \mathbf{F} ,Z_\mathrm{s}^{-1} \mathbf{E}\rangle_{ \Gamma_R }
        &=(\mathbf{F},\mathbf{J})_{\Omega_h} -  \langle \mathbf{F} , \mathbf{g} \rangle_{ \Gamma_R },
   \end{align}

where the Robin boundary condition is applied to the electric field without
loss of generality, and :math:`\hat{\mathbf{E}}`, :math:`\hat{\mathbf{H}}` are
the trace unknowns living on the mesh skeleton. Each product above is
understood as complex-conjugate in its first argument:

.. math::

   \begin{align}
   (\mathbf{v},\mathbf{u})_{\Omega_h} &:= \sum_{K\in\Omega_h} \int_{K} \mathbf{v}^* \cdot \mathbf{u}\,\mathrm{d}x, &
   \langle \mathbf{v},\mathbf{u} \rangle_{\partial\Omega_h} &:= \sum_{K\in\Omega_h} \int_{\partial K} \mathbf{v}^* \cdot \mathbf{u}\, \mathrm{d}s.
   \end{align}

The trial and test spaces of the ultra-weak formulation are

.. math::

   \begin{align}
     \mathbf{E},\mathbf{H} &\in \mathcal{X}_{\mathrm{em}},
       & \mathcal{X}_{\mathrm{em}} &:= (L^2(\Omega_h))^3, \\
     \hat{\mathbf{E}},\hat{\mathbf{H}} &\in \hat{\mathcal{X}}_{\mathrm{em}},
       & \hat{\mathcal{X}}_{\mathrm{em}} &:= \bigl\{ \hat{\mathbf{x}} \in H^{-1/2}(\mathrm{curl},\Omega_h) : \mathbf{n}\times \hat{\mathbf{x}} = \mathbf{n}\times \mathbf{x}_\mathrm{D} \text{ on } \Gamma_D \bigr\}, \\
     \mathbf{F},\mathbf{I} &\in \mathcal{Y}_{\mathrm{em}},
       & \mathcal{Y}_{\mathrm{em}} &:= \bigl\{ \mathbf{y} \in H(\mathrm{curl},\Omega_h) : \mathbf{n}\times \mathbf{y} = 0 \text{ on } \partial\Omega \bigr\}.
   \end{align}

For an ultra-weak formulation, the Neumann boundary condition is, in fact, a
Dirichlet condition on the flux variable — here, the magnetic field
:math:`\mathbf{H}`. This is why only a Dirichlet-type condition needs to be
imposed on the trace space :math:`\hat{\mathcal{X}}_{\mathrm{em}}`, which
covers both the electric and magnetic trace unknowns: the Dirichlet and
Neumann relations above are then imposed directly on the relevant electric or
magnetic trace unknown.

To complete the discrete formulation, the test space is equipped with the
adjoint graph norm induced by the operator associated with the bilinear form
of the system:

.. math::

   \begin{align}
   \|(\mathbf{F}, \mathbf{I})\|^2_V &= \|  \nabla \times \mathbf{F} - i \omega \mu_r \mathbf{I} \|^2_{\Omega_h} + \| \nabla \times \mathbf{I} + i \omega \varepsilon_{r,\mathrm{eff}} \mathbf{F} \|^2_{\Omega_h} \\
     &\quad + \alpha\left(\|\mathbf{F}\|^2_{\Omega_h} + \|\mathbf{I}\|^2_{\Omega_h}\right)
      + \| \mathbf{n} \times \mathbf{I} + Z_\mathrm{s}^{-1} \mathbf{n} \times (\mathbf{F} \times \mathbf{n}) \|^2_{\Gamma_R},
   \end{align}

where :math:`\alpha \sim 1` is a strictly positive constant (:math:`\alpha=1`
in practice) introduced to ensure that the test norm remains localizable once
the test space is broken — this locality is what makes the element-wise
computation of the optimal test functions possible. The boundary term on
:math:`\Gamma_R` enforces the Robin boundary condition and is required to
obtain convergence of the method when this type of boundary condition is
present.


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
* **Impedance condition (Robin)**: :math:`\mathbf{n} \times \mathbf{H} + Z_\mathrm{s}^{-1} \mathbf{E}_\parallel = \mathbf{g}`, with :math:`Z_\mathrm{s}` the surface impedance of the boundary.
* **Absorbing / radiation condition (Silver–Müller)**: :math:`\mathbf{n} \times \mathbf{H} + \sqrt{\varepsilon_r/\mu_r}\, \mathbf{E}_\parallel = 0`. This is the impedance condition above with :math:`Z_\mathrm{s} = \sqrt{\mu_r/\varepsilon_r}` (the impedance of the exterior medium) and no source term.
* **Lossy / imperfect conductor**: :math:`\mathbf{n} \times \mathbf{H} + \sqrt{\varepsilon_{r,\mathrm{eff},2}/\mu_{r,2}}\, \mathbf{E}_\parallel = 0`, the same form as the radiation condition above, using the (possibly complex) effective properties :math:`\varepsilon_{r,\mathrm{eff},2}`, :math:`\mu_{r,2}` of the adjacent conducting medium.
* **Waveguide inlet port, TE**\ :math:`_{mn}` **mode**: :math:`\mathbf{n} \times \mathbf{H} + \dfrac{\mathbf{k} \cdot \mathbf{n}}{\omega \mu_r} \mathbf{E}_\parallel = \mathbf{n} \times \mathbf{H}_{inc} + \dfrac{\mathbf{k} \cdot \mathbf{n}}{\omega \mu_r} \mathbf{E}_{inc,\parallel}`.
* **Waveguide inlet port, TM**\ :math:`_{mn}` **mode**: :math:`\mathbf{n} \times \mathbf{H} + \dfrac{\omega \varepsilon_{r,\mathrm{eff}}}{\mathbf{k} \cdot \mathbf{n}} \mathbf{E}_\parallel = \mathbf{n} \times \mathbf{H}_{inc} + \dfrac{\omega \varepsilon_{r,\mathrm{eff}}}{\mathbf{k} \cdot \mathbf{n}} \mathbf{E}_{inc,\parallel}`.
* **Plane-wave inlet**: :math:`\mathbf{n} \times \mathbf{H} + \dfrac{\mathbf{k} \cdot \mathbf{n}}{\omega \mu_r} \mathbf{E}_\parallel = \mathbf{n} \times \mathbf{H}_{inc} + \sqrt{\varepsilon_r/\mu_r}\, \mathbf{E}_{inc,\parallel}`.

In each waveguide/plane-wave case, :math:`\mathbf{H}_{inc}` and
:math:`\mathbf{E}_{inc,\parallel}` denote the incident field of the excited
mode or plane wave, and :math:`\mathbf{k}` its wavevector.

