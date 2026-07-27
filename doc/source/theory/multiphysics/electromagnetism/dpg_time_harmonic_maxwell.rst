===============================================
DPG formulation for time-harmonic Maxwell problems
===============================================

The discontinuous Petrov-Galerkin (DPG) method is a residual-minimizing finite
element strategy that is particularly effective for wave propagation and
high-frequency electromagnetic problems. Rather than using the same space for
trial and test functions, DPG uses different spaces and seeks the discrete
solution by minimizing the residual of the PDE in a norm induced by the test
space.

A Brief Introduction to DPG
----------------------------

For an abstract problem

.. math::

   b(v,u) = l(v) \quad \forall v \in V,

DPG seeks an approximate trial function :math:`u_h \in U_h` that minimizes the
residual in the dual norm of the test space,

.. math::

   u_h = \arg\min_{w_h \in U_h} \frac{1}{2} \| l - B w_h \|_{V'}^2,

where :math:`B : U \to V'` is the operator associated with the bilinear form
:math:`b`. To make this practical, one introduces the Riesz operator
:math:`R_V : V \to V'` and the corresponding error representation function
:math:`\Psi = R_V^{-1}(l - B u_h)`. The quantity :math:`\Psi` is a test-space
function that measures the residual, in the sense that
:math:`\|l - Bu_h\|_{V'} = \|\Psi\|_V`.

The error representation function satisfies the relation

.. math::

   (v,\Psi)_V = l(v) - b(v,u_h), \quad \forall v \in V,

which leads to the following equivalent **mixed formulation** of the DPG
method: find :math:`u_h \in U_h` and :math:`\Psi \in V` such that

.. math::

   \begin{align}
     (v,\Psi)_V + b(v,u_h) &= l(v),   & \forall v \in V,     \\
     b(w_h,\Psi)           &= 0,      & \forall w_h \in U_h.
   \end{align}

In practice, the test space is chosen as an enriched broken space so that the
Riesz operator becomes block-diagonal and can be inverted locally on each mesh
element. Trace unknowns on the mesh skeleton are then introduced to recover the
conformity of the solution across element interfaces. This is the key idea
behind ultra-weak DPG formulations and is the setting used for time-harmonic
Maxwell problems.

Introducing the trace unknowns :math:`\hat{u}_h \in \hat{U}_h` modifies the
bilinear form to :math:`b(v,(u_h,\hat{u}_h)) = b(v,u_h) + \langle
v,\hat{u}_h \rangle`, where :math:`\langle \cdot,\cdot \rangle` denotes the
pairing between the trace unknowns and the test functions on the mesh
skeleton. The DPG problem posed on a discretized domain :math:`\Omega_h` then
reads: find :math:`u_h \in U_h`, :math:`\hat{u}_h \in \hat{U}_h` and
:math:`\Psi \in V_r(\Omega_h)` such that

.. math::

   \begin{align}
     (v,\Psi)_V + b(v,u_h) + \langle v, \hat{u}_h \rangle & = l(v),  & \forall v \in V_r(\Omega_h),          \\
     b(w_h, \Psi)                                          & = 0,    & \forall w_h \in U_h,                  \\
     \langle \hat{w}_h, \Psi \rangle                       & = 0,    & \forall \hat{w}_h \in \hat{U}_h,
   \end{align}

where :math:`V_r \subset V` is a finite-dimensional test space enriched
relative to the trial space, meaning that the polynomial degree used for
:math:`V_r` is higher than the one used for :math:`U_h`. This enrichment is
what makes it possible to approximate the (otherwise infinite-dimensional)
optimal test functions in practice.

Ultra-weak Formulation for the Maxwell Problem
------------------------------------------------

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

Numerical Implementation and Discretization
----------------------------------------------

The DPG solver discretizes the ultra-weak formulation using hexahedral finite
elements. To exploit the broader range of features available for real-valued
problems, the trial functions (:math:`\mathbf{E}`, :math:`\mathbf{H}`,
:math:`\hat{\mathbf{E}}`, :math:`\hat{\mathbf{H}}`) and test functions
(:math:`\mathbf{F}`, :math:`\mathbf{I}`) are decomposed into their real (e.g.,
:math:`\mathbf{E}_{\mathrm{re}}`) and imaginary (e.g.,
:math:`\mathbf{E}_{\mathrm{im}}`) parts, and the resulting system is solved
entirely with real (double-precision) arithmetic.

The discrete trial, trace, and test spaces are

.. math::

   \begin{align}
     U_h     &= \mathbf{E}_\mathrm{re} \times \mathbf{E}_\mathrm{im} \times \mathbf{H}_\mathrm{re} \times \mathbf{H}_\mathrm{im}, \\
     \hat{U}_h &= \hat{\mathbf{E}}_\mathrm{re} \times \hat{\mathbf{E}}_\mathrm{im} \times \hat{\mathbf{H}}_\mathrm{re} \times \hat{\mathbf{H}}_\mathrm{im}, \\
     V_r     &= \mathbf{F}_\mathrm{re} \times \mathbf{F}_\mathrm{im} \times \mathbf{I}_\mathrm{re} \times \mathbf{I}_\mathrm{im}.
   \end{align}

Each component of the interior trial space :math:`u_h \in U_h` is discretized
using discontinuous finite elements built from tensor products of Lagrange
polynomials of degree :math:`p` in each spatial direction, giving a discrete
approximation of :math:`L^2(\Omega_h)`. The trace unknowns :math:`\hat{u}_h
\in \hat{U}_h` are discretized using Nédélec elements of the first kind of
degree :math:`p`, to which a tangential trace operator is applied (discarding
the interior degrees of freedom) so that the resulting space approximates
:math:`H^{-1/2}(\mathrm{curl}, \Omega_h)`. The test space :math:`v \in V_r` is
also built from Nédélec elements of the first kind, but without the trace
operator, so that it remains conforming to :math:`H(\mathrm{curl},
\Omega_h)`; its polynomial degree is set to :math:`p+\Delta p`, where
:math:`\Delta p \in \mathbb{N}^+` is the enrichment parameter of the method (a
value of :math:`\Delta p = 1` is used in practice).

Expanding the trial, trace, and test functions in terms of their basis
functions,

.. math::

   \begin{align}
      u_h^{(k)}(\mathbf{x})     &= \sum_{i} \mathbf{w}^{(k)}_{i} \boldsymbol{\phi}_{i}(\mathbf{x}), \\
      \hat{u}_h^{(k)}(\mathbf{x}) &= \sum_{i} \hat{\mathbf{w}}^{(k)}_{i} \hat{\boldsymbol{\phi}}_{i}(\mathbf{x}), \\
      v^{(k)}(\mathbf{x})       &= \sum_{i} \mathbf{q}^{(k)}_{i} \boldsymbol{\psi}_{i}(\mathbf{x}),
   \end{align}

where the superscript :math:`k\in \{1,2,3,4\}` indexes the four field
components (:math:`\mathbf{E}_\mathrm{re}, \mathbf{E}_\mathrm{im},
\mathbf{H}_\mathrm{re}, \mathbf{H}_\mathrm{im}`), leads to the following
linear system:

.. math::

   \begin{equation}
   \begin{bmatrix}
     G               & B & \hat{B} \\
     B^\dagger       & 0 & 0       \\
     \hat{B}^\dagger & 0 & 0
   \end{bmatrix}
   \begin{bmatrix}
     \Psi^r \\
     u_h    \\
     \hat{u}_h
   \end{bmatrix}
   =
   \begin{bmatrix}
     l \\
     0 \\
     0
   \end{bmatrix},
   \end{equation}

where :math:`B` is the interior matrix, :math:`\hat{B}` the interface matrix,
:math:`G` the Gram matrix, and :math:`l` the load vector, each assembled from
the bilinear and linear forms of the ultra-weak formulation and the test norm
defined above.

Because the Gram matrix is block-diagonal and local to each element — a
consequence of the broken test space — the error representation function
:math:`\Psi^r` can be eliminated locally, leaving

.. math::

   \begin{equation}
   \begin{bmatrix}
     B^\dagger G^{-1}B       & B^\dagger G^{-1}\hat{B}       \\
     \hat{B}^\dagger G^{-1}B & \hat{B}^\dagger G^{-1}\hat{B}
   \end{bmatrix}
   \begin{bmatrix}
     u_h \\
     \hat{u}_h
   \end{bmatrix}
   =
   \begin{bmatrix}
     B^\dagger G^{-1}l \\
     \hat{B}^\dagger G^{-1}l
   \end{bmatrix}.
   \end{equation}

A second condensation step then eliminates the interior unknowns
:math:`u_h`, leaving a system posed only in terms of the trace unknowns
:math:`\hat{u}_h`. Writing :math:`M_1 = B^\dagger G^{-1}B`, :math:`M_2 =
B^\dagger G^{-1}\hat{B}`, :math:`M_3 = \hat{B}^\dagger G^{-1}\hat{B}`,
:math:`M_4 = B^\dagger G^{-1}`, and :math:`M_5 = \hat{B}^\dagger G^{-1}`, this
final system reads

.. math::

   \begin{equation}
   (M_3 - M_2^\dagger M_1^{-1} M_2) \hat{u}_h = (M_5 - M_2^\dagger M_1^{-1} M_4) l,
   \end{equation}

which is the Schur complement associated with eliminating the interior
unknowns, and is the system actually solved for :math:`\hat{u}_h`. Once solved,
the interior solution and the error representation function are recovered
locally, element by element:

.. math::

   \begin{equation}
   u_h = M_1^{-1} (M_4 l - M_2 \hat{u}_h), \qquad \Psi^r = G^{-1} (B u_h + \hat{B} \hat{u}_h - l).
   \end{equation}

The residual norm on each element, :math:`\|\Psi^r\|_V = \sqrt{(\Psi^r)^\dagger
G \Psi^r}`, is used as the error estimator driving adaptive mesh refinement.

.. note::

   The condensed system in :math:`\hat{u}_h` is solved with a CG solver
   without preconditioning; standard black-box preconditioners (e.g.,
   algebraic multigrid or incomplete LU factorization) were found ineffective
   for this system in practice. The current formulation of the Robin boundary
   condition also leaves certain magnetic-field degrees of freedom
   unconstrained at those boundaries; this has not been observed to affect
   convergence for polynomial degree :math:`p \geq 1` when using an iterative
   linear solver.

.. _common-boundary-conditions:

Common Boundary Conditions
-----------------------------

The general Dirichlet, Neumann, and Robin conditions introduced above cover a
range of standard boundary condition types encountered in electromagnetic
simulations, each corresponding to a particular choice of :math:`\mathbf{E}_D`,
:math:`\mathbf{H}_N`, :math:`Z_\mathrm{s}`, and :math:`\mathbf{g}`. Writing
:math:`\mathbf{E}_\parallel = \mathbf{n} \times (\mathbf{E} \times \mathbf{n})`
for the tangential component of a field (and similarly
:math:`\mathbf{H}_\parallel`), the list below is not exhaustive but covers the
conditions commonly used in practice, including the ones used in the
:ref:`waveguide problem <time-harmonic-maxwell-waveguide-setting>` described
further below.

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

.. _time-harmonic-maxwell-waveguide-setting:

Waveguide Setting and Boundary Conditions
--------------------------------------------

A typical application is the propagation of a guided mode in a rectangular
waveguide. The geometry is assumed to be a box of dimensions :math:`a` and
:math:`b` in the :math:`x` and :math:`y` directions, respectively, and the
fields are excited by a transverse electric (TE) or transverse magnetic (TM)
mode. The TE\ :math:`_{mn}` (or TM\ :math:`_{mn}`) mode is characterized by

.. math::

   k_x = \frac{m\pi}{a}, \qquad k_y = \frac{n\pi}{b},

and the cutoff wavenumber is

.. math::

   k_c = \sqrt{k_x^2 + k_y^2}.

The longitudinal wavenumber is then

.. math::

   k_z = \sqrt{\omega^2 \mu \varepsilon_{eff} - k_x^2 - k_y^2}.

The waveguide problem uses three of the boundary condition types introduced
above. On :math:`\Gamma_1` one imposes a perfect electric conductor (PEC)
condition,

.. math::

   \mathbf{n} \times \mathbf{E} = 0,

on :math:`\Gamma_2` a waveguide inlet port condition prescribes the incident
mode through an impedance-type relation,

.. math::

   \mathbf{n} \times \mathbf{H} + \frac{k_z}{\omega \mu_r} \mathbf{E}_{\parallel} = g,

and on :math:`\Gamma_3` an absorbing (Silver–Müller radiation) condition is
used with :math:`g = 0`. These boundary operators are naturally handled by the
ultra-weak DPG formulation, since the boundary terms are treated explicitly
through the trace unknowns.

For completeness, the incident TE fields can be written in compact form as

.. math::

   \mathbf{E}_{inc} = i\frac{\omega \mu_r}{k_c^2}
   \begin{bmatrix}
   -k_y \cos(k_x x) \sin(k_y y) \\
   k_x \sin(k_x x) \cos(k_y y) \\
   0
   \end{bmatrix} e^{i k_z z},

.. math::

   \mathbf{H}_{inc} =
   \begin{bmatrix}
   - i \dfrac{k_z k_x}{k_c^2} \sin(k_x x) \cos(k_y y) \\
   - i \dfrac{k_z k_y}{k_c^2} \cos(k_x x) \sin(k_y y) \\
   \cos(k_x x) \cos(k_y y)
   \end{bmatrix} e^{i k_z z},

and the equivalent TM incident fields are

.. math::

   \mathbf{E}_{inc} =
   \begin{bmatrix}
   i \dfrac{k_z k_x}{k_c^2} \cos(k_x x) \sin(k_y y) \\
   i \dfrac{k_z k_y}{k_c^2} \sin(k_x x) \cos(k_y y) \\
   \sin(k_x x) \sin(k_y y)
   \end{bmatrix} e^{i k_z z},

.. math::

   \mathbf{H}_{inc} = i\frac{\omega \varepsilon_{r,\mathrm{eff}}}{k_c^2}
   \begin{bmatrix}
   k_y \sin(k_x x) \cos(k_y y) \\
   -k_x \cos(k_x x) \sin(k_y y) \\
   0
   \end{bmatrix} e^{i k_z z}.

.. note::

   The TM\ :math:`_{m0}` and TM\ :math:`_{0n}` modes do not exist, so at least
   one of :math:`m,n` must be nonzero for a TM mode to be defined.

These expressions are useful for defining the incoming wave on the port and
for comparing numerical solutions with the expected guided-mode behavior.

Waveguide Port Power
~~~~~~~~~~~~~~~~~~~~~~

Most of the time, the power carried by the guided wave is prescribed at the
waveguide port rather than the field amplitude itself. This power is related
to the electric and magnetic field amplitudes through the Poynting vector
:math:`\mathbf{S} = \frac{1}{2} \mathbf{E} \times \mathbf{H}^*`. The
time-averaged power flowing through the port area :math:`A` is

.. math::

   \langle P \rangle_{avg} = \frac{1}{2} \int_A \mathbf{S} \cdot \mathbf{n}\, \mathrm{d}A
   = \frac{1}{2} \int_A \Re{(\mathbf{E} \times \mathbf{H}^*)} \cdot \mathbf{n}\, \mathrm{d}A.

For a rectangular waveguide port excited by a TE\ :math:`_{mn}` mode, this
integral can be evaluated analytically:

.. math::

   \langle P \rangle_{avg} =
   \begin{cases}
     \dfrac{\omega \Re(\mu k_{z}^*)\, ab\, H_{0}^{2}}{8\left ( k_{x}^{2}+k_{y}^{2} \right )}, & m,n\neq 0, \\[8pt]
     \dfrac{\omega \Re(\mu k_{z}^*)\, ab\, H_{0}^{2}}{4\left ( k_{x}^{2}+k_{y}^{2} \right )}, & m=0 \text{ or } n=0,
   \end{cases}

where :math:`H_0 = E_0 / Z_0` is the amplitude of the magnetic field at the
waveguide port. The equivalent expression for a TM\ :math:`_{mn}` mode (recall
that TM\ :math:`_{m0}` and TM\ :math:`_{0n}` modes do not exist, so only
:math:`m,n\neq 0` applies) is

.. math::

   \langle P \rangle_{avg} = \frac{\omega \Re(\varepsilon^* k_{z})\, ab\, E_{0}^{2}}{8\left ( k_{x}^{2}+k_{y}^{2} \right )}, \quad m,n\neq 0.

.. note::

   These power expressions are given here in *dimensional* form (using
   :math:`\mu`, :math:`\varepsilon`, etc.) since this is how the input power is
   naturally specified by the user. They are the relations used, and
   inverted, to determine the reference field amplitude :math:`E_0` from a
   user-prescribed input power, before non-dimensionalizing the rest of the
   problem — this is exactly what the ``electromagnetic scaling type =
   power`` option does in the time-harmonic Maxwell parameter subsection.