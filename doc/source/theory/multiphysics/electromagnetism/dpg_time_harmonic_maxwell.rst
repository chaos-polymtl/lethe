==================================================
DPG formulation for time-harmonic Maxwell problems
==================================================

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
