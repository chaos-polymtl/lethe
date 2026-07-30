==================================================
DPG formulation for time-harmonic Maxwell problems
==================================================

The discontinuous Petrov-Galerkin (DPG) method is a residual-minimizing finite
element strategy that is particularly effective for wave propagation and
high-frequency electromagnetic problems. Rather than using the same space for
trial and test functions, DPG uses different spaces and seeks the discrete
solution by minimizing the residual of the PDE in a norm induced by the test
space. For a review of the DPG method, see the following `review article <https://pdxscholar.library.pdx.edu/mth_fac/426/>`_ written by the original authors of the method, `Demkowicz <https://users.oden.utexas.edu/~leszek/>`_ and `Gopalakrishnan <https://web.pdx.edu/~gjay/>`_.

A Brief Introduction to DPG
===========================

The Discontinuous Petrov--Galerkin (DPG) method is a class of discontinuous finite
element methods that reformulates the standard variational problem as a minimum
residual problem. Given the abstract problem

.. math::

   \begin{aligned}
   \text{Find } u \in U \text{ such that } \\
   b(v,u)=l(v), \qquad \forall v \in V,
   \end{aligned}

the DPG method seeks an approximation in a finite-dimensional trial space
:math:`U_h \subset U` by minimizing the residual in the dual of the test space
:math:`V'`. The discrete solution is therefore defined as

.. math::

  \text{Find } u_h \in U_h \text{ such that } (R_V^{-1} (l - B u_h), R_V^{-1} B \delta u_h)_V = 0, \quad \forall \delta u_h \in U_h,

where :math:`B : U \rightarrow V'` is the operator satisfying

.. math::

   b(v,u)=\langle v, Bu \rangle_{V \times V'}.

The operator :math:`B` maps a trial function to its corresponding residual
functional. Here, :math:`\langle \cdot,\cdot \rangle_{V \times V'}` denotes the
duality pairing between the test space and its dual.

Since minimizing the residual directly in the dual norm is generally
impractical, the DPG method introduces the Riesz operator

.. math::

   \langle v, R_V w \rangle_{V \times V'}
   =
   (v,w)_V,
   \qquad
   \forall v,w \in V,

whose inverse maps functionals in :math:`V'` back into the test space
:math:`V`. The functional :math:`J` can then be rewritten as

.. math::

   J(w_h)
   =
   \frac{1}{2}
   \|R_V^{-1}(l-Bw_h)\|_V^2
   =
   \frac{1}{2}
   \left(
   R_V^{-1}(l-Bw_h),
   R_V^{-1}(l-Bw_h)
   \right)_V.

Taking the Gâteaux derivative of :math:`J` and setting it equal to zero yields

.. math::
   :label: minimization_problem

   (R_V^{-1}(l-Bu_h),
   R_V^{-1}B\delta u_h)_V
   =0,
   \qquad
   \forall \delta u_h \in U_h.

This naturally introduces the **error representation function**

.. math::

   \Psi = R_V^{-1}(l-Bu_h),

which belongs to the test space :math:`V`. Its norm is equal to the residual
measured in the dual norm,

.. math::

   \|l-Bu_h\|_{V'}
   =
   \|\Psi\|_V,

making it a natural a posteriori error estimator.

The error representation function also satisfies

.. math::

   (v,\Psi)_V
   =
   l(v)-b(v,u_h),
   \qquad
   \forall v \in V,

which, together with Equation :eq:`minimization_problem`, leads to the mixed DPG
formulation

.. math::

  \begin{cases}
     & \text{Find } u_h \in U_h \text{ and } \Psi \in V \text{ such that }     \\
     & (v,\Psi)_V + b(v,u_h) = l(v), \quad \forall v \in V,                    \\
     & (R_V^{-1} B w_h, \Psi)_V = b(w_h, \Psi) = 0, \quad \forall w_h \in U_h.
  \end{cases}

In practice, computing the exact optimal test functions would require inverting
the Riesz operator over the infinite-dimensional space :math:`V`, which is not
feasible. Instead, the test space is approximated by an enriched finite-dimensional
space :math:`V_r \subset V`. The polynomial order of :math:`V_r` is chosen higher
than that of the trial space :math:`U_h`.

Furthermore, the test space is chosen to be **broken**, meaning that its basis
functions are supported independently on each mesh element. This localization
makes the Riesz operator block diagonal, allowing the local problems defining the
optimal test functions to be solved independently on each element.

Breaking the test space introduces additional trace unknowns
:math:`\hat{u}_h \in \hat{U}_h`, defined on the mesh skeleton, to recover
inter-element conformity. The bilinear form becomes

.. math::

   b(v,(u_h,\hat{u}_h))
   =
   b(v,u_h)
   +
   \langle v,\hat{u}_h \rangle,

where :math:`\langle\cdot,\cdot\rangle` denotes the duality pairing between the
trace unknowns and the test functions.

The practical DPG formulation on a discretized domain
:math:`\Omega_h` is therefore

.. math::
   :label: final_dpg_formulation

   \begin{cases}
      & \text{Find } u_h \in U_h, \hat{u}_h \in \hat{U}_h \text{ and } \Psi \in V_r(\Omega_h) \text{ such that } \\
      & (v,\Psi)_V + b(v,u_h) + \langle v, \hat{u}_h \rangle = l(v), \quad \forall v \in V_r(\Omega_h),          \\
      & b(w_h, \Psi) = 0, \quad \forall w_h \in U_h,                                                             \\
      & \langle \hat{w}_h, \Psi \rangle = 0, \quad \forall \hat{w}_h \in \hat{U}_h.
   \end{cases}

This is the formulation implemented in Lethe. The use of an enriched broken test
space allows the optimal test functions to be computed locally, while the trace
unknowns ensure global continuity through the mesh skeleton.


Numerical implementation and discretization
===========================================

To take advantage of the extensive support for real-valued finite element
problems available in deal.II, all complex-valued fields are decomposed into
their real and imaginary components. Consequently, the entire linear system is
assembled and solved using standard double-precision floating-point arithmetic.

Discrete spaces
---------------

The discrete trial, trace, and test spaces are defined as

.. math::

   \begin{aligned}
   U_h
   &=
   \mathbf{E}_{\mathrm{re}}
   \times
   \mathbf{E}_{\mathrm{im}}
   \times
   \mathbf{H}_{\mathrm{re}}
   \times
   \mathbf{H}_{\mathrm{im}}, \\
   \hat{U}_h
   &=
   \hat{\mathbf{E}}_{\mathrm{re}}
   \times
   \hat{\mathbf{E}}_{\mathrm{im}}
   \times
   \hat{\mathbf{H}}_{\mathrm{re}}
   \times
   \hat{\mathbf{H}}_{\mathrm{im}}, \\
   V_r
   &=
   \mathbf{F}_{\mathrm{re}}
   \times
   \mathbf{F}_{\mathrm{im}}
   \times
   \mathbf{I}_{\mathrm{re}}
   \times
   \mathbf{I}_{\mathrm{im}}.
   \end{aligned}

Each component of the interior trial space is discretized using discontinuous
tensor-product Lagrange finite elements
:math:`(\mathcal{Q}^{-}_{p}\Lambda^3)`, providing an approximation of
:math:`L^2(\Omega_h)`.

The trace unknowns are discretized using first-kind Nédélec elements
:math:`(\mathcal{Q}^{-}_{p}\Lambda^1)`. The tangential trace operator
:math:`\mathrm{tr}_{\mathrm{curl},\top}` is applied, and the interior degrees
of freedom are discarded so that the resulting finite element space
approximates :math:`H^{-1/2}(\mathrm{curl},\Omega_h)`.

The enriched test space is also constructed from first-kind Nédélec elements,
but without applying the trace operator, making it conforming in
:math:`H(\mathrm{curl},\Omega_h)`. As required by the DPG method, the test space
uses a higher polynomial degree than the trial space,

.. math::

   p_{\mathrm{test}} = p + \Delta p,

where :math:`\Delta p` is the enrichment parameter. Throughout this
implementation,

.. math::

   \Delta p = 1.

The discrete fields are expanded in terms of their respective basis functions,

.. math::

   \begin{aligned}
   u_h^{(k)}(\mathbf{x})
   &=
   \sum_i
   \mathbf{w}^{(k)}_i
   \boldsymbol{\phi}_i(\mathbf{x}), \\
   \hat{u}_h^{(k)}(\mathbf{x})
   &=
   \sum_i
   \hat{\mathbf{w}}^{(k)}_i
   \hat{\boldsymbol{\phi}}_i(\mathbf{x}), \\
   v^{(k)}(\mathbf{x})
   &=
   \sum_i
   \mathbf{q}^{(k)}_i
   \boldsymbol{\psi}_i(\mathbf{x}),
   \end{aligned}

where

* :math:`\boldsymbol{\phi}_i` are the interior trial basis functions,
* :math:`\hat{\boldsymbol{\phi}}_i` are the trace basis functions,
* :math:`\boldsymbol{\psi}_i` are the test basis functions,

and :math:`\mathbf{w}_i`, :math:`\hat{\mathbf{w}}_i`, and
:math:`\mathbf{q}_i` denote their corresponding coefficients. The superscript
:math:`k \in \{1,2,3,4\}` identifies the field component
(:math:`\mathbf{E}_{\mathrm{re}}`,
:math:`\mathbf{E}_{\mathrm{im}}`,
:math:`\mathbf{H}_{\mathrm{re}}`,
:math:`\mathbf{H}_{\mathrm{im}}`).

Discrete linear system
----------------------

Combining the finite element discretization with the practical DPG formulation
:eq:`final_dpg_formulation` yields the block system

.. math::

   \begin{bmatrix}
   G & B & \hat{B} \\
   B^\dagger & 0 & 0 \\
   \hat{B}^\dagger & 0 & 0
   \end{bmatrix}
   \begin{bmatrix}
   \Psi^r \\
   u_h \\
   \hat{u}_h
   \end{bmatrix}
   =
   \begin{bmatrix}
   l \\
   0 \\
   0
   \end{bmatrix}.

The matrices are defined as

.. math::

   \begin{aligned}
   G_{(i,k),(j,l)}
   &= (\psi_i^{(k)},\psi_j^{(l)})_V, \\
   B_{(i,k),(j,l)}
   &= b(\psi_i^{(k)},\phi_j^{(l)}), \\
   \hat{B}_{(i,k),(j,l)}
   &= \langle
   \psi_i^{(k)},
   \hat{\phi}_j^{(l)}
   \rangle,
   \end{aligned}

while the load vector is

.. math::

   l_{(i,k)} = l(\psi_i^{(k)}).

The indices :math:`(i,k)` combine the finite element basis index with the field
component.

Static condensation
-------------------

Because the test space is broken, the Gram matrix :math:`G` is block diagonal,
allowing the error representation function to be eliminated locally. This
produces the condensed system

.. math::

   \begin{bmatrix}
   B^\dagger G^{-1}B &
   B^\dagger G^{-1}\hat{B} \\
   \hat{B}^\dagger G^{-1}B &
   \hat{B}^\dagger G^{-1}\hat{B}
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

A second static condensation is also performed since interior elements are discontinuous, which eliminates the interior unknowns from the linear system. Defining

.. math::

   \begin{aligned}
   M_1 &= B^\dagger G^{-1}B, \\
   M_2 &= B^\dagger G^{-1}\hat{B}, \\
   M_3 &= \hat{B}^\dagger G^{-1}\hat{B}, \\
   M_4 &= B^\dagger G^{-1}, \\
   M_5 &= \hat{B}^\dagger G^{-1},
   \end{aligned}

the final system becomes

.. math::

   (M_3-M_2^\dagger M_1^{-1}M_2)\hat{u}_h
   =
   (M_5-M_2^\dagger M_1^{-1}M_4)l.

This Schur complement system is the one solved by Lethe's DPG solver. 

Recovery of local fields
------------------------

Once the trace unknowns have been computed, the interior fields and the error
representation function are recovered independently on each element,

.. math::

   u_h
   =
   M_1^{-1}(M_4l-M_2\hat{u}_h),

.. math::

   \Psi^r
   =
   G^{-1}(Bu_h+\hat{B}\hat{u}_h-l).

The local residual norm,

.. math::

   \|\Psi^r\|_V
   =
   \sqrt{(\Psi^r)^\dagger G\Psi^r},

provides the DPG a posteriori error estimator used for adaptive mesh
refinement.

Polynomial spaces
=================

The DPG method relies on several Sobolev, broken, and trace spaces. Their
precise definition is important when selecting compatible finite element
spaces for the trial, test, and trace unknowns. This page summarizes the spaces used by the implementation mentioned above. 

Sobolev spaces
--------------

The implementation is based on the standard Sobolev spaces defined over the domain :math:`\Omega`.

.. math::
   \begin{aligned}
   H^1(\Omega)
   &=
   \left\{
   u :
   \Omega \rightarrow \mathbb{R}
   \text{ (or } \mathbb{C}\text{)}
   \;:\;
   u,\nabla u \in (L^2(\Omega))^d
   \right\},\\
   H(\mathrm{curl},\Omega)
   &=
   \left\{
   \mathbf{E}
   :
   \Omega \rightarrow \mathbb{R}^d
   \text{ (or } \mathbb{C}^d\text{)}
   \;:\;
   \mathbf{E},
   \nabla\times\mathbf{E}
   \in (L^2(\Omega))^d
   \right\},\\
   H(\mathrm{div},\Omega)
   &=
   \left\{
   \boldsymbol{\sigma}
   :
   \Omega \rightarrow \mathbb{R}^d
   \text{ (or } \mathbb{C}^d\text{)}
   \;:\;
   \boldsymbol{\sigma}
   \in (L^2(\Omega))^d,
   \;
   \nabla\cdot\boldsymbol{\sigma}
   \in L^2(\Omega)
   \right\},\\
   L^2(\Omega)
   &=
   \left\{
   q :
   \Omega \rightarrow \mathbb{R}
   \text{ (or } \mathbb{C}\text{)}
   \;:\;
   \|q\|_{L^2(\Omega)} < \infty
   \right\}.
   \end{aligned}

Broken spaces
-------------

The current implementation supports three-dimensional problems only.
Accordingly, the broken spaces are defined over the mesh
:math:`\Omega_h` as

.. math::
   \begin{aligned}
   H^1(\Omega_h)
   &=
   \prod_{K\in\Omega_h}H^1(K),\\
   H(\mathrm{curl},\Omega_h)
   &=
   \prod_{K\in\Omega_h}
   H(\mathrm{curl},K),\\
   H(\mathrm{div},\Omega_h)
   &=
   \prod_{K\in\Omega_h}
   H(\mathrm{div},K).
   \end{aligned}

The product :math:`\prod_{K\in\Omega_h}` denotes the Cartesian product over all mesh elements, while :math:`u|_K` denotes the restriction of a function to a single element. Unlike the previous spaces,

.. math::

   L^2(\Omega_h)
   \equiv
   L^2(\Omega),

since :math:`L^2` functions do not require continuity between neighboring
elements.

Trace operators
---------------

The trace spaces are defined using the standard trace operators acting on the
mesh skeleton :math:`\partial\Omega_h`.

.. math::

   \begin{aligned}
   \mathrm{tr}_{\mathrm{grad}}
   &: H^1(\Omega_h)
   \rightarrow
   \prod_{K\in\Omega_h} H^{1/2}(\partial K),
   &
   \mathrm{tr}_{\mathrm{grad}}^K(u)
   &= u|_{\partial K},
   \\[0.8em]
   \mathrm{tr}_{\mathrm{curl},\top}
   &: H(\mathrm{curl},\Omega_h)
   \rightarrow
   \prod_{K\in\Omega_h}
   H^{-1/2}(\mathrm{curl},\partial K),
   &
   \mathrm{tr}_{\mathrm{curl},\top}^K(\mathbf{E})
   &= ((\mathbf{n}\times\mathbf{E})\times\mathbf{n})|_{\partial K},
   \\[0.8em]
   \mathrm{tr}_{\mathrm{curl},\dashv}
   &: H(\mathrm{curl},\Omega_h)
   \rightarrow
   \prod_{K\in\Omega_h}
   H^{-1/2}(\mathrm{div},\partial K),
   &
   \mathrm{tr}_{\mathrm{curl},\dashv}^K(\mathbf{E})
   &= (\mathbf{n}\times\mathbf{E})|_{\partial K},
   \\[0.8em]
   \mathrm{tr}_{\mathrm{div}}
   &: H(\mathrm{div},\Omega_h)
   \rightarrow
   \prod_{K\in\Omega_h}
   H^{-1/2}(\partial K),
   &
   \mathrm{tr}_{\mathrm{div}}^K(\boldsymbol{\sigma})
   &= (\mathbf{n}\cdot\boldsymbol{\sigma})|_{\partial K}.
   \end{aligned}

Here, :math:`\mathbf{n}` denotes the outward unit normal vector on the boundary
of each mesh element.

Trace spaces
------------

The trace operators define the interface spaces used by the DPG formulation,

.. math::
   \begin{aligned}
   H^{1/2}(\partial\Omega_h)
   &=
   \mathrm{tr}_{\mathrm{grad}}
   \left(
   H^1(\Omega)
   \right),\\
   H^{-1/2}(\mathrm{curl},\partial\Omega_h)
   &=
   \mathrm{tr}_{\mathrm{curl},\top}
   \left(
   H(\mathrm{curl},\Omega)
   \right),\\
   H^{-1/2}(\mathrm{div},\partial\Omega_h)
   &=
   \mathrm{tr}_{\mathrm{curl},\dashv}
   \left(
   H(\mathrm{curl},\Omega)
   \right),\\
   H^{-1/2}(\partial\Omega_h)
   &=
   \mathrm{tr}_{\mathrm{div}}
   \left(
   H(\mathrm{div},\Omega)
   \right).
   \end{aligned}

These interface spaces provide the mathematical foundation for the trace
unknowns introduced by the practical DPG formulation.
