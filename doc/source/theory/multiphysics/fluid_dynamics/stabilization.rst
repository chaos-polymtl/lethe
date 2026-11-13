==============================
On the Need for Stabilization
==============================

It is well known that, in the absence of stabilization, the choice of the velocity and pressure finite element spaces must be done to ensure that the `Ladyzenskaya-Babuska-Brezzi (LBB) inf-sup <https://en.wikipedia.org/wiki/Ladyzhenskaya%E2%80%93Babu%C5%A1ka%E2%80%93Brezzi_condition>`_ condition is met. It is also known that the Galerkin approximation of the Navier–Stokes equations may fail in convection dominated flows, for which there are boundary layers where the velocity solution and its gradient exhibit rapid variation over short length scales. In these regions, the classical Galerkin approach may lead to numerical oscillations which may greatly pollute the solution over the entire domain. Stabilization of the elements, which is used in Lethe, can circumvent the limitations of the classical Galerkin approach and alleviate the  need to use LBB stable elements. This notably allows for the use of equal order elements (such as P1-P1).


-----------------------------------
SUPG/PSPG Monolithic Formulation
-----------------------------------

This approach consists of two new terms that are added to the classic weak formulation. The first one is known as
SUPG (Streamline-Upwind/ Petrov-Galerkin) formulation and it is in charge of stabilizing the formulation for convection-dominated flows. The second one is known as PSPG (Pressure-Stabilizing/ Petrov-Galerkin) and as its name indicates, it reduces pressure oscillations when using equal order elements (e.g. Q1 elements for the velocity and Q1 elements for the pressure). Thus, the new weak form of the Navier-Stokes equation is given by:

.. math::

  &\int_{\Omega}  q  \partial_l u_l \mathrm{d}\Omega + \sum_{K} \int_{\Omega_K} \left( \partial_t u_k + u_l \partial_l u_k + \partial_k p - \nu \partial_l \partial_l u_k - f_k \right) \cdot \left(\tau_{PSPG} \partial_l q \right) \mathrm{d}\Omega_K  = 0
  \\
  &\int_{\Omega}  v_k \left(\partial_t u_k+ u_l \partial_l u_k - f_k \right) \mathrm{d}\Omega - \int_{\Omega} \left( \partial_k \right) v_k p \mathrm{d}\Omega  + \nu \int_{\Omega} \left( \partial_l v_k \right) \left( \partial_l u_k  \right) \mathrm{d}\Omega   
  \\
  & + \sum_{K} \int_{\Omega_K} \left( \partial_t u_k + u_l \partial_l u_k + \partial_k p - \nu \partial_l \partial_l u_k - f_k \right) \cdot \left(\tau_{SUPG} u_l \partial_l v_k \right) \mathrm{d}\Omega_K =0

This formulation is consistent as the added terms involve the residual and if we substitute the exact solution, the stabilization terms vanish. In the literature, one can find different definitions of the stabilization parameters :math:`\tau_{SUPG}` and :math:`\tau_{PSPG}`. Lethe uses the definition by `Tezduyar (1992) <https://doi.org/10.1016/0045-7825(92)90141-6>`_ for both stabilization parameters. In the case of a transient problem:

.. math::

   \tau = \left[ \left( \frac{1}{\Delta t} \right)^{2} + \left( \frac{2 |\mathrm{u}|}{h_{conv}} \right)^{2} + 9 \left( \frac{4 \nu}{h^2_{diff}} \right)^{2} \right]^{-1/2}

where :math:`\Delta t` is the time step, and :math:`h_{conv}` and :math:`h_{diff}` are the size of the element related to the convection transport and diffusion mechanism, respectively. In Lethe, both element sizes are set to the diameter of a sphere having a volume equivalent to that of the cell. In the case of stationary problems, the following expression is used: 

.. math::

   \tau = \left[ \left( \frac{2 |\mathrm{u}|}{h_{conv}} \right)^{2} + 9 \left(\frac{4 \nu}{h^2_{diff}} \right)^{2} \right]^{-1/2}

To solve the non-linear problem, Lethe uses again the Newton-Raphson method, however, the Jacobian and the residual are now of the following form: 

.. math::
    
  \mathbf{\mathcal{J}} &= \left[ \begin{matrix} 	A^* & B^{*T}  \\[0.3em]	B^* & S^* \end{matrix} \right] \\
  \mathbf{\mathcal{R}} &=  \left[ \begin{matrix} \mathbf{\mathcal{R}}_{v}^*   \\[0.3em]		\mathbf{\mathcal{R}}_{q}^*  \end{matrix} \right]
  
As it can be seen, there is an additional matrix :math:`S^*` that appears in the linear system and eliminates the zero block on the diagonal of the Jacobian matrix.


------------------------------------
Grad-Div Block Formulation
------------------------------------

This approach builds on the work of `Heister et al. (2012) <https://onlinelibrary.wiley.com/doi/10.1002/fld.3654>`_. and can be seen as a natural parallel extension of the `Step-57 of deal.ii <https://www.dealii.org/current/doxygen/deal.II/step_57.html>`_. An additional term is added to the classic weak form of the Navier-Stokes equations: 

.. math::

  &\int_{\Omega}  q  \partial_l u_l \mathrm{d}\Omega =0 
  \\
  &\int_{\Omega}  v_k \left(\partial_t u_k+ u_l \partial_l u_k - f_k \right) \mathrm{d}\Omega  - \int_{\Omega} \left( \partial_k \right) v_k p \mathrm{d}\Omega  
  \\
  &+ \nu \int_{\Omega} \left( \partial_l v_k \right) \left( \partial_l u_k  \right) \mathrm{d}\Omega  + \sum_K \gamma \int_{\Omega_K} \partial_l u_l \partial_k v_k \mathrm{d}\Omega_K = 0

where :math:`\gamma` is an additional parameter that can be related to the augmented lagrangian formulation. The additional stabilization term improves the numerical accuracy of the solution and helps reduce oscillations for convection-dominated flows. In general, the optimal value for :math:`\gamma` depends on the solution on each element and it is therefore, problem dependent. In Lethe the value of :math:`\gamma` is equal to :math:`1`. In this case, the linear system to be solved in each non-linear iteration has the same structure as the one obtained with the classic weak formulation. Therefore, a good preconditioning is necessary to solve the linear system at each nonlinear iteration. More on this topic is found in the linear solvers section.


-----------------------------------
Galerkin Least-Squares Formulation
-----------------------------------

The GLS formulation is built as a generalization of the stabilization procedure in the SUPG/PSPG formulation. It consists of both the SUPG and PSPG terms as well as two additional terms: a term based on a Least-Squares term based on the momentum equation and a term based on a Least-Squares term based on the incompressibility constraint. All terms can be seen in the following weak form:

.. math::

  &\int_{\Omega}  q  \partial_l u_l \mathrm{d}\Omega + \sum_{K} \int_{\Omega_K} \left( \partial_t u_k + u_l \partial_l u_k + \partial_k p - \nu \partial_l \partial_l u_k - f_k \right) \cdot \left(\tau_{PSPG} \partial_l q \right) \mathrm{d}\Omega_K  = 0
  \\
  &\int_{\Omega}  v_k \left(\partial_t u_k+ u_l \partial_l u_k - f_k \right) \mathrm{d}\Omega - \int_{\Omega} \left( \partial_k \right) v_k p \mathrm{d}\Omega  + \nu \int_{\Omega} \left( \partial_l v_k \right) \left( \partial_l u_k  \right) \mathrm{d}\Omega   
  \\
  & + \sum_{K} \int_{\Omega_K} \left( \partial_t u_k + u_l \partial_l u_k + \partial_k p - \nu \partial_l \partial_l u_k - f_k \right) \cdot \left(\tau_{SUPG} u_l \partial_l v_k \right) \mathrm{d}\Omega_K
  \\
  & - \sum_{K} \int_{\Omega_K} \left( \partial_t u_k + u_l \partial_l u_k + \partial_k p - \nu \partial_l \partial_l u_k - f_k \right) \cdot \left(\tau_{GLS} \nu \partial_l \partial_l v_k \right) \mathrm{d}\Omega_K
  \\
  & + \sum_{K} \int_{\Omega_K} \tau_{LSIC} (\partial_l v_l) (\partial_l u_l) \mathrm{d}\Omega_K = 0

This is the version of the GLS stabilization if the finite element method is only used for the spatial discretization and no time-space finite element formulation is used as is the case in Lethe. The stabilization parameters are taken to be the same for the SUPG, PSPG and GLS terms, and given by the same :math:`\tau` expressions presented in the SUPG/PSPG section. In the case of the LSIC term, the stabilization parameter is defined as:

.. math::

   \tau_{LSIC} = \frac{|\mathrm{u}| h}{2}

The non-linear problem is solved in the same fashion and the structure of the Jacobian is the same one as that of the SUPG/PSPG formulation.


-----------------------------------------------------------
Residual-Based Variational Multiscale (RBVMS) Formulation
-----------------------------------------------------------

The residual-based variational multiscale (RBVMS) formulation of Bazilevs *et al.* [#bazilevs2007]_ rests on the same principle as the SUPG/PSPG and GLS formulations. Terms proportional to the residuals of the momentum and continuity equations are added to the Galerkin weak form, so that the exact solution still satisfies it (consistency). It differs from them in two ways:

* The added terms are not postulated. They are derived from a split of the solution into resolved and unresolved scales, the variational multiscale (VMS) method [#hughes1998]_. Besides the SUPG, PSPG and LSIC terms, this derivation yields a cross-stress term and a Reynolds-stress term, which model the effect of the unresolved (subgrid) scales on the resolved ones. No eddy viscosity is introduced: RBVMS acts as an implicit large-eddy simulation (LES) model of turbulence.
* The stabilization parameters are built from the element metric tensor instead of a single scalar element size. They therefore account for element stretching and shearing, which matters on the strongly anisotropic meshes used for wall-bounded turbulence.

The following subsections introduce the notation, explain the multiscale split and how it leads to the weak form, define the stabilization parameters and describe the two models of the fine-scale velocity available in Lethe.

Notation
~~~~~~~~

In addition to the index notation used in the rest of this page, the RBVMS terms are written with the element-wise inner products

.. math::

   \left(\mathbf{a}, \mathbf{b}\right) = \sum_K \int_{\Omega_K} a_k b_k \, \mathrm{d}\Omega_K, \qquad
   \left(\mathbf{A}, \mathbf{B}\right) = \sum_K \int_{\Omega_K} A_{kl} B_{kl} \, \mathrm{d}\Omega_K

for vectors :math:`\mathbf{a}, \mathbf{b}` and second-order tensors :math:`\mathbf{A}, \mathbf{B}`, and with the conventions

.. math::

   \left(\mathbf{a} \otimes \mathbf{b}\right)_{kl} = a_k b_l, \qquad
   \left(\nabla \mathbf{v}\right)_{kl} = \partial_l v_k, \qquad
   \left(\mathbf{a} \cdot \nabla \mathbf{v}\right)_k = a_l \partial_l v_k

The order of the factors of a tensor product matters, because the two products of the velocity :math:`\mathbf{u}` and of the fine-scale velocity :math:`\mathbf{u}'` introduced below are different terms:

.. math::

   \left(\nabla \mathbf{v}, \mathbf{u}' \otimes \mathbf{u}\right) = \sum_K \int_{\Omega_K} u'_k \, u_l \partial_l v_k \, \mathrm{d}\Omega_K = \left(\mathbf{u} \cdot \nabla \mathbf{v}, \mathbf{u}'\right), \qquad
   \left(\nabla \mathbf{v}, \mathbf{u} \otimes \mathbf{u}'\right) = \sum_K \int_{\Omega_K} u_k \, u'_l \partial_l v_k \, \mathrm{d}\Omega_K

In the first product, the velocity differentiates the test function along the streamlines: this is the SUPG term. The second product is the cross-stress term.

The Variational Multiscale Idea
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A finite element space can only represent the scales that are larger than the mesh size. The VMS method makes this explicit by splitting the exact velocity and pressure into a coarse-scale part, which belongs to the finite element space, and a fine-scale part, which does not:

.. math::

   \mathbf{u} = \mathbf{u}^h + \mathbf{u}', \qquad p = p^h + p'

The split is made unique by choosing a projector onto the finite element space, which defines what is meant by "coarse scales". The coarse scales play the role of the filtered field of classical LES, and the fine scales that of the subgrid field.

Substituting the split into the weak form and testing with the finite element functions :math:`(\mathbf{v}^h, q^h)` gives an equation for the coarse scales in which the fine scales still appear. Since :math:`\mathbf{u}'` and :math:`p'` cannot be represented on the mesh, they must be modelled. This is the closure problem, the same one that LES subgrid models address. The residual-based approach of Bazilevs *et al.* [#bazilevs2007]_ models the fine scales as proportional to the residuals of the coarse-scale solution:

.. math::

   \mathbf{u}' \approx -\tau_M \boldsymbol{\mathcal{R}}, \qquad p' \approx -\tau_C \, \partial_l u^h_l, \qquad
   \mathcal{R}_k = \partial_t u^h_k + u^h_l \partial_l u^h_k + \partial_k p^h - \nu \partial_l \partial_l u^h_k - f_k

where :math:`\boldsymbol{\mathcal{R}}` is the strong residual of the momentum equation and :math:`\partial_l u^h_l` that of the continuity equation, both evaluated element by element. The stabilization parameters :math:`\tau_M` and :math:`\tau_C` are defined below.

The idea behind this closure is simple. Where the mesh resolves the flow, the coarse-scale solution nearly satisfies the equations, the residuals are small, and so are the fine scales. Where it does not, the residuals, and therefore the modelled fine scales, grow. Consistency follows directly: the exact solution has zero residuals, hence zero fine scales, and all the added terms vanish.

From the Multiscale Split to the Weak Form
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since the exact velocity is divergence-free, the convective term can be written in conservative form, :math:`u_l \partial_l u_k = \partial_l \left(u_k u_l\right)`. Substitute the split, test with :math:`(\mathbf{v}^h, q^h)`, and integrate by parts the convective term and the fine-scale part of the continuity equation. The coarse-scale equations then contain the following terms:

* The convective term becomes :math:`-\left(\nabla \mathbf{v}^h, (\mathbf{u}^h + \mathbf{u}') \otimes (\mathbf{u}^h + \mathbf{u}')\right)`, which expands into four terms:

  - :math:`-\left(\nabla \mathbf{v}^h, \mathbf{u}^h \otimes \mathbf{u}^h\right)`, the Galerkin convective term;
  - :math:`-\left(\nabla \mathbf{v}^h, \mathbf{u}' \otimes \mathbf{u}^h\right) = -\left(\mathbf{u}^h \cdot \nabla \mathbf{v}^h, \mathbf{u}'\right)`, the SUPG term;
  - :math:`-\left(\nabla \mathbf{v}^h, \mathbf{u}^h \otimes \mathbf{u}'\right)`, the cross-stress term;
  - :math:`-\left(\nabla \mathbf{v}^h, \mathbf{u}' \otimes \mathbf{u}'\right)`, the Reynolds-stress term.

* The continuity equation gives :math:`\left(q^h, \partial_l u'_l\right) = -\left(\nabla q^h, \mathbf{u}'\right)`, the PSPG term.
* The pressure term gives :math:`-\left(\nabla \cdot \mathbf{v}^h, p'\right)`, the LSIC (grad-div) term.

The other fine-scale contributions are removed by the following assumptions of Bazilevs *et al.* [#bazilevs2007]_:

* The fine scales vanish on the boundary :math:`\Gamma` of the domain, so the integrations by parts above produce no boundary term.
* The viscous coupling between the scales, :math:`\nu \left(\nabla \mathbf{v}^h, \nabla \mathbf{u}'\right)`, is zero. This holds when the coarse scales are defined by the projector built from the viscous term (an :math:`\mathcal{H}^1_0`-type projection), for which the fine scales are orthogonal to the coarse ones. For this reason, unlike GLS, RBVMS has no :math:`-\nu \partial_l \partial_l v_k` term in the weight of the residual.
* The time derivative of the fine scales, :math:`\left(\mathbf{v}^h, \partial_t \mathbf{u}'\right)`, does not appear in the formulation. The consequences of neglecting the time dependence of the fine scales are discussed in the subsection on quasi-static and dynamic subscales below.

As in the other formulations of this page, Lethe writes the Galerkin terms with the advective form of the convective term, :math:`\left(\mathbf{v}^h, \mathbf{u}^h \cdot \nabla \mathbf{u}^h\right)`, and with the Laplacian form of the viscous term. Bazilevs *et al.* use the conservative form and the symmetric velocity gradient. These forms coincide for a divergence-free velocity, except for boundary terms on the boundaries where the velocity is not prescribed.

From here on, as in the rest of this page, :math:`\mathbf{u}` and :math:`p` denote the finite element solution and the superscript :math:`h` is dropped. Collecting the terms, the weak form assembled by Lethe is

.. math::

  & \left(q, \partial_l u_l\right) \underbrace{- \left(\nabla q, \mathbf{u}'\right)}_{\text{PSPG}} = 0
  \\
  & \left(\mathbf{v}, \partial_t \mathbf{u} + \mathbf{u} \cdot \nabla \mathbf{u} - \mathbf{f}\right) - \left(\nabla \cdot \mathbf{v}, p\right) + \nu \left(\nabla \mathbf{v}, \nabla \mathbf{u}\right)
  \\
  & \underbrace{- \left(\mathbf{u} \cdot \nabla \mathbf{v}, \mathbf{u}'\right)}_{\text{SUPG}}
    \underbrace{- \left(\nabla \cdot \mathbf{v}, p'\right)}_{\text{LSIC}}
    \underbrace{- \left(\nabla \mathbf{v}, \mathbf{u} \otimes \mathbf{u}'\right)}_{\text{cross stress}}
    \underbrace{- \left(\nabla \mathbf{v}, \mathbf{u}' \otimes \mathbf{u}'\right)}_{\text{Reynolds stress}} = 0

Lethe assembles the PSPG, SUPG, cross-stress and Reynolds-stress terms in this form, from the fine-scale velocity. Substituting :math:`\mathbf{u}' = -\tau_M \boldsymbol{\mathcal{R}}` and :math:`p' = -\tau_C \partial_l u_l` gives, in index notation,

.. math::

  &\int_{\Omega}  q  \partial_l u_l \mathrm{d}\Omega + \sum_{K} \int_{\Omega_K} \mathcal{R}_k \left(\tau_{M} \partial_k q \right) \mathrm{d}\Omega_K  = 0
  \\
  &\int_{\Omega}  v_k \left(\partial_t u_k+ u_l \partial_l u_k - f_k \right) \mathrm{d}\Omega - \int_{\Omega} \left( \partial_k v_k \right) p \mathrm{d}\Omega  + \nu \int_{\Omega} \left( \partial_l v_k \right) \left( \partial_l u_k  \right) \mathrm{d}\Omega
  \\
  & + \sum_{K} \int_{\Omega_K} \mathcal{R}_k \left(\tau_{M} u_l \partial_l v_k \right) \mathrm{d}\Omega_K + \sum_{K} \int_{\Omega_K} \tau_{C} \left( \partial_k v_k \right) \left( \partial_l u_l \right) \mathrm{d}\Omega_K
  \\
  & + \sum_{K} \int_{\Omega_K} \tau_{M} \, u_k \left( \partial_l v_k \right) \mathcal{R}_l \, \mathrm{d}\Omega_K - \sum_{K} \int_{\Omega_K} \tau_{M}^2 \left( \partial_l v_k \right) \mathcal{R}_k \mathcal{R}_l \, \mathrm{d}\Omega_K = 0

The first two lines are the continuity equation with the PSPG term and the Galerkin momentum equation. The third line adds the SUPG and LSIC terms, as in the GLS formulation. The last line contains the cross-stress and Reynolds-stress terms, which distinguish RBVMS from GLS. The Reynolds-stress term is quadratic in the residual. The table below summarizes the role of each added term.

.. list-table:: Terms added to the Galerkin formulation by RBVMS
   :header-rows: 1
   :widths: 16 30 54

   * - Term
     - Expression
     - Role
   * - PSPG
     - :math:`-\left(\nabla q, \mathbf{u}'\right)`
     - Fine-scale part of the continuity equation. Fills the zero pressure block of the Jacobian and stabilizes the pressure of equal-order elements.
   * - SUPG
     - :math:`-\left(\mathbf{u} \cdot \nabla \mathbf{v}, \mathbf{u}'\right)`
     - Transport of the fine-scale momentum by the coarse velocity. Adds dissipation along the streamlines and stabilizes convection-dominated flows.
   * - LSIC
     - :math:`-\left(\nabla \cdot \mathbf{v}, p'\right)`
     - Fine-scale pressure. Penalizes the divergence of the velocity and improves mass conservation.
   * - Cross stress
     - :math:`-\left(\nabla \mathbf{v}, \mathbf{u} \otimes \mathbf{u}'\right)`
     - Transport of the coarse momentum by the fine-scale velocity. Absent from the SUPG/PSPG and GLS formulations.
   * - Reynolds stress
     - :math:`-\left(\nabla \mathbf{v}, \mathbf{u}' \otimes \mathbf{u}'\right)`
     - Transport of the fine-scale momentum by the fine-scale velocity. Absent from the SUPG/PSPG and GLS formulations.

The names of the last two terms come from LES. Splitting the convective flux :math:`u_k u_l` into resolved and unresolved velocities produces the products of a resolved and an unresolved velocity, called cross stresses, and the product of two unresolved velocities, called the Reynolds (subgrid) stress. In the terminology of Bazilevs *et al.*, both mixed products, :math:`-\left(\nabla \mathbf{v}, \mathbf{u}' \otimes \mathbf{u}\right)` and :math:`-\left(\nabla \mathbf{v}, \mathbf{u} \otimes \mathbf{u}'\right)`, are cross-stress terms. Since the first one is the classical SUPG term, only the second one is usually called the cross-stress term. From this perspective, the SUPG/PSPG and GLS formulations are a first step towards the full variational multiscale formulation, which RBVMS completes with the two terms that model the subgrid stress.

Stabilization Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

The stabilization parameters are computed from the element contravariant metric tensor :math:`G_{ij}` and metric vector :math:`g_i`,

.. math::

   G_{ij} = \sum_k \frac{\partial \xi_k}{\partial x_i}\frac{\partial \xi_k}{\partial x_j}, \qquad  g_i = \sum_j \frac{\partial \xi_j}{\partial x_i}

where :math:`\boldsymbol{\xi}(\mathbf{x})` is the mapping from the physical element to the reference element, so that :math:`\partial \xi_k / \partial x_i` are the entries of the inverse of the mapping Jacobian. Bazilevs *et al.* define :math:`\boldsymbol{\xi}` on the bi-unit cell :math:`[-1,1]^d`.

The metric tensor measures the size of the element in every direction. On a rectangular cell of sides :math:`h_x \times h_y`, the mapping is :math:`\xi_i = 2 x_i / h_i + \text{constant}`, so that

.. math::

   \mathbf{G} = \begin{bmatrix} (2/h_x)^2 & 0 \\ 0 & (2/h_y)^2 \end{bmatrix}, \qquad
   u_i G_{ij} u_j = \left(\frac{2 u_x}{h_x}\right)^2 + \left(\frac{2 u_y}{h_y}\right)^2, \qquad
   G_{ij} G_{ij} = \left(\frac{2}{h_x}\right)^4 + \left(\frac{2}{h_y}\right)^4, \qquad
   \mathbf{g} = \left(\frac{2}{h_x}, \frac{2}{h_y}\right)

Each velocity component is therefore compared with the length of the element in its own direction. Consider a cell stretched along a wall (:math:`h_x \gg h_y`) with the flow along :math:`x`: the advective term is then :math:`(2|\mathbf{u}|/h_x)^2`, based on the length that the flow actually crosses. The diffusive term is dominated by the smallest size, :math:`h_y`, across which diffusion is fastest. A single element size, such as the diameter of the sphere of equal volume used by the SUPG/PSPG and GLS formulations, lies between :math:`h_x` and :math:`h_y` and cannot represent both directions. On sheared or curved cells, the off-diagonal entries of :math:`\mathbf{G}` account for the orientation of the element in the same way.

deal.II uses the reference cell :math:`[0,1]^d`, on which :math:`\partial \xi_k / \partial x_i` is half as large. Lethe therefore multiplies :math:`\partial \xi_k / \partial x_i` by :math:`2k`, where :math:`k` is the polynomial degree of the velocity interpolation. The factor 2 maps to the bi-unit cell. The factor :math:`k` makes the lengths :math:`h_x` and :math:`h_y` above the sides of the cell divided by :math:`k`, consistently with the element size of the other stabilizations, which is also divided by :math:`k`.

The momentum and continuity stabilization parameters are then

.. math::

   \tau_{M} = \left( \frac{4}{\Delta t^{2}} + u_i G_{ij} u_j + C_I \nu^2 G_{ij} G_{ij} \right)^{-1/2}, \qquad \tau_{C} = \frac{1}{\tau_{M}\, g_i g_i}

where :math:`C_I` is a positive constant, independent of the mesh size. It comes from an element-wise inverse estimate, an inequality that bounds the second derivatives of a finite element function on an element by its first derivatives divided by the element size. Lethe uses :math:`C_I = 9`: since the metric tensor already includes the degree scaling, this value recovers, in 1D, the diffusion-dominated limit :math:`\tau = h^2/(12\nu)` of the SUPG/PSPG formulation. Both the SUPG and the PSPG terms use :math:`\tau_M`, as there is no separate pressure stabilization parameter. In steady simulations, the transient :math:`4/\Delta t^{2}` term is dropped.

The three terms under the square root of :math:`\tau_M` are squared rates, that is squared inverse time scales: :math:`2/\Delta t` for the time step, :math:`(u_i G_{ij} u_j)^{1/2}`, of the order of :math:`|\mathbf{u}|/h`, for the advection across the element, and :math:`\sqrt{C_I}\, \nu\, (G_{ij} G_{ij})^{1/2}`, of the order of :math:`\nu/h^2`, for the diffusion across it. The largest rate, that is the shortest time scale, dominates. :math:`\tau_M` can thus be read as the time during which a fine-scale perturbation lives in an element, before it is advected out of it, diffused or overtaken by the next time step. The continuity parameter :math:`\tau_C` has the units of a kinematic viscosity. Since :math:`g_i g_i` is of the order of :math:`1/h^2`, :math:`\tau_C` is of the order of :math:`h^2/\tau_M`, that is :math:`|\mathbf{u}| h` in advection-dominated flows. Bazilevs *et al.* derive it from the fine-scale Schur complement of the pressure.

The table below compares the limits of the parameters on a Cartesian cell in :math:`d` dimensions with those of the SUPG/PSPG and GLS formulations.

.. list-table:: Limits of the stabilization parameters on a Cartesian cell in :math:`d` dimensions
   :header-rows: 1
   :widths: 40 25 35

   * - Regime
     - SUPG/PSPG and GLS
     - RBVMS (quasi-static subscale)
   * - Momentum parameter, small time step
     - :math:`\tau \rightarrow \Delta t`
     - :math:`\tau_M \rightarrow \Delta t / 2`
   * - Momentum parameter, advection-dominated
     - :math:`\tau \rightarrow h / (2|\mathbf{u}|)`
     - :math:`\tau_M \rightarrow h / (2|\mathbf{u}|)`
   * - Momentum parameter, diffusion-dominated
     - :math:`\tau \rightarrow h^2 / (12\nu)`
     - :math:`\tau_M \rightarrow h^2 / (12\sqrt{d}\,\nu)`
   * - Continuity parameter, advection-dominated
     - :math:`\tau_{LSIC} = |\mathbf{u}| h / 2` (GLS only)
     - :math:`\tau_C \rightarrow |\mathbf{u}| h / (2d)`
   * - Continuity parameter, small time step
     - :math:`\tau_{LSIC}` does not depend on :math:`\Delta t`
     - :math:`\tau_C \rightarrow h^2 / (2d\, \Delta t)`, which grows without bound

The expressions have the same form in both columns, but :math:`h` is the length scale of each formulation: the side of the cell divided by :math:`k` for RBVMS, and the diameter of the sphere of equal volume divided by :math:`k` for SUPG/PSPG and GLS. On a square or cubic cell, the latter is about 13% larger in 2D and 24% larger in 3D, so the values coincide only in 1D. In the advection-dominated limit, :math:`\tau_C` matches the LSIC parameter of the GLS formulation in 1D only, and is :math:`d` times smaller in :math:`d` dimensions. In the diffusion-dominated limit, :math:`\tau_M` is :math:`\sqrt{d}` times smaller than the SUPG/PSPG parameter.

Solving the Non-linear Problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As in the SUPG/PSPG and GLS formulations, the Jacobian of the Newton-Raphson method is assembled with the stabilization parameters :math:`\tau_M` and :math:`\tau_C` frozen: their dependence on the velocity is not linearized. All the other dependencies on the solution are linearized. Since the fine-scale velocity is affine in the residual, its increment is :math:`\delta \mathbf{u}' = -\tau_M\, \delta \boldsymbol{\mathcal{R}}`, where :math:`\delta \boldsymbol{\mathcal{R}}` is the linearization of the strong residual. The product rule is then applied to the SUPG weight :math:`\mathbf{u} \cdot \nabla \mathbf{v}` and to the cross-stress and Reynolds-stress terms, which are products of :math:`\mathbf{u}` and :math:`\mathbf{u}'`. The Jacobian has the same block structure as that of the SUPG/PSPG formulation, including the additional block :math:`S^*` produced by the PSPG term.

Fine-Scale Velocity: Quasi-Static and Dynamic Subscales
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The closure :math:`\mathbf{u}' = -\tau_M \boldsymbol{\mathcal{R}}` can be understood as the solution of a model equation for the fine scales. Following Codina *et al.* [#codina2007]_, the fine-scale momentum equation is approximated, in each element, by

.. math::

   \partial_t \mathbf{u}' + \frac{1}{\tau_s}\mathbf{u}' = -\boldsymbol{\mathcal{R}}, \qquad \tau_{s} = \left( u_i G_{ij} u_j + C_I \nu^2 G_{ij} G_{ij} \right)^{-1/2}

where the advection and diffusion of the fine scales inside the element are replaced by the algebraic term :math:`\mathbf{u}'/\tau_s`. The steady stabilization parameter :math:`\tau_s` is :math:`\tau_M` without the :math:`4/\Delta t^2` term. This equation describes a relaxation: the coarse-scale residual drives the fine scales, which relax towards :math:`-\tau_s \boldsymbol{\mathcal{R}}` over the time scale :math:`\tau_s`. The ``velocity subscale`` parameter of the ``stabilization`` subsection selects how the time derivative :math:`\partial_t \mathbf{u}'` is treated.

**Quasi-static subscale.** With the default ``quasi_static`` subscale, :math:`\partial_t \mathbf{u}'` is not tracked, and :math:`\mathbf{u}' = -\tau_M \boldsymbol{\mathcal{R}}` is computed from the current residual only. The time step enters through the :math:`4/\Delta t^2` term of :math:`\tau_M`. Since :math:`\tau_M^{-2} = (2/\Delta t)^2 + \tau_s^{-2}`, the rate of the time step is combined with the steady rate :math:`1/\tau_s`, much as :math:`\beta_0` is in the parameter :math:`\tau_t` of the dynamic subscale below, but without the memory of the previous fine-scale velocities. This has two consequences:

* A steady state reached by time stepping depends on the time step, since :math:`\tau_M` does.
* When :math:`\Delta t` becomes small, :math:`\tau_M \rightarrow \Delta t / 2`, so the terms that use :math:`\tau_M` vanish, in particular the PSPG term that stabilizes the pressure of equal-order elements. At the same time, :math:`\tau_C = 1/(\tau_M\, g_i g_i)` grows without bound. Bazilevs *et al.* point out this degeneracy and refer to the dynamic subscales of Codina *et al.* to address it.

**Dynamic subscale.** With the ``dynamic`` subscale, following Codina *et al.* [#codina2007]_, the evolution equation above is solved for :math:`\mathbf{u}'`. It is discretized with the same BDF scheme as the velocity, :math:`\partial_t \mathbf{u}' \approx \beta_0 \mathbf{u}'_{n+1} + \sum_{j \geq 1} \beta_j \mathbf{u}'_{n+1-j}`, which gives

.. math::

   \mathbf{u}'_{n+1} = -\tau_t \left( \boldsymbol{\mathcal{R}} + \mathbf{m}' \right), \qquad \tau_t = \left( \beta_0 + \frac{1}{\tau_s} \right)^{-1}, \qquad \mathbf{m}' = \sum_{j \geq 1} \beta_j \mathbf{u}'_{n+1-j}

The memory :math:`\mathbf{m}'` requires the fine-scale velocity of the previous time steps, which Lethe stores at the quadrature points of the cells. The fine-scale velocity thus becomes a history variable, like the internal variables of inelastic constitutive laws in solid mechanics. The continuity stabilization parameter is built from :math:`\tau_s`, :math:`\tau_C = 1/(\tau_s\, g_i g_i)`, so it does not depend on the time step either.

The memory changes the behaviour of the stabilization in both regimes:

* At a steady state, :math:`\mathbf{u}'` no longer changes. Since the BDF coefficients sum to zero, :math:`\beta_0 \mathbf{u}' + \mathbf{m}' = 0`, so that :math:`\mathbf{u}' = -\tau_s \boldsymbol{\mathcal{R}}`: the steady solution no longer depends on the time step.
* When the time step is small, :math:`\tau_t` is also of the order of :math:`\Delta t`, but the memory compensates for it. For BDF1, :math:`\mathbf{u}'_{n+1} \approx \mathbf{u}'_n - \Delta t\, \boldsymbol{\mathcal{R}}`, so :math:`\mathbf{u}'` keeps its magnitude and relaxes towards :math:`-\tau_s \boldsymbol{\mathcal{R}}` over the time scale :math:`\tau_s`.

.. list-table:: Comparison of the two models of the fine-scale velocity
   :header-rows: 1
   :widths: 30 35 35

   * -
     - ``quasi_static``
     - ``dynamic``
   * - Fine-scale velocity
     - :math:`\mathbf{u}' = -\tau_M \boldsymbol{\mathcal{R}}`
     - :math:`\mathbf{u}'_{n+1} = -\tau_t \left( \boldsymbol{\mathcal{R}} + \mathbf{m}' \right)`
   * - Momentum parameter
     - :math:`\tau_M`, which contains :math:`4/\Delta t^2`
     - :math:`\tau_t = \left( \beta_0 + 1/\tau_s \right)^{-1}`
   * - Continuity parameter
     - :math:`\tau_C = 1/(\tau_M\, g_i g_i)`
     - :math:`\tau_C = 1/(\tau_s\, g_i g_i)`
   * - Steady state reached by time stepping
     - Depends on :math:`\Delta t`
     - Independent of :math:`\Delta t`
   * - Small time step
     - :math:`\mathbf{u}' \rightarrow 0` and :math:`\tau_C \rightarrow \infty`
     - :math:`\mathbf{u}'` keeps its magnitude
   * - Storage
     - None
     - :math:`\mathbf{u}'` of the previous time steps at each quadrature point

In steady simulations, there is no time derivative and the two models are identical.

Implementation Notes and Limitations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* The RBVMS stabilization is currently available for constant-density (incompressible) Newtonian flows with the monolithic solvers (``lethe-fluid`` and ``lethe-fluid-matrix-free``).
* The dynamic subscale is currently available with the matrix-free solver (``lethe-fluid-matrix-free``), for steady simulations and BDF time-stepping methods. It supports neither mesh adaptation in transient simulations nor the mortar method.
* With the dynamic subscale, the fine-scale velocity starts from zero, and the term :math:`\left(\mathbf{v}, \partial_t \mathbf{u}'\right)` is neglected in the momentum equation. The Jacobian is assembled with :math:`\tau_t` and :math:`\mathbf{m}'` frozen, so that :math:`\delta \mathbf{u}' = -\tau_t\, \delta \boldsymbol{\mathcal{R}}`. On the levels of the geometric multigrid preconditioner, the memory is neglected, which only affects the preconditioner.
* The history of the fine-scale velocity at the quadrature points is attached to the cells of the triangulation when a checkpoint is written, so that a restarted simulation continues as an uninterrupted one would. The conditions for restarting are given in the :doc:`stabilization parameters <../../../parameters/cfd/stabilization>` page.


----------
References
----------

.. [#bazilevs2007] \Y. Bazilevs, V. M. Calo, J. A. Cottrell, T. J. R. Hughes, A. Reali, and G. Scovazzi, “Variational multiscale residual-based turbulence modeling for large eddy simulation of incompressible flows,” *Comput. Methods Appl. Mech. Eng.*, vol. 197, no. 1–4, pp. 173–201, 2007, doi: `10.1016/j.cma.2007.07.016 <https://doi.org/10.1016/j.cma.2007.07.016>`_\.

.. [#hughes1998] \T. J. R. Hughes, G. R. Feijóo, L. Mazzei, and J.-B. Quincy, “The variational multiscale method—a paradigm for computational mechanics,” *Comput. Methods Appl. Mech. Eng.*, vol. 166, no. 1–2, pp. 3–24, 1998, doi: `10.1016/S0045-7825(98)00079-6 <https://doi.org/10.1016/S0045-7825(98)00079-6>`_\.

.. [#codina2007] \R. Codina, J. Principe, O. Guasch, and S. Badia, “Time dependent subscales in the stabilized finite element approximation of incompressible flow problems,” *Comput. Methods Appl. Mech. Eng.*, vol. 196, no. 21–24, pp. 2413–2430, 2007, doi: `10.1016/j.cma.2007.01.002 <https://doi.org/10.1016/j.cma.2007.01.002>`_\.
