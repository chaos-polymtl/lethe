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

Metric-tensor definition
~~~~~~~~~~~~~~~~~~~~~~~~~~

The definition above relies on a single, isotropic element size :math:`h`, which discards the directional information of the element and becomes a poor estimate for deformed, stretched or high-aspect-ratio cells. As an alternative, Lethe can compute :math:`\tau` from the element covariant metric tensor following the variational multiscale (VMS) approach of `Bazilevs et al. (2007) <https://doi.org/10.1016/j.cma.2007.07.016>`_. The covariant metric tensor is defined as

.. math::

   G_{ij} = \sum_{k} \frac{\partial \xi_k}{\partial x_i} \frac{\partial \xi_k}{\partial x_j} = J^{-T} J^{-1}

where :math:`\xi` are the reference coordinates, :math:`x` the physical coordinates and :math:`J` the mapping Jacobian. The momentum stabilization parameter is then

.. math::

   \tau = \left[ \left( \frac{1}{\Delta t} \right)^{2} + u_i G_{ij} u_j + C_I \nu^2 G_{ij} G_{ij} \right]^{-1/2}

where :math:`C_I` is a positive constant stemming from an inverse estimate (:math:`C_I = 36` for linear elements), and the grad-div (LSIC) parameter is :math:`\tau_{LSIC} = \left( \tau\, g_i g_i \right)^{-1}` with :math:`g_i = \sum_j \partial \xi_j / \partial x_i`. To remain consistent with the polynomial degree :math:`p` (which divides the element size in the isotropic definition), the metric tensor and the vector :math:`g` are scaled so that :math:`G` carries a factor :math:`p^2` and :math:`g` a factor :math:`p`, making the definition robust for higher-order elements. Because the metric tensor carries the anisotropic size information of the element, this definition is more robust on deformed meshes. It is selected with ``set stabilization parameter = metric_tensor`` and is currently only supported by the matrix-free solver.

.. note::
   The anisotropic metric-tensor definition strongly reduces :math:`\tau` on high-aspect-ratio cells (e.g. boundary-layer cells). Since the same parameter also multiplies the PSPG term, which is responsible for the pressure stabilization of equal-order elements, using it for PSPG would under-stabilize the pressure and can cause the linear solver to diverge. For this reason, the ``metric_tensor`` definition is applied only to the SUPG and LSIC terms, while the PSPG term always uses the isotropic (Tezduyar) definition.

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

