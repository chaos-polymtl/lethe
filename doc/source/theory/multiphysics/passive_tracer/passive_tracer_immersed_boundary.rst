=========================================
Boundary Conditions on Immersed Particles
=========================================

With the sharp-interface immersed boundary method of ``lethe-fluid-sharp`` (see :doc:`../../multiphase/cfd_dem/resolved_cfd-dem`), a boundary condition can be imposed on the tracer at the surface :math:`\Gamma` of the particles. A Dirichlet boundary condition

.. math::
    C = C_\Gamma \quad \text{on } \Gamma

models, for example, a surface reaction that is much faster than the transport of the tracer towards the particle. The reaction is then limited by the diffusion of the tracer, and :math:`C_\Gamma = 0`. A Neumann or a Robin boundary condition instead imposes the flux of tracer into the particle (see `Flux Boundary Conditions`_), for example the rate of a slower surface reaction.

The boundary condition is imposed as the no-slip condition of the fluid [#barbeau2022]_. Since the mesh does not conform to the particles, each cell is classified with the signed distance to the particles at the support points of its degrees of freedom:

* a fluid cell has all its support points outside the particles;
* a cell inside a particle has all its support points inside the particles;
* the other cells are cut by the surface of a particle.

The weak form of the tracer equation is only assembled on the fluid cells. The degrees of freedom shared by a fluid cell and a cut cell therefore form the boundary of the computational domain, and the equation of each of them is replaced by a stencil that imposes the boundary condition.

Let :math:`\mathbf{x}_i` be such a degree of freedom, :math:`\mathbf{x}_\Gamma` the point of the surface closest to it and :math:`d = \| \mathbf{x}_i - \mathbf{x}_\Gamma \|` its distance to the surface. The stencil uses the polynomial :math:`C_h^S` of a fluid cell :math:`S` that contains :math:`\mathbf{x}_i` and lies on the side of the fluid. It evaluates :math:`C_h^S` at :math:`p+1` points of the normal to the surface that passes through :math:`\mathbf{x}_i`: the degree of freedom itself and :math:`p` points further from the surface, the farthest one at a distance :math:`d/L` from the degree of freedom. The polynomial of degree :math:`p` that interpolates these values is extrapolated to :math:`\mathbf{x}_\Gamma`, where it must be equal to the boundary value:

.. math::
    \sum_{k=0}^{p} a_k \: C_h^S(\mathbf{x}_k) = C_\Gamma(\mathbf{x}_\Gamma),

where :math:`a_k` are the values of the Lagrange basis at the surface, :math:`p` is the ``stencil degree`` and :math:`L` is the ``length ratio`` of the subsection ``extrapolation function`` (see :doc:`../../../parameters/sharp-immersed-boundary/sharp-immersed-boundary`). Since :math:`C_h^S` is a linear combination of the degrees of freedom of :math:`S`, this equation is a row of the matrix. When no such stencil exists, for example when the cell :math:`S` is itself cut by a particle, the boundary value is imposed directly on the degree of freedom: :math:`C_i = C_\Gamma(\mathbf{x}_\Gamma)`.

The tracer is not solved inside the particles. The degrees of freedom that only belong to cut cells and to cells inside the particles take the harmonic extension of the solution: a Laplace equation is assembled on these cells. The equations of the other degrees of freedom do not involve them, so this extension does not modify the solution in the fluid. It however gives meaningful values to the degrees of freedom that a moving particle uncovers, since these values enter the time derivative of the following time steps.

Flux Boundary Conditions
------------------------

The Neumann and Robin boundary conditions impose the flux of tracer into the particle,

.. math::
    J = D \frac{\partial C}{\partial n} \quad \text{on } \Gamma,

where :math:`D` is the tracer diffusivity and :math:`\mathbf{n}` is the normal pointing out of the particle. A positive flux removes tracer from the fluid. The Neumann boundary condition imposes :math:`J = q`, and the Robin boundary condition imposes the rate of a surface reaction of order :math:`n`,

.. math::
    J = \alpha \: C_\text{eff}^{n-1} \: C + q \quad \text{with} \quad C_\text{eff} = \sqrt{C^2 + \varepsilon^2},

where :math:`\alpha` is the reaction constant, :math:`\varepsilon` a threshold that keeps the reaction rate differentiable when the order is lower than 1, and :math:`q` an additional flux.

The normal derivative at :math:`\mathbf{x}_\Gamma` is the derivative of the same extrapolated polynomial. On the stencil, the coordinate :math:`s` increases towards the surface, and a unit of :math:`s` is a distance :math:`d/L`. Therefore,

.. math::
    \frac{\partial C}{\partial n}(\mathbf{x}_\Gamma) = -\frac{L}{d} \sum_{k=0}^{p} b_k \: C_h^S(\mathbf{x}_k),

where :math:`b_k` are the derivatives of the Lagrange basis at the surface with respect to :math:`s`. The equation of the degree of freedom becomes

.. math::
    -D \frac{L}{d} \sum_{k=0}^{p} b_k \: C_h^S(\mathbf{x}_k) - J\left(\sum_{k=0}^{p} a_k \: C_h^S(\mathbf{x}_k), \mathbf{x}_\Gamma\right) = 0.

It is nonlinear for a Robin boundary condition of order :math:`n \neq 1`, and its derivative with respect to the degrees of freedom of :math:`S` is assembled in the Jacobian matrix of Newton's method. When no stencil exists, the flux is computed with the value and the gradient of the solution at the degree of freedom, in a fluid cell :math:`F` that contains it:

.. math::
    D \: \nabla C_h^F(\mathbf{x}_i) \cdot \mathbf{n} - J(C_i, \mathbf{x}_\Gamma) = 0.

This approximation is only first-order accurate. The number of degrees of freedom without stencil is printed when the ``verbosity`` of the tracer non-linear solver is ``verbose``.

Flux into the Particles
-----------------------

The flux of tracer into a particle :math:`P`, :math:`\Phi_P = \int_{\Gamma_P} J \: \mathrm{d}\Gamma`, is the rate at which the particle consumes the tracer, for example the rate of a surface reaction. When ``calculate tracer flux on particles`` is enabled in the subsection ``post-processing`` (see :doc:`../../../parameters/cfd/post_processing`), it is computed from the residual of the tracer equation on the fluid cells rather than from the gradient of the solution. This consistent flux is more accurate than the gradient, which converges one order slower than the solution.

Let :math:`\mathcal{I}_P` be the degrees of freedom whose equation is replaced by a stencil and whose closest particle is :math:`P`, :math:`\varphi_i` their test functions, and :math:`w_P = \sum_{i \in \mathcal{I}_P} \varphi_i`. The function :math:`w_P` is one on the boundary :math:`\partial \Omega_{h,P}` of the fluid cells around :math:`P`, and it is zero beyond the first layer of fluid cells. Let :math:`R_i(C_h)` be the residual of the weak form of the tracer equation for the test function :math:`\varphi_i`, assembled on the fluid cells :math:`\Omega_h`. Integrating its diffusion term by parts gives

.. math::
    \sum_{i \in \mathcal{I}_P} R_i(C_h) = \int_{\Omega_h} r(C_h) \: w_P \: \mathrm{d}\Omega + \int_{\partial \Omega_{h,P}} D \: \nabla C_h \cdot \mathbf{n}_h \: \mathrm{d}\Gamma,

where :math:`r` is the strong residual of the tracer equation and :math:`\mathbf{n}_h` is the normal pointing out of the fluid cells, towards the particle. The last term is the opposite of the diffusive flux into the particle through :math:`\partial \Omega_{h,P}`, so the flux into the particle is computed as

.. math::
    \Phi_P = -\sum_{i \in \mathcal{I}_P} R_i(C_h),

with the residuals evaluated at the solution of the time step, before the equations of these degrees of freedom are replaced by the stencils. This flux crosses :math:`\partial \Omega_{h,P}` instead of the surface of the particle. The two differ by the accumulation, advection, reaction and sources of tracer in the cut cells between them, which vanish for a steady diffusion without sources, and whose volume vanishes as the mesh is refined. The degrees of freedom of :math:`\mathcal{I}_P` constrained by a boundary condition of the domain are not counted, so the flux into a particle that touches the boundary of the domain is underestimated.

In a transient simulation, the flux is the balance of the time-stepping scheme over the fluid cells next to the particle. After a discontinuous initial condition, for example a particle that starts consuming the tracer at :math:`t = 0`, the flux of the first time steps can therefore oscillate with BDF2, and even become negative, while their sum remains consistent. For a moving particle, the degrees of freedom that change category from one time step to the next also make the flux of a single time step noisy: it should be averaged over several time steps.

.. [#barbeau2022] \L. Barbeau, S. Étienne, C. Béguin, and B. Blais, “Development of a high-order continuous Galerkin sharp-interface immersed boundary method and its application to incompressible flow problems,” *Computer & Fluids*, vol. 239, p. 105415, May 2022, doi: `10.1016/j.compfluid.2022.105415 <https://doi.org/10.1016/j.compfluid.2022.105415>`_\.
