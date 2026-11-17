=========================================
Boundary Conditions on Immersed Particles
=========================================

With the sharp-interface immersed boundary method of ``lethe-fluid-sharp`` (see :doc:`../../multiphase/cfd_dem/resolved_cfd-dem`), a boundary condition can be imposed on the tracer at the surface :math:`\Gamma` of the particles. A Dirichlet boundary condition

.. math::
    C = C_\Gamma \quad \text{on } \Gamma

models, for example, a surface reaction that is much faster than the transport of the tracer towards the particle. The reaction is then limited by the diffusion of the tracer, and :math:`C_\Gamma = 0`.

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

.. [#barbeau2022] \L. Barbeau, S. Étienne, C. Béguin, and B. Blais, “Development of a high-order continuous Galerkin sharp-interface immersed boundary method and its application to incompressible flow problems,” *Computer & Fluids*, vol. 239, p. 105415, May 2022, doi: `10.1016/j.compfluid.2022.105415 <https://doi.org/10.1016/j.compfluid.2022.105415>`_\.
