===========================
Resolved CFD-DEM
===========================

In resolved CFD-DEM, the incompressible Navier-Stokes equations are solved on a mesh that is significantly finer than the particle size. The particle-fluid coupling is then obtained directly by imposing a no-slip boundary condition on the surface of the moving particles [#berard2020]_.

Sharp-interface immersed boundary method
----------------------------------------

Lethe uses a Sharp-Interface Immersed Boundary Method (SIBM) to impose the no-slip boundary condition on moving solids while coupling the solid dynamics and fluid flow implicitly. This implicit formulation keeps the resolved CFD-DEM solver stable even for large fluid-to-solid density ratios [#hu1992]_.

The SIBM retains the order of convergence of the underlying finite element discretization (up to fourth order) for both stationary and moving boundaries [#barbeau2022]_ [#daunais2023]_ and its full algorithmic details are presented in [#barbeau2024]_.

Hydrodynamic forces and lubrication correction
----------------------------------------------

The hydrodynamic force :math:`\mathbf{F}^\mathrm{pf}` and torque :math:`\mathbf{M}^\mathrm{pf}` acting on each particle are evaluated by integrating the fluid stress tensor over the particle surface (see the particle momentum balance in :doc:`../dem/dem`). When two particles approach, the finite grid resolution can under-resolve the lubrication force, so Lethe adds the correction proposed by [#hori2022]_:

.. math::
    \mathbf{F}^\mathrm{l_c}_{ij} &= \mathbf{F}^\mathrm{l}_{ij} - \mathbf{F}^\mathrm{l}_{ij} \big|_{l=\varepsilon_0} \\
    \mathbf{F}^\mathrm{l}_{ij} &= \frac{3}{2} \pi \mu_f \left (\frac{d_{p_i} d_{p_j}}{d_{p_i} + d_{p_j}} \right )^2 \frac{1}{l} \left (\mathbf{v}_{ij} \cdot \mathbf{e}_{ij} \right ) \mathbf{e}_{ij}

where:

* :math:`\mu_f` is the fluid dynamic viscosity;
* :math:`d_{p}` is the particle diameter;
* :math:`l` is the gap between the particles;
* :math:`\varepsilon_0` is a minimum gap used to regularize the lubrication force;
* :math:`\mathbf{v}_{ij}` is the relative velocity between particles :math:`i` and :math:`j`;
* :math:`\mathbf{e}_{ij}` is the unit vector pointing from particle :math:`i` to :math:`j`.

The model is derived for spheres; it can be enabled for non-spherical solids, but results should be interpreted with care.

Signed distance functions for complex solids
--------------------------------------------

The SIBM supports multiple solid descriptions (spheres, cylinders, boxes, CAD geometries in ``step``/``iges``/``stl`` formats, and Radial Basis Function (RBF) surfaces) [#guevremont2024]_. Each solid is represented through a Signed Distance Function (SDF) :math:`\lambda`. For every node adjacent to a cut cell where :math:`\lambda(\mathbf{x})=0`, the method extrapolates along the outward normal

.. math::
    \mathbf{n}(\mathbf{x}) = \nabla \lambda (\mathbf{x})

using Lagrange polynomials. The extrapolated values are used both to enforce the no-slip boundary condition at the moving interface and to evaluate stresses on the particle surface.

.. [#berard2020] \A. Bérard, G. S. Patience, and B. Blais, “Experimental methods in chemical engineering: Unresolved CFD-DEM,” *Can. J. Chem. Eng.*, vol. 98, no. 2, pp. 424–440, 2020, doi: `10.1002/cjce.23686 <https://doi.org/10.1002/cjce.23686>`_\.
.. [#hu1992] \H. H. Hu, N. A. Patankar, and M. Y. Zhu, “Direct numerical simulations of fluid–solid systems using the distributed Lagrange multiplier/fictitious domain method,” *Theor. Comput. Fluid Dyn.*, vol. 3, no. 3, pp. 285–306, 1992.
.. [#barbeau2022] \L. Barbeau, S. Étienne, C. Béguin, and B. Blais, “Development of a high-order continuous Galerkin sharp-interface immersed boundary method and its application to incompressible flow problems,” *Comput. Fluids*, vol. 239, p. 105415, May 2022, doi: `10.1016/j.compfluid.2022.105415 <https://doi.org/10.1016/j.compfluid.2022.105415>`_\.
.. [#daunais2023] \Y. Daunais, L. Barbeau, C. Béguin, and B. Blais, “Extensive verification of a high-order sharp-interface immersed boundary method for incompressible flows,” *Int. J. Numer. Meth. Fluids*, vol. 95, no. 2, pp. 100–123, 2023.
.. [#barbeau2024] \L. Barbeau, Y. Daunais, C. Béguin, and B. Blais, “An implicit sharp-interface immersed boundary method for fluid–solid interaction with moving boundaries,” 2024.
.. [#hori2022] \N. Hori, M. E. Rosti, and S. Takagi, “An Eulerian-based immersed boundary method for particle suspensions with implicit lubrication model,” *Comput. Fluids*, vol. 236, p. 105278, 2022.
.. [#guevremont2024] \V. Guévremont, L. Barbeau, C. Béguin, and B. Blais, “Pore-resolved CFD digital twin of additive manufactured heat exchangers using radial basis function geometry descriptions,” 2024.
