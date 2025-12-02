===========================
Resolved CFD-DEM
===========================

In resolved CFD-DEM, the incompressible Navier-Stokes equations are solved on a mesh that is significantly finer than the particle size. The particle-fluid coupling is then obtained directly by imposing a no-slip boundary condition on the surface of the moving particles [#berard2020]_.

Sharp-interface immersed boundary method
----------------------------------------

Lethe uses a Sharp-Interface Immersed Boundary Method (SIBM) to impose the no-slip boundary condition on moving solids while coupling the solid dynamics and fluid flow implicitly. This implicit formulation keeps the resolved CFD-DEM solver stable even for large fluid-to-solid density ratios.

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

The SIBM supports multiple solid descriptions: spheres, cylinders, boxes, CAD geometries in ``step``/``iges``/``stl`` formats, Radial Basis Function (RBF) surfaces [#guevremont2024]_, and composites defined using boolean operations. Each solid is represented through a Signed Distance Function (SDF) :math:`\lambda`. For every node adjacent to a cut cell where :math:`\lambda(\mathbf{x})=0`, the method extrapolates along the outward normal

.. math::
    \mathbf{n}(\mathbf{x}) = \nabla \lambda (\mathbf{x})

using Lagrange polynomials. The extrapolated values are used both to enforce the no-slip boundary condition at the moving interface and to evaluate stresses on the particle surface.

.. [#berard2020] \A. Bérard, G. S. Patience, and B. Blais, “Experimental methods in chemical engineering: Unresolved CFD-DEM,” *Canadian Journal of Chemical Engineering*, vol. 98, no. 2, pp. 424–440, 2020, doi: `10.1002/cjce.23686 <https://doi.org/10.1002/cjce.23686>`_\.
.. [#barbeau2022] \L. Barbeau, S. Étienne, C. Béguin, and B. Blais, “Development of a high-order continuous Galerkin sharp-interface immersed boundary method and its application to incompressible flow problems,” *Computer & Fluids*, vol. 239, p. 105415, May 2022, doi: `10.1016/j.compfluid.2022.105415 <https://doi.org/10.1016/j.compfluid.2022.105415>`_\.
.. [#daunais2023] C.-A. Daunais, L. Barbeau, and B. Blais, "An extensive study of shear thinning flow around a spherical particle for power-law and Carreau fluids," *Journal of Non-Newtonian Fluid Mechanics*, vol. 311, p. 104951, 2023. doi: `10.1016/j.jnnfm.2022.104951 <https://doi.org/10.1016/j.jnnfm.2022.104951>`_
.. [#barbeau2024] C. Barbeau, S. Golshan, J. Deng, S. Étienne, C. Béguin, and B. Blais, "High-order moving immersed boundary and its application to a resolved CFD-DEM model," *Computer & Fluids*, vol. 268, p. 106094, 2024. doi: `10.1016/j.compfluid.2023.106094 <https://doi.org/10.1016/j.compfluid.2023.106094>`_
.. [#hori2022] N. Hori, M. E. Rosti, and S. Takagi, "An Eulerian-based immersed boundary method for particle suspensions with implicit lubrication model," *Computer & Fluids*, vol. 236, p. 105278, 2022. doi: `10.1016/j.compfluid.2021.105278 <https://doi.org/10.1016/j.compfluid.2021.105278>`_
.. [#guevremont2024] O. Guévremont, L. Barbeau, V. Moreau, F. Galli, N. Virgilio, and B. Blais, "Robust pore-resolved CFD through porous monoliths reconstructed by micro-computed tomography: From digitization to flow prediction," *Chemical Engineering Journal*, vol. 504, p. 158577, 2025. doi: `10.1016/j.cej.2024.158577 <https://doi.org/10.1016/j.cej.2024.158577>`_

