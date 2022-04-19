Unresolved CFD-DEM coupling
############################

Unresolved CFD-DEM coupling is a technique with high potential for design and analysis of multiphase flows involving particles and fluid, such as fluidized beds, stirred and floculation tanks. In this approach, we apply Newton's second law to each particle individually such that their movement is described at a micro scale (as in DEM simulations). For the fluid, the Volume Average Navier-Stokes (VANS) equations are used to describe the fluid at a meso scale. The meso-micro scale allows for particle-fluid simulations involving large numbers of particles with reasonable computational cost and highly detailed results (in both time and space).

In this guide, we summarize the theory behind Unresolved CFD-DEM. For further details, we refer the reader to the articles by Bérard, Patience & Blais (2019) `[1] <https://doi.org/10.1002/cjce.23686>`_, and Zhou et al. (2010) `[2] <https://doi.org/10.1017/S002211201000306X>`_.

Applying Newton's second law on the particle :math:`i` surrounded by fluid a :math:`f`, we find:

.. math::
    m_i \frac{d \mathbf{v}_i}{dt} = \sum_{j}\mathbf{f}_{c,ij} + \sum_{j}\mathbf{f}_{nc,ij} + \mathbf{f}_{pf,i} + \mathbf{f}_{g,i}
    
.. math::    
    I_i \frac{d\mathbf{\omega}_i}{dt} = \sum_{j}\left ( \mathbf{M}_{c,ij} + \mathbf{M}_{r,ij} \right ) + \sum_{w}\left ( \mathbf{M}_{c,iw} + \mathbf{M}_{r,iw} \right )

where:

* :math:`m_i` is the mass of the particle :math:`i`;
* :math:`v_i` is the velocity vector of the particle :math:`i`;
* :math:`\mathbf{f}_{c,ij}` are the contact forces between particles :math:`i` and :math:`j` (detailed in the DEM section of this guide);
* :math:`\mathbf{f}_{nc,ij}` are the non-contact forces between particles :math:`i` and :math:`j`, such as lubrication `[3] <https://doi.org/10.1002/aic.690400418>`_.
* :math:`\mathbf{f}_{pf,i}` is the force exerted by the surrounding fluid over particle :math:`i`;
* :math:`\mathbf{f}_{g,i}` is the gravity force over the particle :math:`i`.

The :math:`\mathbf{f}_{pf,i}` term is responsible for the momentum transport between phases, written as:

.. math::
    \mathbf{f}_{pf,i} = 



.. math::
    \nabla \cdot \mathbf{u} &= 0   \\
    \frac{\partial \mathbf{u}}{\partial t}  + \mathbf{u} \cdot \nabla \mathbf{u} &= -\frac{1}{\rho} \nabla p  + \nu \nabla^2 \mathbf{u} +\mathbf{f}


where:

* :math:`\mathbf{u}` is the velocity of the fluid. :math:`\mathbf{u}` is a vector such that :math:`\mathbf{u}=[u,v]^T` in 2D and :math:`\mathbf{u}=[u,v,w]^T` in 3D.

* :math:`p` is the pressure

* :math:`\nabla` is the `del operator <https://en.wikipedia.org/wiki/Del>`_

* :math:`\rho` is the density of the fluid.

* :math:`\nu` is the `kinematic viscosity <https://en.wikipedia.org/wiki/Viscosity>`_ of the fluid.

* :math:`\mathbf{f}` is a momentum source term.

In `Einstein notation <https://en.wikipedia.org/wiki/Einstein_notation>`_, they can be written:


.. math::
    \partial_l u_l &= 0 

    \partial_t u_k + u_l \partial_l u_k &= -\frac{1}{\rho} \partial_k p + \nu \partial_l \partial_l u_k + f_k


In Lethe, instead of solving for pressure directly, we solve for a reduced pressure :math:`p^*=\frac{p}{\rho}`. This has two advantages, first the scaling of the pressure is more adequate with respect to that of the velocity and, secondly, for single phase simulations, this means that we require a single physical property: the kinematic viscosity. Consequently, unless specified otherwise, you should assume that pressure actually refers to the reduced pressure :math:`p^*`. Although it may seem confusing at first, this is standard in incompressible fluid dynamics.

By default, Lethe does not use any units for Length (:math:`L`), Time (:math:`T`), Mass (:math:`M`) and Temperature (:math:`T`). Although this may seem confusing at first, this is relatively easy to work with. By default, most cases can be formulated using SI units (meters, seconds, kilograms and Celsius or Kelvin). However, in some cases, this may lead to a large imbalance in the scaling of different variables (e.g. the numerical value being 1000x that of the velocity). In these cases, the user must manually rescale the problem to a better balanced set of units (i.e. by using centimeter units instead of meters).


Reference
---------------
[1] Bérard, Patience, and Blais. Experimental methods in chemical engineering: Unresolved CFD‐DEM. The Canadian Journal of Chemical Engineering, v. 98, n. 2, p. 424-440, 2020. `DOI <https://doi.org/10.1002/cjce.23686>`_.
[2] Zhou, Kuang, Chu, and Yu, Discrete particle simulation of particle–fluid flow: model formulations and their applicability, Journal of Fluid Mechanics, vol. 661, pp. 482–510, 2010. `DOI <https://doi.org/10.1017/S002211201000306X>`_.
[3] Kim, Sangtae, and Karrila. Microhydrodynamics: principles and selected applications. Courier Corporation, 2013. `DOI <https://doi.org/10.1002/aic.690400418>`_.
