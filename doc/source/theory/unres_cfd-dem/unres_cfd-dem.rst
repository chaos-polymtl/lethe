Unresolved CFD-DEM coupling
############################

Unresolved CFD-DEM is a technique with high potential for designing and analyzing multiphase flows involving particles and fluid. Some examples of these systems are fluidized beds, stirred tanks, and flocculation processes. In this approach, we apply Newton's second law to each particle individually such that their movement is described at a micro-scale (as in DEM simulations). On the other hand, the Volume Average Navier-Stokes (VANS) equations describe the fluid at a meso-scale. The micro-meso scale approach allows for particle-fluid simulations involving large numbers of particles with reasonable computational cost and highly detailed results (in both time and space). As a counterpart, the interchanged momentum between phases needs to be modeled, i.e., it is not resolved. The following image represents the micro-meso scale approach applied in unresolved CFD-DEM simulations, where the rectangles represent subdomains of the geometry and the gray spots represent the particles.

.. image:: images/schematic_unresolve_cfd-dem.jpg
    :alt: Schematic represantion of micro-meso scale approach in unresolved CFD-DEM
    :align: center
    :name: geometry

In this guide, we summarize the theory behind Unresolved CFD-DEM. For further details, we refer the reader to the articles by Bérard, Patience & Blais (2019) `[1] <https://doi.org/10.1002/cjce.23686>`_, and Zhou et al. (2010) `[2] <https://doi.org/10.1017/S002211201000306X>`_.

Applying Newton's second law on the particle :math:`i` surrounded by fluid a :math:`f`, we find:

.. math::
    m_i \frac{d \mathbf{v}_i}{dt} = \sum_{j}\mathbf{f}_{c,ij} + \sum_{j}\mathbf{f}_{nc,ij} + \mathbf{f}_{pf,i} + \mathbf{f}_{g,i} \\
    I_i \frac{d\mathbf{\omega}_i}{dt} = \sum_{j}\left ( \mathbf{M}_{c,ij} + \mathbf{M}_{r,ij} \right ) + \sum_{w}\left ( \mathbf{M}_{c,iw} + \mathbf{M}_{r,iw} \right )

where:

* :math:`m_i` is the mass of the particle;
* :math:`v_i` is the velocity vector;
* :math:`\mathbf{f}_{c,ij}` are the contact forces between particles :math:`i` and :math:`j` (detailed in the DEM section of this guide);
* :math:`\mathbf{f}_{nc,ij}` are the non-contact forces between particles :math:`i` and :math:`j`, such as lubrication `[3] <https://doi.org/10.1002/aic.690400418>`_.
* :math:`\mathbf{f}_{pf,i}` is the force exerted by the surrounding fluid over particle :math:`i`;
* :math:`\mathbf{f}_{g,i}` is the gravity force;
* :math:`I_i` is the momentum of inertia;
* :math:`\mathbf{\omega}_i` is the angular velocity;
* :math:`\mathbf{M}_{c,ij}` is the torque between particles :math:`i` and :math:`j`;
* :math:`\mathbf{M}_{r,ij}` is the rolling friction between particles :math:`i` and :math:`j`;
* :math:`\mathbf{M}_{c,iw}` is the torque between particle :math:`i` and walls :math:`w`;
* :math:`\mathbf{M}_{c,iw}` is the rolling friction between particle :math:`i` and walls :math:`w`;

Apart from :math:`\mathbf{f}_{pf,i}`, all the other terms of the previous equations are detailed in the DEM section of this theory guide (**under construction**). The momentum transport between phases :math:`\mathbf{f}_{pf,i}` can be written as:

.. math::
    \mathbf{f}_{pf,i} = \mathbf{f}_{\nabla p,i} + \mathbf{f}_{\nabla \cdot \mathbf{\tau},i} + \mathbf{f}_{d,i} + \mathbf{f}_{Ar,i} + \mathbf{f}_{g,i} + \mathbf{f}''_{i}

where:

* :math:`\mathbf{f}_{\nabla p,i}` is the force due to the pressure gradient;
* :math:`\mathbf{f}_{\nabla \cdot \tau,i}` is the force due to the shear stress;
* :math:`\mathbf{f}_{d,i}` is the drag force;
* :math:`\mathbf{f}_{Ar,i}` is the buoyant (Archimedes) force;
* :math:`\mathbf{f}_{g,i}` is the force due to gravity;
* :math:`\mathbf{f}''_{i}` are the remaining forces, including virtual mass, Basset, Lift, and Magnus (currently not implemented in Lethe).

Since pressure in Lethe does not account for the hydrostatic pressure, we explicitly insert :math:`\mathbf{f}_{Ar,i}` in :math:`\mathbf{f}_{pf,i}`.

In unresolved CFD-DEM drag is calculated using correlations (frequently called drag models). The drag models implemented in Lethe are described in the `unresolved CFD-DEM parameters guide <https://lethe-cfd.github.io/lethe/parameters/unresolved_cfd-dem/cfd_dem.html>`_.

Since we represent the fluid at a meso scale, the quantities calculated for the subdomains are averages among its volume. Additionally, as the volume of fluid is a fraction of the subdomain, the porosity (or void fraction) is taken into account. To do this, we apply the Volume Average Navier-Stokes (VANS) equations to represent the fluid phase. Mainly, the VANS equations are presented in two different formulations, so called Model A (or Set II) and Model B (or Set I).

For both models, considering incompressible fluid, the continuity is:

.. math::
    \frac{\partial \varepsilon_f}{\partial t} + \nabla \cdot \left ( \varepsilon_f \mathbf{u} \right ) = 0

where:

* :math:`\mathbf{u}` is the the fluid velocity vector;
* :math:`\varepsilon_f` is the void fraction.

Models A and B differ from each other in the way the momentum equation is calculated. In Model A, we consider that pressure is in both phases, while for Model B the pressure is only in the fluid:

Model A:

.. math:: 
    \rho_f \left ( \frac{\partial \left ( \varepsilon_f \mathbf{u} \right )}{\partial t} + \nabla \cdot \left ( \varepsilon_f \mathbf{u} \otimes \mathbf{u} \right ) \right ) = -\varepsilon \nabla p + \varepsilon \nabla \cdot \tau + \mathbf{F}_{fp}^A

Model B:

.. math:: 
    \rho_f \left ( \frac{\partial \left ( \varepsilon_f \mathbf{u} \right )}{\partial t} + \nabla \cdot \left ( \varepsilon_f \mathbf{u} \otimes \mathbf{u} \right ) \right ) = -\nabla p + \nabla \cdot \tau + \mathbf{F}_{fp}^B

where:

* :math:`\rho_f` is the density of the fluid;
* :math:`p` is the pressure;
* :math:`\tau` is the shear stress;
* :math:`\mathbf{F}_{fp}^A` and :math:`\mathbf{F}_{fp}^B` are the source terms representing the forces applied back in the fluid due to the interaction with particles for Models A and B, respectively.

For Model A, since the pressure term corresponds to a 'fluid fraction of the pressure', we can write the interaction term as:

.. math:: 
    \mathbf{F}_{fp}^A = \frac{1}{\Delta V}\sum_{i}^{n_p}\left ( \mathbf{f}_{pf, i} - \mathbf{f}_{\nabla p, i} - \mathbf{f}_{\nabla \cdot \tau, i} \right )

while for Model B, since the pressure is totally in the fluid, we write:

.. math:: 
    \mathbf{F}_{fp}^B = \frac{1}{\Delta V}\sum_{i}^{n_p}\left ( \mathbf{f}_{pf, i} \right )

where :math:`n_p` is the number of particles inside the subdomain with volume :math:`\Delta V`.



Reference
---------------
`[1] <https://doi.org/10.1002/cjce.23686>`_ Bérard, Patience, and Blais. Experimental methods in chemical engineering: Unresolved CFD‐DEM. The Canadian Journal of Chemical Engineering, v. 98, n. 2, p. 424-440, 2020.

`[2] <https://doi.org/10.1017/S002211201000306X>`_ Zhou, Kuang, Chu, and Yu, Discrete particle simulation of particle–fluid flow: model formulations and their applicability, Journal of Fluid Mechanics, vol. 661, pp. 482–510, 2010.

`[3] <https://doi.org/10.1002/aic.690400418>`_ Kim, Sangtae, and Karrila. Microhydrodynamics: principles and selected applications. Courier Corporation, 2013.
