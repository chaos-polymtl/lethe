===========================
Unresolved CFD-DEM 
===========================

Unresolved CFD-DEM is a technique with high potential for designing and analyzing multiphase flows involving particles and fluid. Some examples of these systems are fluidized beds, stirred-tanks, and flocculation processes. In this approach, we apply Newton's second law of motion to each particle individually such that their movement is described at a micro-scale (as in DEM simulations). On the other hand, the fluid is represented at a meso-scale by a mesh of cells, to which we apply the Volume Average Navier-Stokes (VANS) equations.

The micro-meso scale approach allows for particle-fluid simulations involving large numbers of particles with reasonable computational cost and highly detailed results (in both time and space). As a counterpart, the interchanged momentum between phases needs to be modeled, i.e., it is not obtained straightforwardly from the application of a no-slip boundary condition at the interface between the fluid and the particles. The following image represents the micro-meso scale approach applied in unresolved CFD-DEM simulations, where the rectangles represent mesh of the geometry and the gray spots represent the particles.

.. image:: images/schematic_unresolve_cfd-dem.png
    :alt: Schematic represantion of micro-meso scale approach in unresolved CFD-DEM
    :align: center
    :name: geometry
    :scale: 40

In this guide, we summarize the theory behind Unresolved CFD-DEM. For further details, we refer the reader to the articles by Bérard *et al.* [#berard2020]_, and Zhou *et al.* [#zhou2010]_.

Particles
----------

Applying Newton's second law and Euler's law of angular motion on the particle :math:`i` surrounded by fluid a :math:`f`, we find:

.. math::
    m_i \frac{\mathrm{d}\mathbf{v}_i}{\mathrm{d}t} = \mathbf{F}_{\mathrm{fp},i} + \mathbf{F}_{\mathbf{g},i} + \sum_{j} \left( \mathbf{F}_{\mathrm{c},ij} + \mathbf{F}_{\mathrm{nc},ij} \right) \\
    I_i \frac{\mathrm{d}\mathbf{\omega}_i}{\mathrm{d}t} = \mathbf{M}_{\mathrm{fp},i} + \sum_{j}\left ( \mathbf{M}_{\mathrm{c},ij} + \mathbf{M}_{\mathrm{r},ij} \right ) 

where:

* :math:`m_i` is the mass of the particle :math:`i`;
* :math:`\mathbf{v}_i` is the velocity of the particle :math:`i`;
* :math:`\mathbf{F}_{\mathrm{fp},i}` is the force exerted by the surrounding fluid over particle :math:`i`. Here, the subscript :math:`fp` indicates the force exerted by the fluid on the particles;
* :math:`\mathbf{F}_{\mathbf{g},i}` is the gravitational force;
* :math:`\mathbf{F}_{\mathrm{c},ij}` are the contact forces (normal and tangential) between particles :math:`i` and :math:`j` (detailed in the DEM section of this guide);
* :math:`\mathbf{F}_{\mathrm{nc},ij}` are the non-contact forces between particles :math:`i` and :math:`j`, such as cohesive and lubrication forces [#nitsche1994]_;
* :math:`I_i` is the moment of inertia;
* :math:`\mathbf{\omega}_i` is the angular velocity;
* :math:`\mathbf{M}_{\mathrm{fp},i}` is the torque exerted by the surrounding fluid over the particle :math:`i`.
* :math:`\mathbf{M}_{\mathrm{c},ij}` is the contact torque between particles :math:`i` and :math:`j`. Since only spherical particles are currently supported the unresolved CFD-DEM implementation of Lethe, this torque is only due to tangential forces;
* :math:`\mathbf{M}_{\mathrm{r},ij}` is the rolling friction between particles :math:`i` and :math:`j`;

The momentum transfer between phases, :math:`\mathbf{F}_{\mathrm{fp},i}`, can be written as:

.. math::
    \mathbf{F}_{fp,i} = \mathbf{F}_{\nabla p,i} + \mathbf{F}_{\nabla \cdot \mathbf{\tau},i} + \mathbf{F}_{\mathrm{d},i} + \mathbf{F}_{\mathrm{Ar},i} + \mathbf{F}_{\mathrm{S},i} + \mathbf{F}_{\mathrm{M},i} + \mathbf{F}''_{i}

where:

* :math:`\mathbf{F}_{\nabla p,i}` is the force due to the pressure gradient;
* :math:`\mathbf{F}_{\nabla \cdot \tau,i}` is the force due to the shear stress;
* :math:`\mathbf{F}_{\mathrm{d},i}` is the drag force;
* :math:`\mathbf{F}_{\mathrm{Ar},i}` is the buoyancy (Archimedes) force;
* :math:`\mathbf{F}_{\mathrm{S},i}` is the Saffman lift force;
* :math:`\mathbf{F}_{\mathrm{M},i}` is the Magnus lift force;
* :math:`\mathbf{F}''_{i}` are the remaining forces, including virtual mass, Basset which are currently not implemented in Lethe.

.. note::
    Since the pressure in Lethe does not account for the hydrostatic pressure, i.e., the gravity term is not taken into account in the Navier-Stokes equations (see :doc:`../../multiphysics/fluid_dynamics/navier-stokes`), we explicitly insert :math:`\mathbf{f}_{Ar,i}` in :math:`\mathbf{f}_{pf,i}`.  

In unresolved CFD-DEM, the drag force is calculated using correlations (frequently called drag models). The drag models implemented in Lethe are described in the `unresolved CFD-DEM parameters guide <../../../parameters/unresolved-cfd-dem/cfd-dem>`_.

Volume-Averaged Navier-Stokes
------------------------------

Since we represent the fluid at a meso-scale, the quantities calculated for the cells are filtered. Additionally, variations in time and space of the volume occupied by the fluid are accounted for by the void fraction (or porosity). The resulting equations are the Volumed-Averaged Navier-Stokes (VANS) equations, which are solved to obtain the velocity and the pressure of the continuous fluid phase. In Lethe, the VANS equations are presented in two different formulations, so-called Model A (or Set II) and Model B (or Set I) `[2] <https://doi.org/10.1017/S002211201000306X>`_.

The VANS equations can be derived by averaging the continuity and momentum conservation equations using a weighting function. This approach was developed by Anderson and Jackson [#Anderson1967]_ and is briefly presented here. In this approach, any weighting function :math:`k_r(\lVert \mathbf{x} \rVert)` may be chosen, provided that it satisfies a set of conditions, including :math:`k_r(\lVert \mathbf{x} \rVert) > 0, \: \forall \mathbf{x}`, and :math:`\int_{\Omega} k_r(\lVert \mathbf{x} \rVert)\,\mathrm{d}\mathbf{x} = 1`, where :math:`\Omega` refers to the whole domain.

Using this function, local mean values of any fluid- or solid-phase property, :math:`\bar{a}(\mathbf{x},t)` and :math:`\bar{b}(\mathbf{x},t)`, respectively, can be defined at a scale larger than the particle diameter yet smaller than the macroscopic length scale of the problem, as per the following:

.. math:: 
    :label: eq:average_definition

        \begin{align}
            \varepsilon_f(\mathbf{x},t) \bar{a}(\mathbf{x},t) &= \int_{\Omega_f} a(\mathbf{y},t) k_r \left (\lVert \mathbf{x} - \mathbf{y} \rVert \right ) \mathrm{d}\mathbf{y} \\
            \left[1-\varepsilon_f(\mathbf{x},t)\right] \bar{b}(\mathbf{x},t) &= \int_{\Omega_s} b(\mathbf{y},t) k_r \left (\lVert \mathbf{x} - \mathbf{y} \rVert \right ) \mathrm{d}\mathbf{y} \approx \sum_{i} V_{\mathrm{p},i} \bar{b}_i(t)k_r(\lVert \mathbf{x}-\mathbf{x_i} \rVert)
        \end{align}

where :math:`\varepsilon_f(\mathbf{x},t)` is the void fraction at position :math:`\mathbf{x}` and time :math:`t`, the subscripts :math:`f` and :math:`s` denote the fluid and solid, respectively, and :math:`a` and :math:`b` are the point values of the corresponding properties. The index :math:`i` in the summation denotes again each particle in the system, with volume :math:`V_{\mathrm{p},i}` and position :math:`\mathbf{x}_i`. :math:`\bar{b}_i(t)` represents the average of the point property :math:`b` over the volume of particle :math:`i`. For the solid-phase property :math:`b(\mathbf{x},t)`, the approximation in the second equation stems from the assumption that the weighting function varies little over the volume of a single particle. This condition must be met when selecting the weighting kernel radius. Filtering the mass and momentum conservation equations using the weighting function leads to the VANS equations presented in the following.

Considering an incompressible flow, the continuity equation for both models A and B is:

.. math::
    \frac{\partial \varepsilon_f}{\partial t} + \nabla \cdot \left ( \varepsilon_f \mathbf{u} \right ) = 0

where:

* :math:`\mathbf{u}` is the filtered fluid velocity vector;
* :math:`\varepsilon_f` is the void fraction.

Models A and B differ from each other in the way the momentum equation is calculated:

Model A:

.. math:: 
    \rho_f \left ( \frac{\partial \left ( \varepsilon_f \mathbf{u} \right )}{\partial t} + \nabla \cdot \left ( \varepsilon_f \mathbf{u} \otimes \mathbf{u} \right ) \right ) = -\varepsilon_f \nabla p + \varepsilon_f \nabla \cdot \tau + \mathbf{f}_{\mathrm{pf},\mathrm{A}}

Model B:

.. math:: 
    \rho_f \left ( \frac{\partial \left ( \varepsilon_f \mathbf{u} \right )}{\partial t} + \nabla \cdot \left ( \varepsilon_f \mathbf{u} \otimes \mathbf{u} \right ) \right ) = -\nabla p + \nabla \cdot \tau +  \mathbf{f}_{\mathrm{pf},\mathrm{B}}

where:

* :math:`\rho_f` is the density of the fluid;
* :math:`p` is the pressure;
* :math:`\tau` is the deviatoric stress tensor;
* :math:`\mathbf{f}_\mathrm{A}^{\mathrm{pf}}` and :math:`\mathbf{f}_\mathrm{B}^{\mathrm{pf}}` denote the interphase momentum exchange terms resulting from fluid–particle interactions in Models A and B, respectively. The term :math:`\mathbf{f}_{\mathrm{pf}}`, which has units of a volumetric force, is obtained by regularizing the particle-fluid interaction forces using a kernel function :math:`k_r` :


.. math::
     \mathbf{f}_{\mathrm{pf},\mathrm{A}} &= \sum_{i} k_r \left (\lVert \mathbf{x} - \mathbf{x}_i \rVert \right ) \left(- \left (\mathbf{F}_{\mathrm{fp},i} - \mathbf{F}_{\nabla p,i} - \mathbf{F}_{\nabla \cdot \mathbf{\tau},i} \right ) \right) \\
     \mathbf{f}_{\mathrm{pf},\mathrm{B}} &= \sum_{i} k_r \left (\lVert \mathbf{x} - \mathbf{x}_i \rVert \right ) \left(-\mathbf{F}_{{\mathrm{fp},i}} \right)            

Different choices of the weighting function :math:`k_r(\lVert \mathbf{x} \rVert)` are possible, and Lethe provides two options. When the Particle-in-Cell (PIC) method is used for particle–fluid coupling, particle properties are regularized over the mesh cells by averaging them within each cell. This procedure is equivalent to using a cell-based top-hat weighting function over the cell containing the point :math:`\mathbf{x}`. However, this approach is not ideal, as it leads to force definitions that depend on the mesh resolution.

To mitigate this issue, Lethe also supports filtering of solid–fluid interaction forces using a spherical top-hat filter centered at the quadrature points, with a user-defined radius that should be larger than the particle diameter. This approach, referred to as the Quadrature-Centered Method (QCM), yields a filtering operation that is independent of the mesh and depends solely on the chosen filter radius.

Lethe is capable of simulating unresolved CFD-DEM cases with both Models A and B (see the :doc:`../../../parameters/unresolved-cfd-dem/cfd-dem` page of this guide).

Void Fraction
--------------
Determining the void fraction is an important step in unresolved CFD-DEM, as can be noted by the VANS equations and the drag models [#rong2013]_. There exist several methods for the calculation of the void fraction in a CFD-DEM simulation. Some are approximations while others are analytical approaches. If one applies the averaging definition in equations :eq:`eq:average_definition` while taking :math:`b(\mathbf{x},t)=1`, one obtains:

.. math::
    1-\varepsilon_f(\mathbf{x},t) \approx \sum_{i} V_{\mathrm{p},i} k_r(\lVert \mathbf{x}-\mathbf{x_i} \rVert)

This equation can be used in practice to calculate the void fraction at a given point :math:`\mathbf{x}`:

.. math:: 
    :label: eq:eps_definition

    \varepsilon_f(\mathbf{x},t) \approx 1 - \sum_{i} V_{\mathrm{p},i} k_r(\lVert \mathbf{x}-\mathbf{x_i} \rVert)

In the finite element method, the void fraction is initially calculated at quadrature points. However, the solution of the VANS equation requires not only the void fraction, but its gradient. Therefore, it is necessary to be able to represent the void fraction using a finite element interpolation (e.g. Lagrange Q elements). The best polynomial representation of the void fraction can be obtained by solving the following minimization problem (which is essentially an  :math:`\mathcal{L}^2` projection  [#larson2013]_):

.. math:: 
    \min_{\varepsilon_{f,j} \in \mathbb{R}} \frac{1}{2} \sum_i \left (\sum_j \varepsilon_{f,j} \varphi_j - \epsilon_{f} \right )^2 \varphi_i

where :math:`\epsilon_{f}` is the void fraction, :math:`\varphi_j` is the finite element shape function of the void fraction, and :math:`\varepsilon_{f,j}` the nodal values of the projected void fraction. Then, we assemble and solve the following:

.. math::
    \int_{\Omega} \varepsilon_{f,j} \varphi_i  \varphi_j d \Omega = \int_{\Omega}  \epsilon_{f} \varphi_i d \Omega


Lethe also has the option of smoothing the void fraction profile, which helps to mitigate sharp discontinuities. This is specifically advantageous when using void fraction schemes that are discontinuous in space and time such as the Particle Centroid Method (PCM) and the Satellite Point Method (SPM). To do so, we add to the left hand side of the previous equation a term similar to a Poisson equation. This acts as an additional parabolic (Gaussian) filter. The resulting equation is:

.. math::
    \int_{\Omega} \varepsilon_{f,j} \varphi_i  \varphi_j d \Omega +  \iint_\Omega L^2 \varepsilon_{f,j} \nabla \varphi_i \nabla \varphi_j d\Omega = \int_{\Omega} \epsilon_{f} \varphi_i d \Omega

where :math:`L` is the smoothing length, used as parameter in Lethe unresolved CFD-DEM simulations. In Lethe, three void fraction schemes are currently supported to calculate :math:`\epsilon_{f}`. They are the particle centroid method, the satellite point method, and the quadrature centered method.

The Particle Centroid Method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The PCM [#peng2014]_ is a simple and commonly used method. It consists of tracking the position of the centroid of each particle and applying the total volume of the particle to the calculation of the void fraction of the cell. This means that in either of the following situations the void fraction of the colored cell is the same:

.. image:: images/void_frac1.png
   :width: 49% 
.. image:: images/void_frac2.png
   :width: 49%

This results in the PCM being discontinuous in space and time. Consequently, the PCM method significantly relies on the smoothing of the void fraction to lead to a sufficiently smooth void fraction field. We refer the reader to [#geitani2023]_  for an extensive discussion on this topic. 

The void fraction in a cell using PCM can be written as:

.. math:: 
    \epsilon_f = 1 - \frac{\sum_{i,\Omega_e}^{n_\mathrm{p}} V_{\mathrm{p},i}}{V_{\Omega_e}}

where :math:`n_\mathrm{p}` is the number of particles whose centroids lie inside the cell :math:`\Omega_e` of volume :math:`V_{\Omega_e}`. Comparing this equation with equation :eq:`eq:eps_definition`, we can see that this is equivalent to using the following weighting function:

.. math::
    k_r \left (\lVert \mathbf{x} - \mathbf{x}_i \rVert \right ) = \frac{1_{\{\mathbf{x}_i \in \Omega_e\}}}{V_{\Omega_e}}

where :math:`\mathbf{1}_{\{\mathbf{x}_i \in \Omega_e\}}` is the indicator function that is equal to 1 if the centroid of particle :math:`i` lies inside cell :math:`\Omega_e` and 0 otherwise, and :math:`V_{\Omega_e}` is the volume of cell :math:`\Omega_e`. It can be verified that using this weighting function for calculating the void fraction and the interaction forces ensures the conservation of mass of both phases, as well as Newton's third law:

.. math::

    \begin{align}
        \int_{\Omega} \mathbf{f}_{\mathrm{pf},\mathrm{A}} d\mathbf{x} &= \sum_{i} -\left(\mathbf{F}_{\mathrm{fp},i}- \mathbf{F}_{\nabla p,i} - \mathbf{F}_{\nabla \cdot \mathbf{\tau},i}\right)\\
        \int_{\Omega} \mathbf{f}_{\mathrm{pf},\mathrm{B}} d\mathbf{x} &= \sum_{i} \left(-\mathbf{F}_{\mathrm{fp},i}\right)
    \end{align}

The Satellite Point Method
~~~~~~~~~~~~~~~~~~~~~~~~~~
This method divides each particle into pseudo-particles where the sum of the volume of all pseudo-particles in a single particle is equal to the volume of the particle. Then, each pseudo-particle is treated similarly to the PCM, that is, the centroid of each pseudo-particle is tracked, and the entire volume of the pseudo-particle is considered in a given cell if its centroid lies within. 

.. image:: images/spm.png
   :width: 49% 
   :align: center
   
The void fraction in a cell using SPM can be written as: 

.. math:: 
       \epsilon_f = 1 - \frac{\sum_{i}^{n_\mathrm{p}}\sum_{j,\Omega_e}^{n_{\mathrm{sp}}} V_{\mathrm{sp},j}}{V_{\Omega_e}}

where :math:`n_{\mathrm{sp}}` is the number of pseudo-particles j belonging to particle i and with a centroid inside the cell :math:`\Omega_e` with volume :math:`V_{\Omega_e}`, and :math:`V_{\mathrm{sp}}` is the volume of the satellite point.  The corresponding weighting function is thus:

.. math::
    k_r \left (\lVert \mathbf{x} - \mathbf{x}_i \rVert \right ) = \frac{\sum_{i}^{n_\mathrm{p}}(1/V_{\mathrm{p},i})\sum_{j,\Omega_e}^{n_{\mathrm{sp}}} V_{\mathrm{sp},j}}{V_{\Omega_e}}

This function is similar to the weighting function used in the Particle-in-Cell (PCM) method, but it varies gradually (albeit discontinuously) at the cell boundaries. Hence, it suffers from the same limitations as the PCM method. However, it also conserves mass. Lethe only supports the SPM for the calculation of the void fraction. For the calculation of the force with SPM, Lethe uses PCM.

The Quadrature Centered Method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The Quadrature Centered Method (QCM) [#geitani2023]_  is an analytical method that decouples the averaging volume from the mesh cells. It constructs an averaging sphere centered at each quadrature point in a given cell, and it calculates the void fraction directly in the averaging volume at the quadrature point. Since the sphere-sphere (particle-averaging sphere) intersection is analytically easier to calculate than sphere-polyhedron (particle-mesh cell), this method is less expensive than other analytical methods as the intersection does not involve the calculation of trigonometric functions at each CFD time step. The advantage of this method is that the void fraction varies within a cell. Additionally, particles in neighboring cells can affect the void fraction of the current cell. This allows the method to be continuous in both space and time. This is advantageous, especially in solid-liquid systems where the term :math:`\rho_f \frac{\partial \epsilon_f}{\partial t}` of the continuity equation is very stiff and unstable, when there exist even small discontinuities in the void fraction, and where it explodes when :math:`\Delta t_{CFD} \to 0`.

An averaging volume sphere is constructed around each quadrature point. All particles lying in the sphere will contribute to the void fraction value of this quadrature point. Therefore, a cell will be affected by the particles lying in it and in its neighboring cells.

.. image:: images/qcm.png
   :width: 49% 
   :align: center

The void fraction at the quadrature point using QCM is:

.. math:: 
      \epsilon_f (\mathbf{x}_q,t) = 1-\sum_{i}{\frac{V_{\mathrm{p},i}}{\sum_{q}{V_{\mathrm{p}\cap S_q}}}V_{\mathrm{p}\cap S_Q}\frac{1}{w_q|J|}}
    
where :math:`S_q` is the sphere centered at :math:`\mathbf{x}_q`, :math:`V_{p\cap S_q}` is the intersection volume of the particle with the sphere, :math:`w_q|J|` is the weight at the quadrature point located at :math:`\mathbf{x}_q`, multiplied by the volume of the cell containing the quadrature point. The division by the product of the latter term and :math:`\sum_{q}{V_{\mathrm{p}\cap S_q}}` is done to ensure mass conservation over the whole domain. The corresponding weighting function is:

.. math:: 
      k_r \left (\lVert \mathbf{x}_q - \mathbf{x}_i \rVert \right )  = \frac{1}{\sum_{q}{V_{p\cap S_q}}}V_{p\cap S_q}\frac{1}{w_q|J|}

This weighting function is also used to compute the fluid–particle momentum exchange term. Integrating this term over the computational domain demonstrates that the formulation also satisfies Newton’s third law (the corresponding expression is shown for Model B without loss of generality):

.. math:: 
      \int_{\Omega} \mathbf{f}_{\mathrm{pf}}  d\mathbf{x} = \sum_{q} \sum_{i}{\left(-\mathbf{F}_{{\mathrm{fp},i}} \right)\frac{1}{\sum_{q}{V_{p\cap S_q}}}V_{p\cap S_q}\frac{1}{w_q|J|}} \times w_q|J| = \sum_{i} \left(-\mathbf{F}_{{\mathrm{fp},i}} \right)

The averaging spheres must be sufficiently large to ensure that all particles in the domain contribute to the averaged void fraction, which can be achieved by selecting an adequate number of quadrature points together with an appropriate averaging radius, :math:`R_s`. However, in Lethe, the radius of the averaging spheres is currently constrained by the mesh resolution. Since each cell only has access to its direct neighboring cells, the averaging volume cannot extend beyond this local neighborhood. Consequently, the radius of the averaging spheres must satisfy:

.. math:: 
    R_s \leq h_{\Omega}

where :math:`h_{\Omega}` denotes the characteristic cell size. This limitation may be relaxed in future versions of Lethe.

References
-----------
.. [#berard2020] \A. Bérard, G. S. Patience, and B. Blais, “Experimental methods in chemical engineering: Unresolved CFD-DEM,” *Can. J. Chem. Eng.*, vol. 98, no. 2, pp. 424–440, 2020, doi: `10.1002/cjce.23686 <https://doi.org/10.1002/cjce.23686>`_\.

.. [#zhou2010] \Z. Y. Zhou, S. B. Kuang, K. W. Chu, and A. B. Yu, “Discrete particle simulation of particle–fluid flow: model formulations and their applicability,” *J. Fluid Mech.*, vol. 661, pp. 482–510, Oct. 2010, doi: `10.1017/S002211201000306X <https://doi.org/10.1017/S002211201000306X>`_\.

.. [#Anderson1967] \T. B. Anderson and R. Jackson, “A Fluid Mechanical Description of Fluidized Beds,” *Ind. Eng. Chem. Fundam.*, vol. 6, no. 4 pp. 527–539, Nov. 1967, doi: `10.1021/i160024a007 <https://doi.org/10.1021/i160024a007>`_\.

.. [#nitsche1994] \L. C. Nitsche, “Microhydrodynamics: Principles and selected applications. By Sangtae Kim and Seppo J. Karrila, Butterworth-Heinemann, Boston, 1991” *AIChE J.*, vol. 40, no. 4, pp. 739–743, 1994, doi: `10.1002/aic.690400418 <https://doi.org/10.1002/aic.690400418>`_\.

.. [#rong2013] \L. W. Rong, K. J. Dong, and A. B. Yu, “Lattice-Boltzmann simulation of fluid flow through packed beds of uniform spheres: Effect of porosity,” *Chem. Eng. Sci.*, vol. 99, pp. 44–58, Aug. 2013, doi: `10.1016/j.ces.2013.05.036 <http://dx.doi.org/10.1016/j.ces.2013.05.036>`_\.

.. [#peng2014] \Z. Peng, E. Doroodchi, C. Luo, and B. Moghtaderi, “Influence of void fraction calculation on fidelity of CFD-DEM simulation of gas-solid bubbling fluidized beds,” *AIChE J.*, vol. 60, no. 6, pp. 2000–2018, 2014, doi: `10.1002/aic.14421 <https://doi.org/10.1002/aic.14421>`_\.

.. [#larson2013] \M. G. Larson and F. Bengzon, *The Finite Element Method: Theory, Implementation, and Applications*. Springer Science & Business Media, 2013. https://link.springer.com/book/10.1007/978-3-642-33287-6\.

.. [#geitani2023] \T. El Geitani and B. Blais, “Quadrature-Centered Averaging Scheme for Accurate and Continuous Void Fraction Calculation in Computational Fluid Dynamics–Discrete Element Method Simulations,” *Ind. Eng. Chem. Res.*, vol. 62, no. 12, pp. 5394–5407, Mar. 2023, doi: `10.1021/acs.iecr.3c00172 <https://doi.org/10.1021/acs.iecr.3c00172>`_\.
