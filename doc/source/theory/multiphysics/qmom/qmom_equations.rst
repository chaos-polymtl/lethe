==============
QMOM Equations
==============

Lethe can solve population balance equations with the quadrature method of moments (QMOM) [#mcgraw1997]_ [#marchisio2003]_. The evolution of the moments in phase space, as well as the computation of the moments of particle size distributions, are done by the `QMOMKit <https://github.com/polycfd/QMOMKit>`_ library. The parameters of the QMOM physics are described in :doc:`../../../parameters/cfd/qmom`.

Transported Quantities
----------------------

The QMOM physics transports:

* :math:`n_s` user-defined scalar fields;
* :math:`K` moments :math:`M_0, \dots, M_{K-1}` of the particle size distribution. The quadrature has :math:`N = \lceil K/2 \rceil` nodes, since :math:`N` nodes are defined by :math:`K=2N` moments.

The moment of order :math:`k` of a particle size distribution is:

.. math::

  M_k = n \, \mathrm{E}[L^k]

where :math:`n` is the number density and :math:`L` the size of the particles. In particular, :math:`M_0` is the number density, and :math:`M_3` is proportional to the volume of the particles per unit volume.

All the scalars and moments are discretized with the same scalar discontinuous Galerkin finite element, on a single set of degrees of freedom, and each of them is stored in its own vector.

Operator Splitting
------------------

Each time step :math:`\Delta t` of the fluid dynamics is split into the transport of the scalars and moments in physical space, :math:`T`, and their evolution in phase space, :math:`P`. Two splitting schemes are defined:

* Lie splitting: :math:`T(\Delta t)` followed by :math:`P(\Delta t)`. This scheme is first-order accurate in time. The transport :math:`T` solves the moments first, then the scalars, so that the phase-space evolution uses the scalars at the end of their transport.
* Strang splitting: :math:`P(\Delta t / 2)`, then :math:`T(\Delta t)`, then :math:`P(\Delta t / 2)`. This scheme is second-order accurate in time when both parts are at least second-order accurate.

.. warning::

  Only the Lie splitting is currently implemented.

Transport in Physical Space
---------------------------

Every component :math:`\zeta_c`, scalar or moment, is transported by the velocity :math:`\mathbf{u}` of the fluid at the end of the time step, which is solved before the QMOM physics, with a source term :math:`S_c`:

.. math::

  \frac{\partial \zeta_c}{\partial t} + \nabla \cdot (\mathbf{u} \zeta_c) = S_c

The current implementation assumes a purely advective transport: the scalars and moments do not diffuse.

The components are ordered as the scalars, followed by the moments: :math:`\zeta_c` is the scalar :math:`c` for :math:`c < n_s`, and the moment :math:`M_{c - n_s}` otherwise.

Linear Systems
^^^^^^^^^^^^^^

Since all the scalars and moments share the same element and the same velocity, they share the same linear system matrix: it is assembled and preconditioned once per time step, then solved for the right-hand side of every moment, and then of every scalar. The linear systems are solved with the GMRES method and an ILU or AMG preconditioner.

The equations are linear, so no non-linear solver is used. Since the moments have very different magnitudes, the tolerance of each scalar and moment is a relative residual times its initial residual, computed with the solution at the beginning of the time step as initial guess, but at least a minimum residual. The minimum residual should therefore be much smaller than the initial residual of the smallest moment.

Evolution in Phase Space
------------------------

The evolution in phase space is done by QMOMKit, independently at every degree of freedom: since the element is nodal and discontinuous, the moments of a degree of freedom are the moments at its support point. They are inverted into a quadrature of :math:`N` nodes, of weights :math:`w_i` and abscissas (sizes) :math:`\xi_i`, with the Wheeler algorithm, the quadrature is evolved over the time step, and the moments are rebuilt from it:

.. math::

  M_k = \sum_{i=1}^{N} w_i \xi_i^k

Three mechanisms evolve the quadrature.

Growth
^^^^^^

The growth moves the abscissas with the growth rate :math:`G`, the velocity of the particle size in phase space, while the weights are unchanged:

.. math::

  \frac{\mathrm{d} \xi_i}{\mathrm{d} t} = G(\xi_i, t)

The trajectories are integrated with an adaptive Runge-Kutta method. Two particular growth rates have simple solutions:

* a constant growth rate :math:`G = a` increases all the sizes by the same amount :math:`a t`;
* a linear growth rate :math:`G = a \xi` multiplies all the sizes by :math:`\exp(a t)`, so that a log-normal distribution stays log-normal, with the mean of the logarithm of the size increased by :math:`a t`, and :math:`M_k(t) = M_k(0) \exp(k a t)`.

Nucleation
^^^^^^^^^^

The nucleation creates particles at the rate :math:`J`, the number of particles created per unit volume and time, whose sizes follow a log-normal distribution of mean :math:`m` and standard deviation :math:`\sigma` [#limpert2001]_ [#marchisio2013]_:

.. math::

  \frac{\mathrm{d} M_k}{\mathrm{d} t} = J \exp\left(k \mu + \frac{k^2 s^2}{2}\right), \quad s^2 = \ln\left(1 + \frac{\sigma^2}{m^2}\right), \quad \mu = \ln(m) - \frac{s^2}{2}

The standard deviation :math:`\sigma` must be strictly positive: nuclei of a single size would have moments that can only be inverted into a single quadrature node.

Coalescence
^^^^^^^^^^^

The coalescence merges two particles of sizes :math:`\xi_i` and :math:`\xi_j` at the rate given by the kernel :math:`\beta`, into a particle of the same volume, of size :math:`(\xi_i^3 + \xi_j^3)^{1/3}`. The moment :math:`M_3` is therefore conserved [#marchisio2003]_ [#marchisio2003b]_:

.. math::

  \frac{\mathrm{d} M_k}{\mathrm{d} t} = \frac{1}{2} \sum_{i=1}^{N} \sum_{j=1}^{N} w_i w_j \beta(\xi_i, \xi_j, t) \left[ \left(\xi_i^3 + \xi_j^3\right)^{k/3} - \xi_i^k - \xi_j^k \right]

With a constant kernel :math:`\beta = a`, the number density decreases as :math:`\mathrm{d} M_0 / \mathrm{d} t = -a M_0^2 / 2` [#friedlander2000]_.

Time Integration in Phase Space
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The changes of the moments due to the nucleation and the coalescence are integrated with a third-order Runge-Kutta method, and combined with the growth within the phase-space step by a Strang splitting, the default of QMOMKit. This inner splitting is distinct from the splitting of the transport and the phase-space evolution. QMOMKit also provides two first-order alternatives (growth then sources, or sources then growth), which are not exposed as parameters.

The rates of the three mechanisms can depend on the scalars of the degree of freedom, which keep their values at the end of the transport during the phase-space step: the models of the generic process do not use them, while those of the lipid nanoparticle production process compute the supersaturation of the lipid from them (see :ref:`qmom-lipid-nanoparticle-production`). Without any mechanism enabled, the quadrature does not evolve and the moments are unchanged, up to the round-off of the inversion.

Realizability
^^^^^^^^^^^^^

Degrees of freedom without particles (:math:`M_0` below :math:`10^{-30}`) have no quadrature: only the nucleation changes their moments, with its rate integrated by the midpoint rule.

The moments must nevertheless be realizable, i.e. be the moments of a distribution of positive sizes: the simulation stops with an error at the first time step where they cannot be inverted, for example because of the oscillations of the transport near steep gradients. Optionally, QMOMKit may instead reduce the number of quadrature nodes when the moments are close to the boundary of the realizable set.

.. _qmom-lipid-nanoparticle-production:

Lipid Nanoparticle Production
-----------------------------

With the lipid nanoparticle production process, the nucleation, growth and coalescence are given by the models of Shin et al. [#shin2025]_, whose parameters are in the ``lipid nanoparticle production`` subsection of :doc:`../../../parameters/cfd/qmom`: the generic growth, nucleation and coalescence models are ignored, but each mechanism is still enabled by its ``enable`` parameter. The two scalars are the ethanol concentration :math:`c_e` and the free lipid concentration :math:`C`, both in :math:`\text{kg} \cdot \text{m}^{-3}`.

Solubility and Supersaturation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At every degree of freedom, the organic (ethanol) volume fraction is :math:`\phi_o = c_e / \rho_o`, clipped to :math:`[0, 1]`, and the aqueous volume fraction is :math:`\phi_a = 1 - \phi_o`. The mole fraction solubility :math:`x` of the lipid in the mixture, its solubility :math:`C_\infty` and the supersaturation :math:`S` are:

.. math::

  \ln x = \phi_a \ln x_a + \phi_o \ln x_o, \quad C_\infty = \frac{M_l x}{1 - x} \left(\frac{\phi_a \rho_a}{M_a} + \frac{\phi_o \rho_o}{M_o}\right), \quad S = \frac{\max(C, 0)}{C_\infty}

where :math:`M_l` is the molar mass of the lipid and :math:`x_a`, :math:`\rho_a`, :math:`M_a`, :math:`x_o`, :math:`\rho_o` and :math:`M_o` are the mole fraction solubilities, densities and molar masses of the aqueous and organic phases. The capillary length is :math:`\alpha = 4 \gamma v_m / (k_B T)`, where :math:`\gamma` is the surface tension, :math:`v_m` the molecular volume of the lipid, :math:`k_B` the Boltzmann constant and :math:`T` the operating temperature.

Growth
^^^^^^

The growth includes Ostwald ripening: particles smaller than their equilibrium size dissolve.

.. math::

  G = k_g \xi C_\infty \left(S - e^{\alpha / \xi}\right)

with :math:`k_g` the growth factor.

Nucleation
^^^^^^^^^^

The nucleation only happens in a supersaturated mixture (:math:`S > 1`), at the rate:

.. math::

  J = A \exp\left(-\frac{16 \pi \gamma^3 v_m^2}{3 (k_B T)^3 (\ln S)^2}\right)

with :math:`A` the nucleation factor. The sizes of the nuclei follow a log-normal distribution of mean the critical size :math:`L_c = \alpha / \ln S` and of standard deviation the nuclei standard deviation.

Coalescence
^^^^^^^^^^^

The coalescence kernel is the Brownian kernel times the attachment efficiency of the Fuchs model, limited by the PEG factor :math:`k_{PEG}`:

.. math::

  \beta = \frac{2 k_B T (\xi_1 + \xi_2)^2}{3 \mu \xi_1 \xi_2} \min\left(\frac{k_{PEG}}{W}, k_{PEG}\right), \quad W = r_s \int_0^\infty \frac{\exp\left(\Phi(a) / (k_B T)\right)}{(a + r_s)^2} \mathrm{d} a

where :math:`\mu` is the viscosity of the mixture, :math:`r_s` the sum of the radii of the particles, :math:`a` the distance between their surfaces and :math:`\Phi` the sum of their van der Waals attraction (Hamaker constant) and electrostatic repulsion (surface potential). The Debye length uses the ion concentration of the mixture and the relative permittivity of the local mixture, :math:`\phi_a \varepsilon_a + \phi_o \varepsilon_o`. The integral is computed with the Simpson rule on 400 logarithmic intervals of :math:`a` between :math:`10^{-13}` and :math:`10^{-5}\ \text{m}`. This model is much more expensive than the constant coalescence model.

Consumption of the Free Lipids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The particles consume the free lipids: at the end of the phase-space step, the free lipid concentration of every degree of freedom loses the mass of lipid gained by its particles:

.. math::

  C \leftarrow C - \rho_l k_v \left(M_3 - M_3^0\right)

where :math:`\rho_l` is the density of the lipid in the particles, :math:`k_v` the volume shape factor and :math:`M_3^0` the moment :math:`M_3` at the beginning of the step. The lipid mass is therefore conserved. Since the supersaturation is constant during the phase-space step, a time step that is too large for the kinetics consumes more free lipids than available: the simulation then stops with an error, and the time step must be reduced.

References
----------
.. [#mcgraw1997] \McGraw, R. (1997). Description of aerosol dynamics by the quadrature method of moments. Aerosol science and Technology, 27(2), 255-265. doi: `10.1080/02786829708965471 <https://doi.org/10.1080/02786829708965471>`_\.

.. [#marchisio2003] \Marchisio, D. L., Pikturna, J. T., Fox, R. O., Vigil, R. D., & Barresi, A. A. (2003). Quadrature method of moments for population‐balance equations. AIChE Journal, 49(5), 1266-1276. doi: `10.1002/aic.690490517 <https://doi.org/10.1002/aic.690490517>`_\.

.. [#limpert2001] \Limpert, E., Stahel, W. A., & Abbt, M. (2001). Log-normal distributions across the sciences: keys and clues. BioScience, 51(5), 341-352. doi: `10.1641/0006-3568(2001)051[0341:LNDATS]2.0.CO;2 <https://doi.org/10.1641/0006-3568(2001)051[0341:LNDATS]2.0.CO;2>`_\.

.. [#marchisio2013] \Marchisio, D. L., & Fox, R. O. (2013). Computational models for polydisperse particulate and multiphase systems. Cambridge University Press. doi: `10.1017/CBO9781139016599 <https://doi.org/10.1017/CBO9781139016599>`_\.

.. [#marchisio2003b] \Marchisio, D. L., Vigil, R. D., & Fox, R. O. (2003). Quadrature method of moments for aggregation–breakage processes. Journal of Colloid and Interface Science, 258(2), 322-334. doi: `10.1016/S0021-9797(02)00054-1 <https://doi.org/10.1016/S0021-9797(02)00054-1>`_\.

.. [#friedlander2000] \Friedlander, S. K. (2000). Smoke, dust, and haze: Fundamentals of aerosol dynamics (2nd ed.). Oxford University Press.

.. [#shin2025] Shin, Devos, Udepurkar, Inguva, Myerson and Braatz, "Mechanistic modeling of lipid nanoparticle (LNP) precipitation via population balance equations (PBEs)", Chemical Engineering Journal 523 (2025) 167786.