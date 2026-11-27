.. _qmom-lipid-nanoparticle-production:

=============================
Lipid Nanoparticle Production
=============================

With the lipid nanoparticle production process, the nucleation, growth and coalescence are given by the models of Shin et al. [#shin2025]_, whose parameters are in the ``lipid nanoparticle production`` subsection of :doc:`../../../parameters/cfd/qmom`: the generic growth, nucleation and coalescence models are ignored, but each mechanism is still enabled by its ``enable`` parameter. The two scalars are the ethanol concentration :math:`c_e` and the free lipid concentration :math:`C`, both in :math:`\text{kg} \cdot \text{m}^{-3}`.

Solubility and Supersaturation
------------------------------

At every degree of freedom, the organic (ethanol) volume fraction is :math:`\phi_o = c_e / \rho_o`, clipped to :math:`[0, 1]`, and the aqueous volume fraction is :math:`\phi_a = 1 - \phi_o`. The mole fraction solubility :math:`x` of the lipid in the mixture, its solubility :math:`C_\infty` and the supersaturation :math:`S` are:

.. math::

  \ln x = \phi_a \ln x_a + \phi_o \ln x_o, \quad C_\infty = \frac{M_l x}{1 - x} \left(\frac{\phi_a \rho_a}{M_a} + \frac{\phi_o \rho_o}{M_o}\right), \quad S = \frac{\max(C, 0)}{C_\infty}

where:

- :math:`x` is the mole fraction solubility of the lipid in the mixture;
- :math:`x_a` and :math:`x_o` are the mole fraction solubilities of the lipid in the aqueous and organic phases;
- :math:`\rho_a` and :math:`\rho_o` are the densities of the aqueous and organic phases (:math:`\text{kg} \cdot \text{m}^{-3}`);
- :math:`M_a` and :math:`M_o` are the molar masses of the aqueous and organic phases (:math:`\text{kg} \cdot \text{mol}^{-1}`);
- :math:`M_l` is the molar mass of the lipid (:math:`\text{kg} \cdot \text{mol}^{-1}`);
- :math:`C_\infty` is the solubility of the lipid in the mixture (:math:`\text{kg} \cdot \text{m}^{-3}`);
- :math:`C` is the free lipid concentration (:math:`\text{kg} \cdot \text{m}^{-3}`);
- :math:`S` is the supersaturation.

The capillary length is :math:`\alpha = 4 \gamma v_m / (k_B T)`, where :math:`\gamma` is the surface tension, :math:`v_m` the molecular volume of the lipid, :math:`k_B` the Boltzmann constant and :math:`T` the operating temperature.

Growth
------

The growth includes Ostwald ripening: particles smaller than their equilibrium size dissolve.

.. math::

  G = k_g \xi C_\infty \left(S - e^{\alpha / \xi}\right)

where:

- :math:`G` is the growth rate;
- :math:`k_g` is the growth factor;
- :math:`\xi` is the particle size;
- :math:`C_\infty` is the solubility of the lipid in the mixture;
- :math:`S` is the supersaturation;
- :math:`\alpha` is the capillary length.

Nucleation
----------

The nucleation only happens in a supersaturated mixture (:math:`S > 1`), at the rate:

.. math::

  J = A \exp\left(-\frac{16 \pi \gamma^3 v_m^2}{3 (k_B T)^3 (\ln S)^2}\right)

where:

- :math:`J` is the nucleation rate;
- :math:`A` is the nucleation factor;
- :math:`\gamma` is the surface tension;
- :math:`v_m` is the molecular volume of the lipid;
- :math:`k_B` is the Boltzmann constant;
- :math:`T` is the operating temperature;
- :math:`S` is the supersaturation.

The sizes of the nuclei follow a log-normal distribution whose standard deviation is the nuclei standard deviation and whose mean is the critical size:

.. math::

  L_c = \frac{\alpha}{\ln S}

where:

- :math:`L_c` is the critical size;
- :math:`\alpha` is the capillary length;
- :math:`S` is the supersaturation.

Coalescence
-----------

The coalescence kernel is the Brownian kernel times the attachment efficiency of the Fuchs model, limited by the PEG factor :math:`k_{PEG}`:

.. math::

  \beta = \frac{2 k_B T (\xi_1 + \xi_2)^2}{3 \mu \xi_1 \xi_2} \min\left(\frac{k_{PEG}}{W}, k_{PEG}\right), \quad W = r_s \int_0^\infty \frac{\exp\left(\Phi(a) / (k_B T)\right)}{(a + r_s)^2} \mathrm{d} a

where:

- :math:`\beta` is the coalescence kernel;
- :math:`\xi_1` and :math:`\xi_2` are the sizes of the two colliding particles;
- :math:`k_B` is the Boltzmann constant;
- :math:`T` is the operating temperature;
- :math:`\mu` is the viscosity of the mixture (:math:`\text{Pa} \cdot \text{s}`);
- :math:`k_{PEG}` is the PEG factor;
- :math:`W` is the Fuchs stability ratio, whose inverse is the attachment efficiency;
- :math:`r_s` is the sum of the radii of the particles;
- :math:`a` is the distance between the surfaces of the particles;
- :math:`\Phi` is the interaction potential, the sum of the van der Waals attraction (Hamaker constant) and of the electrostatic repulsion (surface potential) of the particles.

The Debye length uses the ion concentration of the mixture and the relative permittivity of the local mixture, :math:`\phi_a \varepsilon_a + \phi_o \varepsilon_o`, where :math:`\varepsilon_a` and :math:`\varepsilon_o` are the relative permittivities of the aqueous and organic phases. The integral is computed with the Simpson rule on 400 logarithmic intervals of :math:`a` between :math:`10^{-13}` and :math:`10^{-5}\ \text{m}`.

Consumption of the Free Lipids
------------------------------

The particles nucleation and growth consume and release free lipids; at the end of the phase-space step, the free lipid concentration of every degree of freedom loses the mass of lipid gained by its particles (end of step 2 of :ref:`qmom-operator-splitting`):

.. math::

  C \leftarrow C - \rho_l k_v \left(M_3 - M_3^0\right)

where:

- :math:`C` is the free lipid concentration, at the end of its transport in physical space before the update;
- :math:`\rho_l` is the density of the lipid in the particles (:math:`\text{kg} \cdot \text{m}^{-3}`);
- :math:`k_v` is the volume shape factor;
- :math:`M_3` is the third moment of the size distribution at the end of the step;
- :math:`M_3^0` is the moment :math:`M_3` at the beginning of the phase-space step, i.e. at the end of its transport in physical space.

The lipid mass is therefore conserved. 

.. note::
    Since the supersaturation is constant during the phase-space step, a time step that is too large for the kinetics consumes more free lipids than available: the simulation then stops with an error, and the time step must be reduced.

References
----------
.. [#shin2025] Shin, Devos, Udepurkar, Inguva, Myerson and Braatz, "Mechanistic modeling of lipid nanoparticle (LNP) precipitation via population balance equations (PBEs)", Chemical Engineering Journal 523 (2025) 167786. doi: `10.1016/j.cej.2025.167786 <https://doi.org/10.1016/j.cej.2025.167786>`_