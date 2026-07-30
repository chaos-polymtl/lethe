==================================
Time Harmonic Maxwell's Equations
==================================

Lethe now features a solver for the time-harmonic Maxwell's equations that can be coupled to the heat transfer, enabling the simulation of microwave heating processes. The scope of this theory section is to give a brief description of where the equations come from, for further details, the reader is invited to look at the `Griffiths <https://en.wikipedia.org/wiki/Introduction_to_Electrodynamics (book)>`_ book for an introduction on electromagnetism theory and on the following `article <https://ieeexplore.ieee.org/document/9682427>`_ for a review on computational electromagnetism. 

Maxwell's equations are a set of partial differential equations that form the foundation of classical electromagnetism. This famous set of equations can take multiple forms, but one of the most common ones is the macroscopic volumetric equations in matter :

.. math::
    \begin{align*}
    \nabla \cdot \mathbf{D} &= \rho_f \\
    \nabla \cdot \mathbf{B} &= 0 \\
    \nabla \times \mathbf{E} &= -\frac{\partial \mathbf{B}}{\partial t} \\
    \nabla \times \mathbf{H} &= \mathbf{J}_f + \frac{\partial \mathbf{D}}{\partial t}
    \end{align*}

where:

* :math:`\mathbf{D}` is the electric field flux density;
* :math:`\mathbf{B}` is the magnetic field flux density;
* :math:`\mathbf{E}` is the electric field intensity;
* :math:`\mathbf{H}` is the magnetic field intensity;
* :math:`\rho_f` is the free electric charge density;
* :math:`\mathbf{J}_f` is the free current density, and :math:`t` is time.

For linear medium, the electromagnetic fluxes and current can be related to their intensity with the following constitutive relations:

.. math::
    \begin{align*}
    \mathbf{D} &= \varepsilon_{\mathrm{em}} \mathbf{E} \\
    \mathbf{B} &= \mu_{\mathrm{em}} \mathbf{H}\\
    \mathbf{J}_f &= \sigma_e \mathbf{E} + \mathbf{J}_{\mathrm{ext}}
    \end{align*}
    
where :math:`\varepsilon_{\mathrm{em}}` is the permittivity, :math:`\mu_{\mathrm{em}}` is the permeability, :math:`\sigma_e` is the conductivity of the medium and :math:`\mathbf{J}_{\mathrm{ext}}` is the externally applied current density. When one is confronted with oscillating electromagnetic fields, a simplification that is commonly made is to consider the fields and exciting currents as time harmonic, `i.e.` they can be expressed as [1]_:

.. math::
    \begin{align*}
    \mathbf{E}(\mathbf{x},t) &= \Re{\{\mathbf{E}_{\mathrm{spatial}}(\mathbf{x}) e^{-i\omega t}\}},\\
    \mathbf{H}(\mathbf{x},t) &= \Re{\{\mathbf{H}_{\mathrm{spatial}}(\mathbf{x}) e^{-i\omega t}\}},\\
    \mathbf{J}_{\mathrm{ext}}(\mathbf{x},t) &= \Re{\left\{\mathbf{J}_{\mathrm{ext, spatial}}(\mathbf{x}) e^{-i\omega t}\right\}},
    \end{align*}
    
where :math:`\omega` is the angular frequency of the oscillating fields and :math:`t` is the time.

By substituting these expressions in Maxwell's equations presented above and cleverly combining them, one can obtain the either ones of the following equations:

.. math::
    \begin{align*}
    \nabla \times \left( \frac{1}{\mu_{\mathrm{em}}} \nabla \times \mathbf{E} \right) -\omega^2 \varepsilon_{\mathrm{em},{\mathrm{eff}}} \mathbf{E} &= i \omega \mathbf{J}_{\mathrm{ext}},\\
    \nabla \times \left( \frac{1}{\varepsilon_{\mathrm{em},{\mathrm{eff}}}} \nabla \times \mathbf{H} \right) - \omega^2 \mu_{\mathrm{em}} \mathbf{H} &= \nabla \times \frac{\mathbf{J}_{\mathrm{ext}}}{\varepsilon_{\mathrm{em},{\mathrm{eff}}}},
    \end{align*}
    
where :math:`\varepsilon_{\mathrm{em},{\mathrm{eff}}} = \varepsilon_{\mathrm{em}} + i \frac{\sigma_e}{\omega}` is the effective permittivity of the medium. These equations are the time-harmonic Maxwell's equations for spatially varying permittivity and permeability, which can be rank 2 tensors when the medium is anisotropic. Note that even if the harmonic oscillation is assumed, it is not a restrictive simplification. Indeed, any signal can be obtained by a summation of harmonic frequencies, and the time dependence of the electromagnetic field, if it is a quantity of interest, can be obtained by performing an inverse Fourier transform. 

.. note::
    In Lethe, however, we only support solving for a single frequency and the reconstruction of a multi-frequency signal or its time dependence needs to be done externally. 

For numerical purposes, in Lethe, the time-harmonic Maxwell's equations are solved in a dimensionless coupled system that takes the following form:

.. math::
    \begin{align*}
    \nabla \times \mathbf{E} -i\omega \mu_r \mathbf{H} &= 0,\\
    \nabla \times \mathbf{H} + i\omega \varepsilon_{r,\mathrm{eff}} \mathbf{E} &=  \mathbf{J}_{\mathrm{ext}},
    \end{align*}
    
where :math:`\mu_r` is the relative permeability and :math:`\varepsilon_{r,\mathrm{eff}}` is the effective relative permittivity of the medium. 

.. [1] Lethe uses the physicist convention for the wave propagation, which describes a traveling wave as :math:`e^{\mathbf{k}\cdot \mathbf{x} - i\omega t}`. This has some implications on the sign at multiple places in the equations, and the user should be aware of this when comparing with other references.