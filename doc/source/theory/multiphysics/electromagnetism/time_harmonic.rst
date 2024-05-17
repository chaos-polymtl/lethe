==================================
Time Harmonic Maxwell's Equations
==================================

In a future version of Lethe, there will be a module to solve time-harmonic Maxwell's equations which are used in microwave heating processes. The scope of this theory section is to give a brief description of where the equations come from, for further details, the reader is invited to look at the `Griffiths <https://en.wikipedia.org/wiki/Introduction_to_Electrodynamics (book)>`_. Thereupon, Maxwell's equations are a set of partial differential equations that form the foundation of classical electromagnetism. This famous set of equations can take multiple forms, but one of the most common ones is the macroscopic volumetric equations in matter :

.. math::
    \begin{align*}
    \nabla \cdot \mathbf{D} &= \rho_f \\
    \nabla \cdot \mathbf{B} &= 0 \\
    \nabla \times \mathbf{E} &= -\frac{\partial \mathbf{B}}{\partial t} \\
    \nabla \times \mathbf{H} &= \mathbf{J}_f + \frac{\partial \mathbf{D}}{\partial t}
    \end{align*}

where :

* :math:`\mathbf{D}` is the electric field flux density;
* :math:`\mathbf{B}` is the magnetic field flux density;
* :math:`\mathbf{E}` is the electric field intensity;
* :math:`\mathbf{H}` is the magnetic field intensity;
* :math:`\rho_f` is the free electric charge density;
* :math:`\mathbf{J}_f` is the free current density, and :math:`t` is time.

For linear medium, the electromagnetic fluxes can be related to their intensity with the following constitutive relations and the electric current is proportional to the electric field intensity:

.. math::
    \begin{align*}
    \mathbf{D} &= \varepsilon \mathbf{E} \\
    \mathbf{B} &= \mu \mathbf{H}\\
    \mathbf{J}_f &= \sigma_e \mathbf{E} + \mathbf{J}_{ext}
    \end{align*}

where :math:`\varepsilon` is the permittivity, :math:`\mu` is the permeability, :math:`\sigma_e` is the conductivity of the medium and :math:`\mathbf{J}_{ext}` is the externally applied current density. When one is confronted with oscillating electromagnetic fields a simplification that is commonly made is to consider the fields and exciting current as time harmonic, i.e. they can be expressed as :

.. math::
    \begin{align*}
    \mathbf{E}(\mathbf{x},t) &= \Re{\{\mathbf{E}(\mathbf{x}) e^{i\omega t}\}},\\
    \mathbf{H}(\mathbf{x},t) &= \Re{\{\mathbf{H}(\mathbf{x}) e^{i\omega t}\}},\\
    \mathbf{J}_{ext}(\mathbf{x},t) &= \Re{\{\mathbf{J}_{ext}(\mathbf{x}) e^{i\omega t}\}},
    \end{align*}

where :math:`\omega` is the angular frequency of the oscillating field. By substituting these expressions in Maxwell's equations presented above and cleverly combining them, one can obtain the following set of equations :

.. math::
    \begin{align*}
    \nabla \times \left( \frac{1}{\mu} \nabla \times \mathbf{E} \right) -\omega^2 \varepsilon_{eff} \mathbf{E} &= -i \omega \mathbf{J}_{ext},\\
    \nabla \times \left( \frac{1}{\varepsilon_{eff}} \nabla \times \mathbf{H} \right) - \omega^2 \mu \mathbf{H} &= \nabla \times \frac{\mathbf{J}_{ext}}{\varepsilon_{eff}},
    \end{align*}

where :math:`\varepsilon_{eff} = \varepsilon - i \frac{\sigma_e}{\omega}` is the effective permittivity of the medium. These equations are the time-harmonic Maxwell's equations for spatially varying permittivity and permeability, which can be order 2 tensor when the medium is anisotropic. Note that even if the harmonic oscillation is assumed, it is not a restrictive simplification. Indeed, any signals can be obtained by a summation of harmonic frequencies, and the time dependence of the electromagnetic field, if it is a quantity of interest, can be obtained by performing an inverse Fourier transform.