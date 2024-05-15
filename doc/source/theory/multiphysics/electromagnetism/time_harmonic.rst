==================================
Time Harmonic Maxwell's Equations
==================================

In future version of Lethe, there will be a module to solve time harmonic Maxwell's equations which are used in microwave heating processes. The scope of this theory section is to give a brief description of where the equations comes from, for further details, the reader is invited to look at the `Griffiths <https://en.wikipedia.org/wiki/Introduction_to_Electrodynamics (book)>`_. Thereupon, Maxwell's equations are a set of partial differential equations that form the foundation of classical electromagnetism. This famous set of equations can take multiple forms, but one of the most common one is the macroscopic volumetric equations in matter :

.. math::
    \nabla \cdot \mathbf{D} = \rho_f \\
    \nabla \cdot \mathbf{B} = 0 \\
    \nabla \times \mathbf{E} = -\frac{\partial \mathbf{B}}{\partial t} \\
    \nabla \times \mathbf{H} = \mathbf{J}_f + \frac{\partial \mathbf{D}}{\partial t}

where :

* :math:`\mathbf{D}` is the electric field flux density;
* :math:`\mathbf{B}` is the magnetic field flux density;
* :math:`\mathbf{E}` is the electric field intensity;
* :math:`\mathbf{H}` is the magnetic field intensity;
* :math:`\rho_f` is the free electric charge density;
* :math:`\mathbf{J}_f` is the free current density, and :math:`t` is time.

For linear medium, the electromagnetic fluxes can be related to their intensity with the following constitutive relations and the electric current is proportional to the electric field intensity:

.. math::
    \mathbf{D} = \varepsilon \mathbf{E} \\
    \mathbf{B} = \mu \mathbf{H}\\
    \mathbf{J}_f = \sigma_e \mathbf{E}

where :math:`\varepsilon` is the permittivity, :math:`\mu` is the permeability and :math:`\sigma_e` is the conductivity of the medium. When one is confronted to oscillating electromagnetic fields a simplification that is commonly made is to consider the fields as time harmonic, i.e. they can be expressed as :

.. math::
    \mathbf{E}(\mathbf{x},t) = \Re{\{\mathbf{E}(\mathbf{x}) e^{i\omega t}\}},\\
    \mathbf{H}(\mathbf{x},t) = \Re{\{\mathbf{H}(\mathbf{x}) e^{i\omega t}\}},

where :math:`\omega` is the angular frequency of the oscillating field. By substituting these expressions in the Maxwell's equations presented above and combining them in a clever way, one can obtain the following set of equations :

.. math::
    \nabla \times \left( \frac{1}{\mu} \nabla \times \mathbf{E} \right) -\omega^2 \epsilon_{eff} \mathbf{E} = 0,\\
    \nabla \times \left( \frac{1}{\epsilon_{eff}} \nabla \times \mathbf{H} \right) - \omega^2 \mu \mathbf{H} = 0,

where :math:`\epsilon_{eff} = \epsilon - i \frac{\sigma_e}{\omega}` is the effective permittivity of the medium. These equations are the time harmonic Maxwell's equations for spatially varying permittivity and permeability, which can be order 2 tensor when the medium is anisotropic. Note that even if the harmonic oscillation is assumed, it is not a restrictive simplification. Indeed, any signals can be obtain by a summation of harmonic frequencies and the time dependence of the electromagnetic field, if it is a quantity of interest, can be obtain by performing a Fourier transform.