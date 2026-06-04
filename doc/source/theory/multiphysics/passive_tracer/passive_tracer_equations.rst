=========================
Passive Tracer Equations
=========================

Lethe allows the simulation of a passive tracer, which is governed by the equation:

.. math::
    \frac{\partial T}{\partial t} + \mathbf{u} \cdot \nabla T - \kappa \nabla^2 T + S = 0

where:

* :math:`T` is the tracer concentration;
* :math:`\mathbf{u}` is the velocity of the fluid;
* :math:`\nabla` is the `del operator <https://en.wikipedia.org/wiki/Del>`_;
* :math:`\kappa` is the diffusivity coefficient;
* :math:`S` is the tracer source term per unit of mass.

==========================
Finite Element Formulation
==========================

To obtain the finite element discretization, we start from the strong form previously presented and integrate it over a domain :math:`\Omega` with boundary :math:`\Gamma` such that:

.. math::
    \int_{\Omega} \left( \frac{\partial T}{\partial t} + \mathbf{u} \cdot \nabla T - \kappa \nabla^2 T - S \right) d\Omega = 0

To obtain the corresponding weak form, we multiply the equation above by a test function :math:`q` and apply integration by parts on the Laplacian term:

.. math::
    \int_{\Omega} q \frac{\partial T}{\partial t} d\Omega + \int_{\Omega} q \mathbf{u} \cdot \nabla T d\Omega + \int_{\Omega} \nabla q \kappa \nabla T d\Omega - \int_{\Omega} q S d\Omega = 0

Using the discrete form of the tracer concentration field :math:`T = \sum_j T_j \phi_j` and the test function :math:`q = \phi_i`:

.. math::
    \int_{\Omega} \phi_i \frac{\partial \phi_j}{\partial t} d\Omega + \int_{\Omega} \phi_i \mathbf{u} \cdot \nabla \phi_j d\Omega + \int_{\Omega} \nabla \phi_i \kappa \nabla \phi_j d\Omega - \int_{\Omega} \phi_i S d\Omega = 0

========================
Stabilization Techniques
========================

The GLS stabilization term added is given by:

.. math::
    a_{GLS} = \int_{\Omega} \tau (\mathbf{u} \cdot \nabla q) \mathrm{R}(T),

where :math:`\mathrm{R}(T)` is the strong residual of the governing equation. The stabilization parameter definition herein used is the one presented by `Tezduyar (1992) <https://doi.org/10.1016/0045-7825(92)90141-6>`_ and presented in the `Fluid Dynamics <../fluid_dynamics/stabilization.html>`_ subsection.

====================
Shock-capturing term
====================

**Under construction**
