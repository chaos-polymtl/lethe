================================
Heat Transfer Equations
================================

In Lethe, it is possible to solve heat transfer problems under various loading conditions. The main equation is derived from the energy equation in incompressible flows. Assuming a constant heat capacity :math:`C_p` and dynamic viscosity :math:`\mu`, the equation takes the following form: 

.. math::
    \rho C_p \frac{\partial T}{\partial t} + \rho C_p (\mathbf{u} \cdot \nabla)T - \nabla \cdot (k \nabla T) = \phi_v + Q

where:

* :math:`T` is the temperature;

* :math:`\mathbf{u}` is the velocity of the fluid;

* :math:`\nabla` is the `del operator <https://en.wikipedia.org/wiki/Del>`_;

* :math:`\rho` is the density;

* :math:`C_p` is the isobaric heat capacity;

* :math:`k` is the thermal conductivity;

* :math:`Q` is the energy source term or heat generation per unit of mass;

* :math:`\phi_v` is the viscous dissipation term. For an incompressible fluid it takes the following form: :math:`\phi_v = \mu (\nabla \mathbf{u} + \nabla \mathbf{u}^T):\nabla \mathbf{u}`, where :math:`\mu` is the `dynamic viscosity <https://en.wikipedia.org/wiki/Viscosity>`_;

Depending on the physics involved, the terms :math:`\phi_v` and :math:`Q` can be included or excluded and can take various definitions.

================================
Finite Element Formulation
================================

For the finite element formulation, we start from the strong form of the equation as shown above. We consider a domain :math:`\Omega` with boundary :math:`\Gamma`. The applied boundary conditions can vary depending on the problem being solved, but the main ones are the following:

* Dirichlet boundary condition: :math:`T = T_0` on :math:`\Gamma_D`;

* Neumann boundary condition: :math:`\mathbf{n} \cdot \nabla T = q` on :math:`\Gamma_N`;

* Robin boundary condition: :math:`\mathbf{n} \cdot \nabla T + h(T - T_{\infty}) = 0` on :math:`\Gamma_R`;

where :math:`\mathbf{n}` is the normal vector to the boundary, :math:`q` is the heat flux, :math:`h` is the heat transfer coefficient, and :math:`T_{\infty}` is the ambient temperature.

The weak form of the heat equation is obtained by multiplying the strong form by a test function :math:`v` and integrating over the domain :math:`\Omega`. After applying integration by parts and the Gauss-Ostrogradsky theorem, the weak form of the heat equation is given by the following equation:

.. math::
    \int_{\Omega} \left[ v \rho C_p \frac{\partial T}{\partial t} + v \rho C_p (\mathbf{u} \cdot \nabla)T + k (\nabla v \cdot \nabla T) \right] \;d\Omega - \int_{\Gamma} v k (\nabla T \cdot \mathbf{n}) \;d\Gamma = \int_{\Omega} v \left[ \phi_v + Q \right] \;d\Omega

Note that this formulation treats the thermal conductivity :math:`k` as a constant.

================================
Stabilization Techniques
================================

**Under construction**
