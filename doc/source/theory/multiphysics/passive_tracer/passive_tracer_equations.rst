=========================
Passive Tracer Equations
=========================

Lethe allows the simulation of a passive tracer, which is governed by the following advection-diffusion equation:

.. math::
    \frac{\partial C}{\partial t} + \mathbf{u} \cdot \nabla C - D \nabla^2 C + f = 0,

where:

* :math:`C` is the tracer concentration;
* :math:`\mathbf{u}` is the velocity of the fluid;
* :math:`\nabla` is the `del operator <https://en.wikipedia.org/wiki/Del>`_;
* :math:`D` is the diffusivity coefficient;
* :math:`f` is the tracer source term per unit of mass.

.. important::
    The tracer is advected by the fluid velocity, to which can be added a `drift velocity <../../../parameters/cfd/tracer_drift_velocity.html>`_. A reactive tracer feature also allows the prescription of a tracer consumption rate.

To obtain the finite element discretization, we start from the strong form previously presented, multiply it by a test function :math:`q`, then integrate it over the domain :math:`\Omega` with boundary :math:`\Gamma` such that:

.. math::
    \int_{\Omega} q \: \left( \frac{\partial C}{\partial t} + \mathbf{u} \cdot \nabla C - D \nabla^2 C - f \right) \: d\Omega = 0.

To obtain the corresponding weak form, we then apply integration by parts on the Laplacian term:

.. math::
    F(C) := (q, \partial_t C)_{\Omega} + (q, \mathbf{u} \cdot \nabla C)_{\Omega} + (\nabla q, D \: \nabla C)_{\Omega} - (q, f)_{\Omega} = 0.

To ensure that the integral is well defined in the domain :math:`\Omega`, we use appropriate solution spaces:

.. math::
    C \in \mathcal{C} (\Omega) = \mathcal{H}^1 (\Omega), \qquad q \in \mathcal{Q} (\Omega) = \mathcal{H}^1 (\Omega),

noting that this solution space for :math:`q` is suitable for a Continuous Galerkin formulation. If a Discontinuous Galerkin framework is assumed, the appropriate solution is defined in :math:`\mathcal{L}^2 (\Omega)`. 
Both DG and CG tracer formulations are available in Lethe (see how to assign them in the `FEM parameters <../../../parameters/cfd/fem.html>`_ section).

The weak formulation is written as: find :math:`C \in \mathcal{C} (\Omega) \times (0, T]` such that 

.. math::
    (q, \partial_t C)_{\Omega} + (q, \mathbf{u} \cdot \nabla C)_{\Omega} + (\nabla q, D \: \nabla C)_{\Omega} - (q, f)_{\Omega} = 0 \qquad \forall \: q \in \mathcal{Q}.

========================
Stabilization Techniques
========================

Two stabilization terms are added to the weak form of the passive tracer equation:

.. math::
   F(C)_{\text{stab}} := F(C) + \underbrace{(\tau \: [\mathbf{u} \cdot \nabla q], \: \mathcal{R}(C))}_{a_{\text{GLS}}} + \underbrace{(v_{\text{DCDD}} \nabla q, \: \mathbf{d}_{\text{DCDD}} \cdot \nabla C)}_{a_{\text{DCDD}}} ,

In the GLS stabilization term :math:`a_{\text{GLS}}`, :math:`\mathcal{R}(C)` is the strong residual of the governing equation, and the stabilization parameter :math:`\tau_k` is defined as in `Tezduyar (1992) <https://doi.org/10.1016/0045-7825(92)90141-6>`_:

.. math::

   \tau = \left[ \left( \frac{1}{\Delta t} \right)^{2} + \left( \frac{2 |\mathrm{u}|}{h_{\text{conv}}} \right)^{2} + 9 \left( \frac{4 \nu}{h^2_{\text{diff}}} \right)^{2} \right]^{-1/2},

where :math:`\Delta t` is the time step, and :math:`h_{\text{conv}}` and :math:`h_{\text{diff}}` are the size of the element related to the convection transport and diffusion mechanism, respectively. In Lethe, both element sizes are set to the diameter of a sphere having a volume equivalent to that of the cell. 

.. seealso::
    The same definition of the stabilization parameter is used in the `Fluid Dynamics <../fluid_dynamics/stabilization.html>`_ solver.

The second stabilization term is a Discontinuity-Capturing Directional Dissipation (DCDD) shock-capturing scheme proposed by `Tezduyar (2003) <https://doi.org/10.1002/fld.505>`_ that aims to deal with crosswind oscillations. 

.. seealso::
    The DCDD shock-capturing scheme is also used in the `CLS <../../multiphase/cfd/cls.html>`_ and `Heat Transfer <../heat_transfer/heat_transfer.html>`_ modules. The CLS DCDD implementation is detailed in `this publication <https://doi.org/10.1016/j.cpc.2025.109880>`_, and it is analogous to the tracer module.
