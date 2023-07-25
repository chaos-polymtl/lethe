===========================================
The Incompressible Navier-Stokes Equations
===========================================

Lethe was designed to solve the incompressible Navier-Stokes equations. For the physics in Lethe, we generally try to follow the convention of the `BSL <https://en.wikipedia.org/wiki/Transport_Phenomena_(book)>`_. Assuming a constant density :math:`\rho` and kinematic viscosity :math:`\nu`, these equations take the following form:

.. math::
    \nabla \cdot \mathbf{u} &= 0   \\
    \frac{\partial \mathbf{u}}{\partial t}  + \mathbf{u} \cdot \nabla \mathbf{u} &= -\frac{1}{\rho} \nabla p^*  + \nu \nabla^2 \mathbf{u} +\mathbf{f}


where:

* :math:`\mathbf{u}` is the velocity of the fluid. :math:`\mathbf{u}` is a vector such that :math:`\mathbf{u}=[u,v]^T` in 2D and :math:`\mathbf{u}=[u,v,w]^T` in 3D;

* :math:`p^*` is the pressure;

* :math:`\nabla` is the `del operator <https://en.wikipedia.org/wiki/Del>`_;

* :math:`\rho` is the density of the fluid;

* :math:`\nu` is the `kinematic viscosity <https://en.wikipedia.org/wiki/Viscosity>`_ of the fluid;

* :math:`\mathbf{f}` is a momentum source term.

In `Einstein notation <https://en.wikipedia.org/wiki/Einstein_notation>`_, they can be written:


.. math::
    \partial_l u_l &= 0 

    \partial_t u_k + u_l \partial_l u_k &= -\frac{1}{\rho} \partial_k p^* + \nu \partial_l \partial_l u_k + f_k


In Lethe, instead of solving for pressure directly, we solve for a reduced pressure :math:`p=\frac{p^*}{\rho}`. This has two advantages, first the scaling of the pressure is more adequate with respect to that of the velocity and, secondly, for single phase simulations, this means that we require a single physical property: the kinematic viscosity. Consequently, unless specified otherwise, you should assume that pressure actually refers to the reduced pressure :math:`p`. Although it may seem confusing at first, this is standard in incompressible fluid dynamics.

.. note::
    The choice to use :math:`p^*` for the pressure and :math:`p` for the reduced pressure was made to lighten the writing throughout the documentation.

By default, Lethe does not use any units for Length (:math:`L`), Time (:math:`T`), Mass (:math:`M`) and Temperature (:math:`T`). Although this may seem confusing at first, this is relatively easy to work with. By default, most cases can be formulated using SI units (meters, seconds, kilograms and Celsius or Kelvin). However, in some cases, this may lead to a large imbalance in the scaling of different variables (e.g. the numerical value being 1000x that of the velocity). In these cases, the user must manually rescale the problem to a better balanced set of units (i.e. by using centimeter units instead of meters).

