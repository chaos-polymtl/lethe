The Incompressible Navier-Stokes equations
#############################################

Lethe was designed to solve the incompressible Navier-Stokes equations. For the physics in Lethe, we generally try to follow the convention of the `BSL <https://en.wikipedia.org/wiki/Transport_Phenomena_(book)>`_ Assuming a constant density :math:`\rho` and kinematic viscosity :math:`\nu`, these equations take the following form:

.. math::
    \nabla \cdot \mathbf{u} &= 0   \\
    \frac{\partial \mathbf{u}}{\partial t}  + \mathbf{u} \cdot \nabla \mathbf{u} &= -\frac{1}{\rho} \nabla p  - \nu \nabla^2 \mathbf{u} 


where:

* :math:`\mathbf{u}` is the velocity of the fluid. :math:`\mathbf{u}` is a vector such that :math:`\mathbf{u}=[u,v]^T` in 2D and :math:`\mathbf{u}=[u,v,w]^T` in 3D.

* :math:`p` is the pressure

* :math:`\nabla` is the `del operator <https://en.wikipedia.org/wiki/Del>`_

* :math:`\rho` is the density of the fluid.

* :math:`\nu` is the `kinematic viscosity <https://en.wikipedia.org/wiki/Viscosity>`_ of the fluid.

In `Einstein notation <https://en.wikipedia.org/wiki/Einstein_notation>`_, they can be written:

.. math::
    \partial_j u_j &= 0  \\
    \partial_i+ u_i \partial_i u_j &= -\frac{1}{p} \partial_i p + \nu \partial_j \partial_j u_i


.. note::
    We try to use Einstein's notation to write the majority of the partial differential equations in this documentation. We have found that this way is less error-prone when higher-order tensors are concerned. 


In Lethe, instead of solving for pressure directly, we solve for a reduced pressure :math:`p^*=\frac{p}{\rho}`. This has two advantages, first the scaling of the pressure is more adequate with respect to that of the velocity and, secondly, for single phase simulations, this means that we require a single physical property: the kinematic viscosity. Consequently, unless specified otherwise, you should assume that pressure actually refers to the reduced pressure :math:`p^*`. Although may seem confusing at first, this is standard in incompressible fluid dynamics.

By default, Lethe does not use any units for Length (:math:`L`), Time (:math:`T`), Mass (:math:`M`) and Temperature (:math:`T`). Although this may seem confusing at first, this is relatively easy to work with. By default, most cases can be formulated using SI units (meters, seconds, kilograms and Celsius or Kelvin). However, in some cases, this may lead to a large imbalance in the scaling of different variables (e.g. the numerical value being 1000x that of the velocity). In these cases, the user must manually rescale the problem to a set of unit which are better balanced.



