==================================
Rheology and Non-Newtonian Fluids
==================================

For numerous flows that we encounter in chemical engineering, Newton's equation for viscosity is insufficient to describe the rheological behavior of the fluids. For these problems, the viscosity is not constant anymore. Lethe supports the simulation of non-Newtonian fluids by using generalized Newtonian models. In these models, the viscosity depends on the magnitude of the shear-rate tensor :math:`\dot{\gamma}`. In this case, the viscosity cannot be assumed to be constant anymore and the incompressible Navier-Stokes equations take the following form:

.. math::
    \partial_l u_l &= 0  \\
    \partial_t u_k + u_l \partial_l u_k &= -\frac{1}{\rho} \partial_k p  -  \partial_l \mathbf{\tau}_{lk} + f_k

where :math:`\mathbf{\tau}_{ji}` is the deviatoric stress tensor 

.. math::
    \mathbf{\tau}_{lk} = - \eta \left( \partial_l u_k + \partial_k u_l \right)

The fact that the viscosity may now depend on the shear rate introduces additional non-linearity to the underlying Navier-Stokes equations.

Lethe supports two rheology model:

* The Power-law model

.. math::

  \eta(\dot{\gamma}) = K \dot{\gamma}^{n-1}


where :math:`\eta` is the **kinematic viscosity**, :math:`\dot{\gamma}` is the shear rate, :math:`K` is the consistency index and :math:`n` is the power exponent.

* The Carreau 5 parameter model:

.. math::

  \eta(\dot{\gamma}) =\eta_{\infty} + (\eta_0 - \eta_{\infty}) \left[ 1 + (\dot{\gamma}\lambda)^a\right]^{\frac{n-1}{a}}
  
where :math:`\eta_\infty` is the viscosity in the high-shear rate plateau, :math:`\eta_0` is the zero-shear viscosity, :math:`a` is the Carreau parameter (generally set to 2), :math:`\lambda` is the relaxation time associated to the fluid and :math:`n` is the power exponent. 


.. note::
    Numerically, the Power-law model is ill-posed because it presumes an infinite viscosity in regions of no-shear. Consequently, we generally advise the users to use the Carreau model since it is much more robust.
