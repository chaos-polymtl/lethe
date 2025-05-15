================================
Cahn-Hilliard Method
================================

The Cahn-Hilliard system of equations [#cahn1958]_ is a model used to describe the process of phase separation based on the principle of free energy minimization. The key idea at the heart of the Cahn-Hilliard equation is that the system try to achieve a state where the free energy is minimized while competing with the energy cost associated with creating new interfaces between phases. Let us introduce those concepts formally.

Let :math:`\Omega = \Omega_0 \cup \Omega_1` be the domain formed by two fluids, namely fluid :math:`0` and :math:`1`, with :math:`\Gamma` the boundaries of the system. Like in :doc:`vof`, we define a scalar function :math:`\phi` as a phase indicator such that:

.. math::
  \phi =
  \begin{cases}
    1 \phantom{-} \quad \forall \mathbf{x} \in \Omega_0\\
    -1 \quad \forall \mathbf{x} \in \Omega_1
  \end{cases}
  
The phase indicator transitions smoothly from one extremum to the other with the shape of an hyperbolic tangent function as illustrated below: 

.. image:: images/tanh-solution.png
    :alt: Schematic
    :align: center
    :width: 600
    
The length required to go from :math:`\phi=-0.99` to :math:`\phi=0.99` is about 7.5 times :math:`\varepsilon`, which we will see appear later in the equations.

Let us introduce the free energy functional :math:`\mathcal{F}`:

.. math::
  \mathcal{F} = \int_{\Omega} \Psi \mathrm{d}\Omega
  
:math:`\Psi` corresponds to the free energy density. Its expression is as follows:

.. math::
  \Psi = \lambda\left(F(\phi) + \frac{1}{2}|\nabla \phi|^2\right)
  
It is decomposed in a bulk free energy :math:`F(\phi)` and an interface energy :math:`\frac{1}{2}|\nabla \phi|^2`. :math:`\lambda` is called the mixing energy though not dimensionally consistent to an energy. It can be understood as a scaling factor for the energy of the system and it is proportionnal to the surface tension coefficient. The bulk free energy has a double-well form, its expression is:

.. math::
  F(\phi) = \frac{(1-\phi^2)^2}{4\varepsilon^2}

with :math:`\varepsilon` a parameter related to the thickness of the interface between the phases.
  
To formalize the idea that the system tries to lower its free energy, we introduce a new variable, :math:`\eta`, called the chemical potential. It is defined as the variational derivative of :math:`\mathcal{F}`:

.. math::
  \eta = \frac{\delta\mathcal{F}}{\delta\phi} = \lambda\left(-\frac{\phi(1-\phi^2)}{\varepsilon^2} + \nabla^2\phi\right)
  
  
Then, the phases have to move to satisfy free energy minimization, by going from high chemical potential regions to low chemical potential regions. Let us introduce the flux of phase due to chemical potential differences, denoted by :math:`\mathbf{J}`:

.. math::
  \mathbf{J} = -M(\phi)\nabla\eta
   
With :math:`M(\phi)` the mobility function. Then, we can write the following continuity equation, which also takes into account an external velocity field :math:`\mathbf{u}`:

.. math::
  \frac{\partial \phi}{\partial t} + \mathbf{u}\cdot \nabla \phi = \nabla \cdot (M(\phi)\nabla \eta)
  
To demonstrate that the stable solutions of the Cahn-Hilliard have naturally diffuse interfaces, let us look at a simplified 1D case in equilibrium: :math:`\eta = 0` in the chemical potential equation. We obtain a differential equation on :math:`\phi` whose solution is:

.. math::
  \phi(x) = -\tanh{\left(\frac{x}{\sqrt{2}\varepsilon}\right)}
  
  
Lastly, we need to relate the surface tension coefficient with the mixing energy. Cahn and Hilliard define surface tension as the excess free energy of the system due to the interface, the second term in the expression of :math:`\Psi`. If we consider a system in equilibrium, and suppose a planar interface, we can write:

.. math::
  \sigma = \int_{-\infty}^{+\infty}\lambda \left(\frac{\mathrm{d}\phi}{\mathrm{d}x}\right)^2 \mathrm{d}x
  
Substituting the solution previously found, we obtain the following expression for :math:`\lambda`:

.. math::
  \lambda = \frac{3\varepsilon\sigma}{2\sqrt{2}}
  
For the problem to have a unique solution, we give the following no-flux boundary conditions on :math:`\partial \Omega` for the phase field and chemical potential:

.. math::
  (\nabla \phi \cdot\mathbf{n})_{| \partial \Omega} = 0
  
  (\nabla \eta \cdot \mathbf{n})_{| \partial \Omega} = 0
  
Finite Element Formulation
---------------------------

Let us write the weak formulation. Let :math:`\alpha` and :math:`\beta` be the scalar test functions associated to :math:`\phi` and :math:`\eta`. Let us first introduce the function spaces used to ensure the integrals exist:  

.. math::

  \begin{align}
  & (\phi, \eta) \in \psi(\Omega) = (H^1(\Omega)\times [0,T])\times (H^1(\Omega)\times [0,T])\\
  & (\alpha, \beta) \in \xi(\Omega) = (L^2(\Omega)\times [0,T]) \times (L^2(\Omega)\times [0,T])\\
  \end{align}


We multiply each equation by their test function and integrate over :math:`\Omega`:

.. math::
  \int_\Omega \alpha\left(\frac{\partial \phi}{\partial t} + \mathbf{a}\cdot \nabla \phi -\nabla \cdot (M(\phi)\nabla \eta)\right) \mathrm{d}\Omega = 0
  
  \int_\Omega \beta\lambda\left(\frac{\phi(1-\phi^2)}{\varepsilon^2} - \nabla^2\phi\right)\mathrm{d}\Omega = 0
  
After using the integration by part and Green-Ostrogradski's theorem:

.. math::
  \int_\Omega \alpha\left(\frac{\partial \phi}{\partial t} + \mathbf{a}\cdot \nabla \phi\right)\mathrm{d}\Omega -\int_\Omega M(\phi) \nabla\alpha \cdot\nabla\eta \mathrm{d}\Omega + \cancelto{\mathrm{no-flux}}{\int_{\Gamma} M(\phi)\alpha \nabla \eta \cdot \mathbf{n} \mathrm{d}\Gamma} = 0
  
  \int_\Omega \beta\lambda\left(\frac{\phi(1-\phi^2)}{\varepsilon^2}\right)\mathrm{d}\Omega - \int_\Omega \nabla \beta \cdot \lambda\nabla\phi\mathrm{d}\Omega + \cancelto{\mathrm{no-flux}}{\int_{\Gamma} \alpha \lambda\nabla \phi \cdot \mathbf{n} \mathrm{d}\Gamma} = 0
  
Using Petrov-Galerkin method, the finite element formulation reads:

Find :math:`(\phi^h,\eta^h) \in \psi^h` such that:
  
.. math::  
  \begin{array}{rl}
  \displaystyle \int_\Omega \alpha^h\left(\frac{\partial \phi^h}{\partial t} + \mathbf{a} \cdot \nabla \phi^h \right) \mathrm{d}\Omega - \int_\Omega M(\phi^h) \nabla \alpha^h \cdot \nabla \eta^h \mathrm{d}\Omega &= 0 \\[1em]
  \displaystyle \int_\Omega \beta^h \lambda \left( \frac{\phi^h(1-(\phi^h)^2)}{\varepsilon^2} \right) \mathrm{d}\Omega - \int_\Omega \lambda\nabla \beta^h \cdot \nabla \phi^h \mathrm{d}\Omega &= 0 
  \end{array}
  \quad \forall (\alpha^h, \beta^h) \in \xi^h
  
Stabilization
---------------------------
   
While developping the code, it turned useful to add a numerical diffusion term in the chemical potential form for some example. The new equation is:

.. math::
  \eta = \lambda\left(-\frac{\phi(1-\phi^2)}{\varepsilon^2} + \nabla^2\phi\right) - \xi h^2 \nabla^2 \eta = 0
  
With :math:`h` the local cell size and :math:`\xi` a user-defined smoothing coefficient (in general between 0 and 1). This fonctionnality may be deprecated later.

Coupling to the Navier-Stokes equations
----------------------------------------

Because of the presence of two fluids and the interface, two additional effects must be taken into account in the fluid dynamics equations. 
First, the surface tension forces will deform the interface to minimize the interface energy. The link between the phase field and surface tension force is given by the **Kortoweg stress tensor**:

.. math::
  \mathbf{T_K} = \lambda(\nabla \phi \otimes \nabla \phi) 
  
This tensor is added to the usual viscous stress tensor to take into account the capillary effects. The capillary forces are obtained by taking its divergence:

.. math::
  \begin{align}
   \mathbf{f_\sigma} & = \nabla \cdot (\lambda(\nabla \phi \otimes \nabla \phi))\\
  & = \eta\nabla\phi + \nabla\Psi
  \end{align}

We then define a modified pressure :math:`\hat{p}`, which corresponds to the usual pressure with the additional :math:`\Psi` term. This new pressure is the same in the bulk phases and varies more smoothly in the interface [#lovric2019]_. Then, to take into account the change of momentum of the system due to the diffusive flux of species, we add the following term into the :doc:`usual momentum equation<../../multiphysics/fluid_dynamics/navier-stokes>`:

.. math::
  (\mathbf{\tilde{J}}\cdot \nabla)\mathbf{v} = (\frac{\rho_0-\rho_1}{2}\mathbf{J}\cdot \nabla)\mathbf{v}
  
Finally, the local physical properties (density, viscosity, `etc`.)  are deduced from the phase field by taking a linear approximation:

.. math::
  \begin{align}
  &\rho(\phi) = \frac{1-\phi}{2}\rho_1 + \frac{1+\phi}{2}\rho_0 \\
  &\mu(\phi) = \frac{1-\phi}{2}\mu_1 + \frac{1+\phi}{2}\mu_0 \\
  \end{align}
  
The Cahn-Hilliard-Navier-Stokes momentum equation solved in Lethe is:

.. math::
  \begin{align}
  & \rho(\phi)\left(\frac{\partial\mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u}\right) + \left(\frac{\rho_0-\rho_1}{2}M(\phi)\nabla\eta\cdot \nabla\right)\mathbf{u}  \\
   & - \nabla \cdot \left(\mu(\phi)(\nabla\mathbf{u} + \nabla\mathbf{u}^\mathbf{T})\right) + \nabla \hat{p} - \eta\nabla\phi = 0 \\
  \end{align}
  
With an adequate choice of definition of velocity [#abels2011]_, the velocity field remains divergence-free:

.. math::
  \nabla \cdot \mathbf{u} = 0
  
However, the continuity equation is slightly different than the usual single-phase one:

.. math::
  \frac{\partial\rho}{\partial t} + \nabla \cdot (\rho\mathbf{u} +\mathbf{\tilde{J}}) = 0
  
This is to take into account the diffusion of species in the system.
  
References
-----------

`[1] <https://dx.doi.org/10.1063/1.1744102>`_ J. W. Cahn and J. E. Hilliard, ‘Free Energy of a Nonuniform System. I. Interfacial Free Energy’, The Journal of Chemical Physics, vol. 28, no. 2, pp. 258–267, Feb. 1958, doi: 10.1063/1.1744102.


`[2] <https://doi.org/10.48550/arXiv.1911.06718>`_ A. Lovrić, W. G. Dettmer, and D. Perić, ‘Low Order Finite Element Methods for the Navier-Stokes-Cahn-Hilliard Equations’. arXiv, Nov. 15, 2019. doi: 10.48550/arXiv.1911.06718.

`[3] <https://doi.org/10.48550/arXiv.1104.1336>`_ H. Abels, H. Garcke, and G. Grün, ‘Thermodynamically Consistent, Frame Indifferent Diffuse Interface Models for Incompressible Two-Phase Flows with Different Densities’. arXiv, Apr. 07, 2011. doi: 10.48550/arXiv.1104.1336.




  
.. [#cahn1958] \J. W. Cahn and J. E. Hilliard, ‘Free Energy of a Nonuniform System. I. Interfacial Free Energy’, The Journal of Chemical Physics, vol. 28, no. 2, pp. 258–267, Feb. 1958, doi: `10.1063/1.1744102 <https://dx.doi.org/10.1063/1.1744102>`_\.

.. [#lovric2019] \A. Lovrić, W. G. Dettmer, and D. Perić, ‘Low Order Finite Element Methods for the Navier-Stokes-Cahn-Hilliard Equations’. arXiv, Nov. 15, 2019. doi: `10.48550/arXiv.1911.06718 <https://doi.org/10.48550/arXiv.1911.06718>`_\.

.. [#abels2011] \H. Abels, H. Garcke, and G. Grün, ‘Thermodynamically Consistent, Frame Indifferent Diffuse Interface Models for Incompressible Two-Phase Flows with Different Densities’. arXiv, Apr. 07, 2011. doi: `10.48550/arXiv.1104.1336 <https://doi.org/10.48550/arXiv.1104.1336>`_\.

