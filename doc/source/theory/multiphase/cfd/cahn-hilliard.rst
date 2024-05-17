================================
Cahn-Hilliard Method
================================

**Under construction**

  
The Cahn-Hilliard system of equation is a model used to describe the process of phase separation based on the principle of free energy minimization. The key idea is that the system evolves to a state where the free energy is minimized, which often leads to the formation of distinct phases or regions within the material. This competition between the tendency of the system to minimize its overall free energy and the energy cost associated with creating new interfaces between phases is at the heart of the Cahn-Hilliard equation. Let us introduce those concepts formally.

Let :math:`\Omega = \Omega_0 \cup \Omega_1` be the domain formed by two fluids, namely fluid :math:`0` and :math:`1`, with :math:`\Gamma` denoting their interface and :math:`\partial \Omega`, the remaining boundaries. Like in the VOF method (LINK), we define a scalar function :math:`\phi` as a phase indicator such that:

.. math::
  \phi =
  \begin{cases}
    1 \phantom{-} \quad \forall \mathbf{x} \in \Omega_0\\
    -1 \quad \forall \mathbf{x} \in \Omega_1
  \end{cases}
  

Let us introduce the free energy functional :math:`\mathcal{F}`:

.. math::
  \mathcal{F} = \int_{\Omega} \Psi \mathrm{d}\Omega
  
:math:`\Psi` corresponds to the free energy density. Its expression is as follows:

.. math::
  \Psi = \lambda\left(F(\phi) + \frac{1}{2}|\nabla \phi|^2\right)
  
It is decomposed in a bulk free energy :math:`F(\phi)` and an interface energy :math:`\frac{1}{2}|\nabla \phi|^2`. :math:`\lambda` is called the mixing energy though not dimensionally consistent to an energy, see (LINK TO UNITS SUBSECTION). The bulk free energy has a double-well form, its expression is:

.. math::
  F(\phi) = \frac{(1-\phi^2)^2}{4\epsilon^2}

with :math:`\epsilon` a parameter related to the thickness of the interface between the phases.
  
To formalize the idea that the system tries to lower its free energy, we introduce a new variable, :math:`\eta`, called the chemical potential. It is defined as the variational derivative of :math:`\mathcal{F}`:

.. math::
  \eta = \frac{\delta\mathcal{F}}{\delta\phi} = \lambda\left(-\frac{\phi(1-\phi^2)}{\epsilon^2} + \nabla^2\phi\right)
  
  
Then, the phases have to move to satisfy free energy minimization, by going from high chemical potential regions to low chemical potential regions. Let us introduce the flux of phase due to chemical potential differences, denoted by :math:`\mathbf{J}`:

.. math::
  \mathbf{J} = -M(\phi)\nabla\eta
   
With :math:`M(\phi)` the mobility function. Then, we can write the following continuity equation, which also takes into account an external velocity field :math:`\mathbf{u}`:

.. math::
  \frac{\partial \phi}{\partial t} + \mathbf{u}\cdot \nabla \phi = \nabla \cdot (M(\phi)\nabla \eta)
  
For a 1D case, we obtain the following equilibrium phase field by solving the chemical potential equation with :math:`\eta = 0`:

.. math::
  \phi(x) = -\tanh{\left(\frac{x}{\sqrt{2}\epsilon}\right)}
  
  
Lastly, we need to relate the surface tension coefficient with the mixing energy. Cahn and Hilliard 
`[1] <https://doi.org/10.1063/1.1730447>`_ define surface tension as the excess free energy of the system due to the interface, the second term in the expression of :math:`\psi`. If we consider a system in equilibrium, and suppose a planar interface, we can write:

.. math::
  \sigma = \int_{-\infty}^{+\infty}\lambda \left(\frac{\mathrm{d}\phi}{\mathrm{d}x}\right)^2 \mathrm{d}x
  
Substituting the solution previously found, we obtain the following expression for :math:`\lambda`:

.. math::
  \lambda = \frac{3\epsilon\sigma}{2\sqrt{2}}
  
For the problem to have a unique solution, we give the following no-flux boundary conditions on :math:`\Gamma` for the phase field and chemical potential:

.. math::
  \nabla \phi \cdot\mathbf{n} = 0
  
  \nabla \eta \cdot \mathbf{n} = 0

  
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
  
  \int_\Omega \beta\lambda\left(\frac{\phi(1-\phi^2)}{\epsilon^2} - \nabla^2\phi\right)\mathrm{d}\Omega = 0
  
After using the integration by part and Green-Ostrogradski's theorem:

.. math::
  \int_\Omega \alpha\left(\frac{\partial \phi}{\partial t} + \mathbf{a}\cdot \nabla \phi\mathrm{d}\Omega\right) -\int_\Omega M(\phi) \nabla\alpha \cdot\nabla\eta \mathrm{d}\Omega + \cancelto{\mathrm{no-flux}}{\int_{\Gamma} M(\phi)\alpha \nabla \eta \cdot \mathbf{n} \mathrm{d}\Gamma} = 0
  
  \int_\Omega \beta\lambda\left(\frac{\phi(1-\phi^2)}{\epsilon^2}\right)\mathrm{d}\Omega - \int_\Omega \nabla \beta \cdot \nabla\phi\mathrm{d}\Omega + \cancelto{\mathrm{no-flux}}{\int_{\Gamma} \alpha \nabla \phi \cdot \mathbf{n} \mathrm{d}\Gamma} = 0
  
Using Petrov-Galerkin method, the finite element formulation reads:

Find :math:`(\phi^h,\eta^h) \in \psi^h` such that:
  
.. math::  
  \begin{array}{rl}
  \displaystyle \int_\Omega \alpha^h\left(\frac{\partial \phi^h}{\partial t} + \mathbf{a} \cdot \nabla \phi^h \right) \mathrm{d}\Omega - \int_\Omega M(\phi^h) \nabla \alpha^h \cdot \nabla \eta^h \mathrm{d}\Omega &= 0 \\[1em]
  \displaystyle \int_\Omega \beta^h \lambda \left( \frac{\phi^h(1-(\phi^h)^2)}{\epsilon^2} \right) \mathrm{d}\Omega - \int_\Omega \nabla \beta^h \cdot \nabla \phi^h \mathrm{d}\Omega &= 0 
  \end{array}
  \quad \forall (\alpha^h, \beta^h) \in \xi^h
   
 




  

  

